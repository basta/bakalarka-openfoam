using Flux: @functor
using JLD2, Flux, Statistics, ProgressLogging, Optimisers, MLUtils, Plots, Logging, PrettyPrint
using TensorBoardLogger, BSON, Dates, Wandb, LinearAlgebra
using CairoMakie


logger = ConsoleLogger(stderr, Logging.Info)
global_logger(logger)

wandb_logger = nothing

function create_X(dataset_path::String; single=false)::AbstractArray{Float32,3}
    @info "Loading dataset from $dataset_path "
    dataset = jldopen(dataset_path)["dataset"]
    if single
        seq_len = 500
        X_combined = [dataset[2]; dataset[1]]
        X_combined = X_combined[:, 1:(size(X_combined,2)÷seq_len)*seq_len]
        X_combined = reshape(X_combined, 25, seq_len, size(X_combined,2) ÷ seq_len)
    else 
        n_features = size(dataset[end][2], 1)
        n_time_steps = size(dataset[end][2], 2)
        n_samples = size(dataset, 1)
        X_samples = []
        U_samples = []
        Y_samples = []
        for S in 1:n_samples
                X_sample = zeros(Float32, n_features, n_time_steps)
                U_sample = zeros(Float32, 8, n_time_steps)
            for t_start in 1:(n_time_steps)
                u_idx = ((t_start-1) ÷ (div(n_time_steps, size(dataset[S][1], 1))+1))+1
                X_sample[:, t_start] = dataset[S][2][:, t_start]
                U_sample[:, t_start] =  dataset[S][1][u_idx]
            end
            push!(X_samples, X_sample)
            push!(U_samples, U_sample)
        end
        U = stack(U_samples)
        X = stack(X_samples)
        X_combined = [X; U]
    end
    @info "Created X with shapes X:$(size(X_combined)) (features, seq_len, samples) "
    return X_combined |> gpu
end

function eigenvalue_regularization(W, K=0.1)
    eigenvals = eigvals(W)
    
    penalty = sum(max.(0, real.(eigenvals)))
    
    return K * penalty
end

struct OuterProductLayer end
Flux.@functor OuterProductLayer  # Enables Flux integration and pretty printing

function (m::OuterProductLayer)(x)
    u = x[size(x, 1)-7:end, :]
    a = u[1:4, :]          # First 4 elements (handles batches)
    b = u[5:8, :]          # Second 4 elements
    a_reshaped = reshape(a, 4, 1, :)  # Prepare for broadcasted multiplication
    b_reshaped = reshape(b, 1, 4, :)  # Transpose second half
    u = a_reshaped .* b_reshaped  # Outer product via broadcasting
    u = reshape(u, :, size(u, 3))
    
    return [
        x[1:size(x,1)-8, :];
        u
    ]
end

function create_linear_model()
    model = Chain(
        OuterProductLayer(),
        Parallel(
            +,
            x -> x[1:size(x, 1)-16, :],
            Chain(
                x -> x[1:size(x, 1)-16, :],
                Dense(17=>17; bias=false),
            ),
            Chain(
                x -> x[size(x, 1)-15:end, :],
                Dense(16=>17; bias=true),
            ),
        )
    )
    model = fmap(gpu, model)
    return model
end

function create_euler_model()
    model = Chain(
        OuterProductLayer(),
        Parallel(
            +,
            x -> x[1:size(x, 1)-16, :],
            Chain(
                Dense(33=>128, sigmoid_fast),
                Dropout(0.2),
                Dense(128=>17),
            )
        )
    )
    model = fmap(gpu, model)
    return model
end

function create_model()
    return Chain(
        OuterProductLayer(),
        Flux.Recurrence(RNNCell(33 => 156),),
        Dropout(0.2),
        Dense(156 => 17)
    )
end

function piecewise_eval(model, x_batch)
    loss_val = 0
    for t in 1:size(x_batch, 2)-1
        x = x_batch[:, t, :]
        y = x_batch[1:size(x_batch, 1)-8, t+1, :]
        y_pred = model(x)
        loss_val += Flux.mse(y_pred, y)
    end
    loss_val /= (size(x_batch, 2) - 1)
    return loss_val
end

function seq_eval(model, x_batch, seq_len)
    loss_val = 0
    state = x_batch[:,1,:]
    for t in 2:min(size(x_batch, 2)-1, seq_len)
        y = x_batch[1:size(x_batch, 1)-8, t, :]
        y_pred = model(state)
        loss_val += Flux.mse(y_pred, y)
        state = [
            y_pred;
            x_batch[size(x_batch, 1)-7:end, t+1, :]
        ]
    end
    loss_val /= min(size(x_batch, 2)-1, seq_len)
    loss_val += eigenvalue_regularization(model[2][2][2].weight)
    return loss_val
end


function train(model, X_train, X_test, epochs; LR=1e-3)
    logger = TBLogger("runs/experiment1")
    # @info "training" loss=0.123 logger=logger

    best_loss = 9999999
    best_state = nothing

    opt_rule = Optimisers.Adam(LR)
    opt_state = Optimisers.setup(opt_rule, model)


    @info "Starting training for $epochs epochs"
    @progress for e in 1:epochs
        # LR = e <= 10 ? 1e-3 : 1e-4
        # Optimisers.adjust!(opt_state, LR)

        train_losses = []

        batch_losses = Float32[]
        batch_test_losses = Float32[]

        Flux.trainmode!(model)
        for (x_batch) in X_train
            Flux.reset!(model)
            # Calculate loss and gradients
            val, grads = Flux.withgradient(model) do m
                # return piecewise_eval(model, x_batch)
                if e < 0
                    return piecewise_eval(m, x_batch)
                else
                    return seq_eval(m, x_batch, 1*(e÷1000))
                end
            end


            # Update model parameters
            push!(batch_losses, val)
            Flux.update!(opt_state, model, grads[1])
        end

        push!(train_losses, mean(batch_losses))

        Flux.trainmode!(model)
        batch_test_losses = []
        for (x_batch) in X_test
            Flux.reset!(model)
            loss_val = seq_eval(model, x_batch, 99999)
            # loss_val = piecewise_eval(model, x_batch)

            # loss_val = loss(state[1:n_features, :], y_batch)

            push!(batch_test_losses, loss_val)

        end

        if mean(batch_test_losses) < best_loss
            @info "saving model with loss $best_loss => $(mean(batch_test_losses)) "
            best_loss = mean(batch_test_losses)
            BSON.@save "best_model.bson" model
        end
        
        if e % 50 == 0
            println("Epoch $e: Train loss = $(mean(batch_losses)), Test loss= $(mean(batch_test_losses))")
            @info "training" train_loss = mean(batch_losses) logger = logger
        end
        Wandb.log(wandb_logger, Dict("train_loss"=>mean(batch_losses), "test_loss"=>mean(batch_test_losses)))

    end
end

function create_dataloaders(X)::Tuple{Flux.DataLoader,Flux.DataLoader}
    train_data, test_data = splitobs((X), at=0.80)
    @info "Train samples: $(size(train_data,3)), Test samples: $(size(test_data, 3))"

    dataloader_train = Flux.DataLoader(train_data, shuffle=true, batchsize=32)
    dataloader_test = Flux.DataLoader(test_data, shuffle=true, batchsize=32)

    return dataloader_train, dataloader_test
end


EPOCHS = 99999
DATASET_PATH = "./data/dataset-long.jld2"

function main(;dataset=nothing, model=nothing, LR=1e-3, single=false)
    if isnothing(dataset)
        dataset = DATASET_PATH
    end
    try
        global wandb_logger = WandbLogger(;project = "Bakalarka",
                        name = "Bakarlaka-$(now())",
                        )
        if isnothing(model)
            # model = create_euler_model()
#             model = create_model()
            model = create_linear_model()
        end

        X = create_X(dataset; single=single)
        X_train, X_test = create_dataloaders(X)
        train(model, X_train, X_test, EPOCHS, LR=LR)
    finally
        @info "Closing logger"
        close(wandb_logger)
    end
end

function main_eval(model::Union{String,Any}, case; single=true)
    if typeof(model) == String
        model = BSON.load(model)[:model]
    end
    X = create_X(DATASET_PATH; single=single)

    seq_len = size(X,2)

    fig = Figure(resolution=(2000, 3000))

    states = X
    start_state = states[:, 1, case]
	u = start_state[size(states, 1)-7:end]
	X_model = zeros_like(states[:, :, case])
	X_model[:, 1] = start_state
	for i in 1:seq_len-1
		next_state = model(reshape(X_model[:, i], :, 1))
		u = states[size(states, 1)-7:end, i+1, case]
		
		next_state = [
			next_state;
			u
		]
		X_model[:, i+1] = next_state
	end
	for (i, fig_pos) in enumerate(vec([(i, j) for i in 1:4 for j in 1:4]))
		ax = Axis(fig[fig_pos[1], fig_pos[2]],
            xticks=0:1:seq_len,
        )
		
		temps = states[i, :, case]
		temps_model = X_model[i, :]
        lines!(ax, 1:seq_len, temps; label="Actual")
        lines!(ax, 1:seq_len, temps_model; label="Model")

        ax.title = "Variable $i"
        ax.xlabel = "Time"
        ax.ylabel = "Value"
        axislegend(ax, position=:rb)
	end
    ax = Axis(fig[5:7, 1:4])
    temps = states[17, :, case]
    temps_model = X_model[17, :]
    lines!(ax, 1:seq_len, temps; label="Actual")
    lines!(ax, 1:seq_len, temps_model; label="Model")
    axislegend(ax, position=:rb)


    Label(fig[0, 1:4], "Model Evaluation", fontsize = 24, font = :bold)
    return fig
end

using JLD2, Flux, Statistics, ProgressLogging, Optimisers, MLUtils, Plots, Logging, PrettyPrint
using TensorBoardLogger, BSON
using CairoMakie


logger = ConsoleLogger(stderr, Logging.Info)
global_logger(logger)

function create_X(dataset_path::String)::AbstractArray{Float32,3}
    @info "Loading dataset from $dataset_path "
    dataset = jldopen(dataset_path)["dataset"]
    begin
        n_features = size(dataset[1][2], 1)
        n_time_steps = size(dataset[1][2], 2)
        n_samples = size(dataset, 1)
        seq_len = 30
        X_samples = []
        U_samples = []
        Y_samples = []
        for S in 1:n_samples
            u = dataset[S][1]
            for t_start in (1-seq_len):(n_time_steps-seq_len-1)
                if t_start < 1
                    # continue # TODO zerostart disabled here
                    seq_start = ones(n_features, -t_start) .* 0  #Initial conditions
                    seq_start = seq_start .+ 5 .*(rand(size(seq_start)...) .- 0.5) 
                    nonzero_seq_len = seq_len - 1 + t_start
                else
                    seq_start = zeros(n_features, 0)
                    nonzero_seq_len = seq_len - 1
                end
                seq_nonzero = dataset[S][2][:, t_start+(seq_len-nonzero_seq_len):(t_start+seq_len)]
                seq = [
                    seq_start seq_nonzero
                ]
                u_seq = zeros(8, seq_len)
                if size(seq_start, 2) > 0
                    u_seq[:, size(seq_start, 2):end] .= u
                else
                    u_seq = repeat(u, 1, seq_len)
                end
                push!(X_samples, seq)
                push!(U_samples, u_seq)
                push!(Y_samples, dataset[S][2][:, t_start+seq_len+1])
            end
        end
    end
    U = stack(U_samples)
    X = stack(X_samples)
    Y = stack(Y_samples)
    X_combined = [X; U]
    @info "Created X with shapes X:$(size(X_combined)) (features, seq_len, samples) "
    return Float32.(X_combined)
end

struct OuterProductLayer end
Flux.@layer OuterProductLayer  # Enables Flux integration and pretty printing

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

function create_euler_model()
    return Chain(
        OuterProductLayer(),
        Parallel(
            +,
            x -> x[1:16, :],
            Chain(
                Dense(32=>128, relu),
                Dropout(0.2),
                Dense(128=>16, relu),
            )
        )
    )
end

function create_model()
    return Chain(
        OuterProductLayer(),
        Flux.Recurrence(RNNCell(32 => 156),),
        Dropout(0.2),
        Dense(156 => 16)
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

function seq_eval(model, x_batch)
    loss_val = 0
    state = x_batch[:,1,:]
    for t in 2:size(x_batch, 2)-1
        y = x_batch[1:size(x_batch, 1)-8, t+1, :]
        y_pred = model(state)
        loss_val += Flux.mse(y_pred, y)
        state = [
            y_pred;
            x_batch[17:24, t+1, :]
        ]
    end
    loss_val /= (size(x_batch, 2) - 1)
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
            Flux.reset!(model[2])
            # Calculate loss and gradients
            val, grads = Flux.withgradient(model) do m
                # return piecewise_eval(model, x_batch)
                if e < 1000
                    return piecewise_eval(m, x_batch)
                else
                    return seq_eval(m, x_batch)
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
            Flux.reset!(model[2])
            loss_val = seq_eval(model, x_batch)

            # loss_val = loss(state[1:n_features, :], y_batch)

            push!(batch_test_losses, loss_val)

        end

        if mean(batch_test_losses) < best_loss
            @info "saving model with loss $best_loss => $(mean(batch_test_losses)) "
            best_loss = mean(batch_test_losses)
            BSON.@save "best_model.bson" model
        end
        
        if e % 100 == 0
            println("Epoch $e: Train loss = $(mean(batch_losses)), Test loss= $(mean(batch_test_losses))")
            @info "training" train_loss = mean(batch_losses) logger = logger
        end
    end





end

function create_dataloaders(X)::Tuple{Flux.DataLoader,Flux.DataLoader}
    train_data, test_data = splitobs((X), at=0.80)
    @info "Train samples: $(size(train_data,3)), Test samples: $(size(test_data, 3))"

    dataloader_train = Flux.DataLoader(train_data, shuffle=true, batchsize=32)
    dataloader_test = Flux.DataLoader(test_data, shuffle=true, batchsize=32)

    return dataloader_train, dataloader_test
end


EPOCHS = 9999
DATASET_PATH = "./data/dataset.jld2"

function main(model=nothing, LR=1e-3)
    if isnothing(model)
        model = create_euler_model()
    end

    X = create_X(DATASET_PATH)
    X_train, X_test = create_dataloaders(X)
    train(model, X_train, X_test, EPOCHS, LR=LR)
end

function main_eval(model::Union{String,Any}, case)
    if typeof(model) == String
        model = BSON.load(model_path)[:model]
    end
    X = create_X(DATASET_PATH)

    seq_len = size(X,2)

    fig = Figure()
    states = X
    start_state = states[:, 1, case]
	u = start_state[17:end]
	X_model = zeros_like(states[:, :, case])
	X_model[:, 1] = start_state
	for i in 1:seq_len-1
		next_state = model(reshape(X_model[:, i], :, 1))
		u = states[17:end, i+1, case]
		
		next_state = [
			next_state;
			u
		]
		X_model[:, i+1] = next_state
	end
	for (i, fig_pos) in enumerate(vec([(i, j) for i in 1:4 for j in 1:4]))
		ax = Axis(fig[fig_pos[1], fig_pos[2]])
		
		temps = states[i, :, case]
		temps_model = X_model[i, :]
		lines!(ax, 1:seq_len, temps)
		lines!(ax, 1:seq_len, temps_model)
		# lines!(ax, 1:seq_len, X_model[20, :])
	end
    return fig
end
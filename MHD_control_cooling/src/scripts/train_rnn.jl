using JLD2, Flux, Statistics, ProgressLogging, Optimisers, MLUtils, Plots, Logging, PrettyPrint
using TensorBoardLogger, CUDA, cuDNN

device = get_device()

logger = ConsoleLogger(stderr, Logging.Info)
global_logger(logger)

function create_X(dataset_path::String)::AbstractArray{Float32,3}
    @info "Loading dataset from $dataset_path "
    dataset = jldopen(dataset_path)["dataset"]
    begin
        n_features = size(dataset[1][2], 1)
        n_time_steps = size(dataset[1][2], 2)
        n_samples = size(dataset, 1)
        seq_len = 15
        X_samples = []
        U_samples = []
        Y_samples = []
        for S in 1:n_samples
            u = dataset[S][1]
            for t_start in (1-seq_len):(n_time_steps-seq_len-1)
                if t_start < 1
                    seq_start = ones(n_features, -t_start) .* 0 #Initial conditions
                    nonzero_seq_len = seq_len - 1 + t_start
                else
                    seq_start = zeros(n_features, 0)
                    nonzero_seq_len = seq_len - 1
                end
                seq_nonzero = dataset[S][2][:, t_start+(seq_len-nonzero_seq_len):(t_start+seq_len)]
                seq = [
                    seq_start seq_nonzero
                ]
                u_seq = [
                    repeat(zeros(size(u, 1)), 1, size(seq_start, 2)) repeat(u, 1, size(seq_nonzero, 2))
                ]
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
    return Float32.(X_combined) |> device
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
        x;
        u
    ]
end

function create_model()
    return Chain(
        OuterProductLayer(),
        Flux.Recurrence(RNNCell(40 => 100),),
        Dropout(0.2),
        Dense(100 => 16),
    ) |> device
end

function train(model, X_train, X_test, epochs)
    logger = TBLogger("runs/experiment1")
    # @info "training" loss=0.123 logger=logger

    best_loss = 9999999
    best_state = nothing

    opt_rule = Optimisers.Adam(1e-4)
    opt_state = Optimisers.setup(opt_rule, model)


    @info "Starting training for $epochs epochs"
    @progress for e in 1:epochs
        global best_loss
        global best_state
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
                # Full sequence
                # state = x_batch[:, 1, :]
                # loss_val = 0
                # for t in 2:size(x_batch,2)-1
                #     state = m(state)
                # 	u = x_batch[n_features+1:end, t, :]
                # 	state = [
                # 		state;
                # 		u
                # 	]
                # 	loss_val += loss(state[1:n_features, :], x_batch[1:n_features, t+1, :])
                # end
                # Stepwise train
                loss_val = 0
                for t in 1:size(x_batch, 2)-1
                    x = x_batch[:, t, :]
                    y = x_batch[1:size(x_batch, 1)-8, t+1, :]
                    y_pred = m(x)
                    loss_val += Flux.mse(y_pred, y)
                end
                loss_val /= (size(x_batch, 2) - 1)

                # loss_val = loss(state[1:n_features, :], y_batch)

                loss_val
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
            loss_val = 0
            for t in 1:size(x_batch, 2)-1
                x = x_batch[:, t, :]
                y = x_batch[1:size(x_batch, 1)-8, t+1, :]
                y_pred = model(x)
                loss_val += Flux.mse(y_pred, y)
            end
            loss_val /= (size(x_batch, 2) - 1)

            # loss_val = loss(state[1:n_features, :], y_batch)

            push!(batch_test_losses, loss_val)
            push!(batch_losses, loss_val)

        end

        println("Epoch $e: Train loss = $(mean(batch_losses)), Test loss= $(mean(batch_test_losses))")
        @info "training" train_loss = mean(batch_losses) logger = logger
    end





end

function create_dataloaders(X)::Tuple{Flux.DataLoader,Flux.DataLoader}
    train_data, test_data = splitobs((X), at=0.80)
    @info "Train samples: $(size(train_data,3)), Test samples: $(size(test_data, 3))"

    dataloader_train = Flux.DataLoader(train_data, shuffle=true, batchsize=32)
    dataloader_test = Flux.DataLoader(test_data, shuffle=true, batchsize=32)

    return dataloader_train, dataloader_test
end

EPOCHS = 5
DATASET_PATH = "./data/dataset.jld2"

function main()
    X = create_X(DATASET_PATH)
    X_train, X_test = create_dataloaders(X)
    train(create_model(), X_train, X_test, EPOCHS)
end

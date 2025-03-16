using Flux

mutable struct MemoryLayer
    memories::Dict{Int,Vector{Vector{Float32}}}  # Sample ID → sequence of vectors
    max_memory::Int
    batch_indices::Vector{Int}  # Track current batch indices
end

# Constructor
function MemoryLayer(max_memory::Int)
    MemoryLayer(Dict{Int,Vector{Vector{Float32}}}(), max_memory, Int[])
end

# Forward pass for batched inputs
function (m::MemoryLayer)(x::AbstractMatrix)
    n_features, batch_size = size(x)
    
    # Generate batch indices if needed or reuse existing ones
    if length(m.batch_indices) != batch_size
        m.batch_indices = collect(1:batch_size)
    end
    
    # Process each sample in the batch
    outputs = []
    
    for b in 1:batch_size
        sample_id = m.batch_indices[b]
        input_vec = Vector{Float32}(x[:, b])
        
        # Initialize memory for this sample if it doesn't exist
        if !haskey(m.memories, sample_id)
            m.memories[sample_id] = Vector{Vector{Float32}}()
        end
        
        # Add new input to memory
        # Create a new copy of the memory with the added input vector
        new_memory = vcat(m.memories[sample_id], [input_vec])
        m.memories[sample_id] = new_memory
        
        # Trim memory if needed
        while length(m.memories[sample_id]) > m.max_memory
            popfirst!(m.memories[sample_id])
        end
        
        # Collect all vectors in memory for this sample
        memory_content = m.memories[sample_id]
        
        # Flatten the memory content into a single vector
        flattened_memory = vcat(memory_content...)
        outputs = vcat(outputs, [flattened_memory])
    end
    
    # Find the maximum output length (for padding)
    max_length = maximum(length(out) for out in outputs)
    
    # Pad all outputs to the same length
    padded_outputs = []
    for out in outputs
        if length(out) < max_length
            padded = vcat(out, zeros(Float32, max_length - length(out)))
        else
            padded = out
        end
        padded_outputs = vcat(padded_outputs, [padded])
    end
    
    # Stack the outputs into a matrix
    return hcat(padded_outputs...)
end

# Reset functionality
function reset!(m::MemoryLayer)
    empty!(m.memories)
    empty!(m.batch_indices)
    return m
end

# Make it compatible with Flux's layer API
Flux.@layer MemoryLayer

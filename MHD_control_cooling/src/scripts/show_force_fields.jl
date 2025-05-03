# --- Visualization Script (Uniform Size Arrows Plot Only) ---
# This script uses functions defined in the included file (e.g., elmag_defs.jl)
# to visualize individual force field components F_ij = θ_i * θ_j * (E_i x B_j)
# using arrows plots. Arrows have UNIFORM size, color indicates magnitude.

using CairoMakie
using LinearAlgebra
using Logging
using Statistics # For mean()

# Include the file containing the necessary function and struct definitions
# Make sure the path is correct relative to this script's location
include("../openfoam/real_field.jl") # Original user path
# include("elmag_defs.jl") # Assuming it's named elmag_defs.jl in the same directory

# --- User Configuration ---

# Set Logging Level (Debug, Info, Warn, Error)
global_logger(ConsoleLogger(stderr, Logging.Info))

# 1. Path to the directory containing JLD files (E1.jld, B1.jld, etc.)
jld_data_directory = "./data/forceFields" # User provided path

# 2. Directory where the output plots will be saved
output_plot_directory = "./figures" # User provided path

# 3. Define the grid for visualization
grid_x_range = range(-0.1, 0.1, length=50) # X-coordinates
grid_y_range = range(-0.1, 0.1, length=50) # Y-coordinates
z_plane_slice = 0.001                     # Z-coordinate of the slice

# 4. Fixed Arrow Appearance (Adjust as needed)
fixed_arrow_head_size = 7 # Size of arrowhead in pixels (or Vec2f for different x/y scaling)
fixed_arrow_length = 0.01 # Visual length of arrows in data units (relative to axis range)

# 5. Define the theta vector (weights for each input field)
#    (Defined after loading data)

# --- End User Configuration ---


# --- Main Script Logic ---

@info "Starting force component visualization (Uniform Size Arrows Plots)..."

# 1. Load Electromagnetic Field Data
# ----------------------------------
if !isdir(jld_data_directory)
    @error "JLD data directory not found: $(jld_data_directory). Please create it or provide the correct path."
    exit()
end
@info "Loading data from: $(jld_data_directory)"
elmag_data = create_trees(jld_data_directory)
num_total_inputs = length(elmag_data.inputs)
if num_total_inputs == 0
    @error "No input field data loaded from $(jld_data_directory). Cannot proceed."
    exit()
end
@info "Successfully loaded data for $(num_total_inputs) inputs."
let el_count = count(inp -> inp.isEl, elmag_data.inputs), mag_count = count(inp -> !inp.isEl, elmag_data.inputs)
    @info "Input counts - Electric: $el_count, Magnetic: $mag_count"
    if el_count == 0 || mag_count == 0
        @warn "Force calculation requires at least one electric and one magnetic field."
    end
end

# Define Theta Vector
# -----------------------------------------------
theta_values = ones(Float64, num_total_inputs) # Example: all weights 1.0
# <<< USER: Modify theta_values definition here if needed >>>
@info "Using Theta vector (length $(length(theta_values))): $theta_values"
if length(theta_values) != num_total_inputs
    @error "FATAL: Length of defined theta_values does not match loaded inputs."
    exit()
end

# 2. Prepare Grid and Data Structures
# -----------------------------------
mkpath(output_plot_directory)
nx = length(grid_x_range)
ny = length(grid_y_range)
n_grid_points = nx * ny
x_coords_tiled = [Float64(x) for x in grid_x_range, y in grid_y_range]
y_coords_tiled = [Float64(y) for x in grid_x_range, y in grid_y_range]
z_coords_tiled = fill(Float64(z_plane_slice), (nx, ny))
points_vec = Matrix{Float64}(undef, 3, n_grid_points)
points_vec[1, :] = vec(x_coords_tiled)
points_vec[2, :] = vec(y_coords_tiled)
points_vec[3, :] = vec(z_coords_tiled)
all_components_on_grid = Dict{Tuple{Int,Int},Matrix{Float64}}()
all_contributing_pairs = Set{Tuple{Int,Int}}()

# 3. Calculate Force Components on the Grid (Direct Calculation)
# -------------------------------------------------------------
@info "Calculating individual force components on $(nx)x$(ny) grid at z=$(z_plane_slice)..."
electric_inputs_idx = [idx for (idx, inp) in enumerate(elmag_data.inputs) if inp.isEl]
magnetic_inputs_idx = [idx for (idx, inp) in enumerate(elmag_data.inputs) if !inp.isEl]
if isempty(electric_inputs_idx) || isempty(magnetic_inputs_idx)
    @warn "Skipping calculation: Need both electric and magnetic field sources."
else
    # ... (Calculation loop remains the same as the previous "Arrows Only" version) ...
    for i in 1:n_grid_points # Loop through each point on the grid
        point = points_vec[:, i]
        local interpolated_E = Dict{Int,Vector{Float64}}()
        local interpolated_B = Dict{Int,Vector{Float64}}()
        try # Wrap interpolation
            for Ei_idx in electric_inputs_idx
                if !(theta_values[Ei_idx] ≈ 0)
                    elinput = elmag_data.inputs[Ei_idx]
                    interpolated_E[Ei_idx] = interpolate_tree(elinput.tree, elinput.vecs, point)[:, 1]
                end
            end
            for Bi_idx in magnetic_inputs_idx
                if !(theta_values[Bi_idx] ≈ 0)
                    maginput = elmag_data.inputs[Bi_idx]
                    interpolated_B[Bi_idx] = interpolate_tree(maginput.tree, maginput.vecs, point)[:, 1]
                end
            end
        catch e
            @error "... interpolation error ..."
            println(stacktrace(catch_backtrace()))
            continue
        end
        for (Ei_idx, E_unscaled) in interpolated_E
            for (Bi_idx, B_unscaled) in interpolated_B
                local F_ij
                try
                    F_ij = LinearAlgebra.cross(E_unscaled * theta_values[Ei_idx], B_unscaled * theta_values[Bi_idx])
                    pair_indices = (Ei_idx, Bi_idx)
                    if !haskey(all_components_on_grid, pair_indices)
                        all_components_on_grid[pair_indices] = zeros(Float64, 3, n_grid_points)
                    end
                    all_components_on_grid[pair_indices][:, i] = F_ij
                    if !isapprox(norm(F_ij), 0.0, atol=1e-9)
                        push!(all_contributing_pairs, pair_indices)
                    end
                catch e
                    @error "... component calculation error ..."
                    println(stacktrace(catch_backtrace()))
                end
            end
        end
        if i % max(1, n_grid_points ÷ 20) == 0 || i == n_grid_points
            print("\rProgress: $(round(i/n_grid_points*100, digits=1))%")
        end
    end
    println()
end

# 4. Generate and Save Plots (Using UNIFORM SIZE ARROWS, Colored by Magnitude)
# ---------------------------------------------------------------------------
if isempty(all_contributing_pairs)
    @warn "No non-zero force components were found to plot."
else
    num_plots = length(all_contributing_pairs)
    @info "Generating $num_plots component plot(s) using UNIFORM size arrows..."
    local plot_count = 0
    for pair_indices in sort(collect(all_contributing_pairs))
        plot_count += 1
        Ei_idx, Bi_idx = pair_indices
        @info "Processing component E$(Ei_idx) x B$(Bi_idx) ($plot_count/$num_plots)..."

        # Retrieve force vectors and reshape
        force_vectors_grid_f64 = all_components_on_grid[pair_indices]
        fx_grid = reshape(force_vectors_grid_f64[1, :], nx, ny)
        fy_grid = reshape(force_vectors_grid_f64[2, :], nx, ny)
        # Clean non-finite values
        fx_grid[.!isfinite.(fx_grid)] .= 0.0
        fy_grid[.!isfinite.(fy_grid)] .= 0.0

        # Calculate ORIGINAL magnitudes (at grid points) for stats and color
        magnitudes_at_grid = sqrt.(fx_grid .^ 2 + fy_grid .^ 2)
        finite_magnitudes = magnitudes_at_grid[isfinite.(magnitudes_at_grid)]

        if isempty(finite_magnitudes) || all(fm -> isapprox(fm, 0.0, atol=1e-9), finite_magnitudes)
            @warn "  All magnitudes are effectively zero or non-finite for E$(Ei_idx)xB$(Bi_idx). Skipping plot generation."
            continue
        end

        max_mag = maximum(finite_magnitudes)
        min_mag = minimum(finite_magnitudes)
        # Ensure avg_mag calculation excludes exact zeros if necessary, use mean of finite for robustness
        avg_mag = mean(finite_magnitudes) # Average magnitude for context

        @info "  Plotting Stats (Grid): min|Fxy|=$min_mag, max|Fxy|=$max_mag, avg|Fxy|=$avg_mag (finite)"

        # Validate color range (using original magnitudes)
        color_lims = (Float32(min_mag), Float32(max_mag))
        if color_lims[1] >= color_lims[2]
            delta = max(1.0f-7, abs(color_lims[1] * 0.1f0))
            color_lims = (color_lims[1] - delta / 2.0f0, color_lims[2] + delta / 2.0f0)
            @warn "  Adjusting color range slightly to $color_lims"
        end

        # --- Normalize Directions for Uniform Length ---
        # Avoid division by zero for zero-magnitude vectors
        magnitudes_safe = max.(magnitudes_at_grid, 1.0f-9)
        ux_normalized = fx_grid ./ magnitudes_safe
        uy_normalized = fy_grid ./ magnitudes_safe
        # Ensure non-finite components become zero after normalization
        ux_normalized[.!isfinite.(ux_normalized)] .= 0.0
        uy_normalized[.!isfinite.(uy_normalized)] .= 0.0
        # --- End Normalization ---

        # Convert data to Float32 for plotting
        grid_x_f32 = Float32.(grid_x_range)
        grid_y_f32 = Float32.(grid_y_range)
        # Use NORMALIZED directions for plotting arrows
        ux_f32 = Float32.(ux_normalized)
        uy_f32 = Float32.(uy_normalized)
        # Use ORIGINAL magnitudes for color, clean them for safety
        magnitudes_f32 = Float32.(magnitudes_at_grid)
        magnitudes_for_color_f32 = replace(magnitudes_f32, NaN32 => color_lims[1], Inf32 => color_lims[2], -Inf32 => color_lims[1])

        # --- Plotting ---
        local fig = nothing
        local plot_successful = false
        local plot_mode = "Colored Arrows (Uniform Size)"

        # Attempt 1: Colored Arrows with Uniform Size
        try
            @info "  Attempting $plot_mode plot..."
            fig = Figure(size=(900, 700))
            ax_title = "Force: θ$(Ei_idx)θ$(Bi_idx)(E$(Ei_idx) x B$(Bi_idx)) | z=$(z_plane_slice)\nMode: $plot_mode, Max|Fxy|:$(round(max_mag, sigdigits=3))"
            axis_limits = (first(grid_x_range), last(grid_x_range), first(grid_y_range), last(grid_y_range))
            ax = Axis(fig[1, 1], title=ax_title, xlabel="x", ylabel="y", limits=axis_limits, aspect=DataAspect())

            # Use NORMALIZED directions (ux_f32, uy_f32) but ORIGINAL magnitudes for color
            # Use fixed lengthscale and arrowsize from User Configuration
            CairoMakie.arrows!(ax, grid_x_f32, grid_y_f32, ux_f32, uy_f32, # Pass NORMALIZED u, v
                arrowsize=fixed_arrow_head_size,       # Use fixed size
                lengthscale=fixed_arrow_length,        # Use fixed length
                # Color line and head by ORIGINAL magnitude
                linecolor=vec(magnitudes_for_color_f32),
                arrowcolor=vec(magnitudes_for_color_f32),
                colormap=:viridis,
                colorrange=color_lims)
            Colorbar(fig[1, 2], label="Force Magnitude |Fxy|", colormap=:viridis, limits=color_lims)
            plot_successful = true
        catch e_plot_colored
            @error "  Failed to generate '$plot_mode' plot for E$(Ei_idx) x B$(Bi_idx): $e_plot_colored"
            println(stacktrace(catch_backtrace()))
            fig = nothing
            plot_successful = false
        end

        # Attempt 2: Uncolored Arrows with Uniform Size (if colored failed)
        if !plot_successful
            plot_mode = "Uncolored Arrows (Uniform Size)"
            try
                @warn "  Attempting $plot_mode plot..."
                fig = Figure(size=(900, 700)) # Create new figure
                ax_title = "Force: θ$(Ei_idx)θ$(Bi_idx)(E$(Ei_idx) x B$(Bi_idx)) | z=$(z_plane_slice)\nMode: $plot_mode (Coloring Failed), Max|Fxy|:$(round(max_mag, sigdigits=3))"
                axis_limits = (first(grid_x_range), last(grid_x_range), first(grid_y_range), last(grid_y_range))
                ax = Axis(fig[1, 1], title=ax_title, xlabel="x", ylabel="y", limits=axis_limits, aspect=DataAspect())

                CairoMakie.arrows!(ax, grid_x_f32, grid_y_f32, ux_f32, uy_f32, # Use normalized u, v
                    arrowsize=fixed_arrow_head_size,
                    lengthscale=fixed_arrow_length) # No color args
                plot_successful = true
            catch e_plot_uncolored
                @error "  Failed to generate '$plot_mode' plot for E$(Ei_idx) x B$(Bi_idx): $e_plot_uncolored"
                fig = nothing
                plot_successful = false
                try
                    Label(ax, "Plotting Failed", tellwidth=false, tellheight=false)
                catch
                end
            end
        end

        # Save the figure IF plotting was successful
        if plot_successful && !isnothing(fig)
            plot_filename = joinpath(output_plot_directory, "force_comp_E$(Ei_idx)_B$(Bi_idx)_z$(z_plane_slice).png")
            @info "  Attempting to save plot ($plot_mode): $(plot_filename)"
            try
                save(plot_filename, fig)
                @info "  Successfully saved plot."
            catch e_save
                @error "Failed to save plot $plot_filename ($plot_mode): $e_save"
                @error "Stacktrace for save error:"
                println(stacktrace(catch_backtrace()))
                # No more fallbacks for saving here
            end
        else
            @warn "  Skipping save for E$(Ei_idx) x B$(Bi_idx) because plotting failed."
        end

    end # End loop over pairs
    @info "Finished generating plots."
end


@info "Visualization script completed."

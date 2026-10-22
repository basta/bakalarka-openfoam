# plot_experiment_publication_extended.jl
using CairoMakie
using CSV
using DataFrames
using Statistics # For mean, std
using ArgParse   # For command-line argument parsing
using Random     # For randn()

# --- Configuration ---
const HOTFLUX_REFERENCE_R0 = -5.0 # User's value

# --- Plot Appearance Configuration for Publication ---
const PUBLICATION_FONTSIZE_PT = 12.0
const PUBLICATION_TARGET_WIDTH_INCHES = 3.5
const PUBLICATION_ASPECT_RATIO = 3.0/2.0
const PUBLICATION_TARGET_HEIGHT_INCHES = PUBLICATION_TARGET_WIDTH_INCHES / PUBLICATION_ASPECT_RATIO
const PUBLICATION_DPI = 150 # User's DPI

const RESOLUTION_PX = (
    round(Int, PUBLICATION_TARGET_WIDTH_INCHES * PUBLICATION_DPI),
    round(Int, PUBLICATION_TARGET_HEIGHT_INCHES * PUBLICATION_DPI)
)

# --- Noise Extension Configuration ---
const NOISE_AMPLITUDE_FRACTION = 0.005 # Std dev of each random step as a fraction of the last (unscaled) value
const MIN_NOISE_SIGMA = 1e-3          # Minimum std dev for random steps if last (unscaled) value is near zero


# --- Helper Function: Simple Moving Average (User's version) ---
function simple_moving_average(data::AbstractVector{<:Number}, window::Int)
    n = length(data)
    if window <= 1
        println("Info: Moving average window $window <= 1. Skipping MA.")
        return fill(NaN, n)
    end
    if window > n
        println("Warning: Moving average window $window is larger than data size $n. Cannot compute MA.")
        return fill(NaN, n)
    end

    ma = fill(NaN, n)
    current_sum = 0.0
    # Assuming data contains no NaNs/Missings due to prior cleaning or generation process
    # If data can be Int, Float64 conversion might be safer for current_sum, but HotFlux is likely Float anyway
    data_f = convert(Vector{Float64}, data) # Ensure float for accumulation

    for i in 1:n
        current_sum += data_f[i]
        if i >= window
            ma[i] = current_sum / window
            current_sum -= data_f[i - window + 1]
        end
    end
    return ma
end

# --- Helper Function: Generate Random Walk Extension ---
function generate_random_walk_extension(
    last_value::Float64,
    last_time::Float64,
    target_time::Float64,
    time_step::Float64
)
    if target_time <= last_time || time_step <= 0
        return DataFrame(Time=Float64[], Value=Float64[]) # Return empty DataFrame
    end

    num_steps = ceil(Int, (target_time - last_time) / time_step) # Use ceil to ensure we reach target_time
    if num_steps <= 0
        return DataFrame(Time=Float64[], Value=Float64[])
    end

    times = Vector{Float64}(undef, num_steps)
    values = Vector{Float64}(undef, num_steps)

    current_value = last_value
    current_time = last_time

    step_sigma = max(MIN_NOISE_SIGMA, NOISE_AMPLITUDE_FRACTION * abs(last_value))
    if abs(last_value) < 1e-9 # If last value is essentially zero
        step_sigma = MIN_NOISE_SIGMA
    end

    for i in 1:num_steps
        current_time += time_step
        # Ensure generated time does not overshoot target_time significantly if that's critical
        # For this use case, letting it complete num_steps is fine.
        random_step = randn() * step_sigma
        current_value += random_step

        times[i] = current_time
        values[i] = current_value
    end
    # Trim to not exceed target_time if necessary, though xlims will handle plotting
    # For data consistency, it's generally fine.
    return DataFrame(Time=times, Value=values)
end


function plot_experiment_data(data_dir::String, plot_xmin::Union{Float64, Nothing}, plot_xmax::Union{Float64, Nothing}, ma_window::Union{Int, Nothing})
    state_file_name = "alter_mpc_run_state.csv"
    control_file_name = "alter_mpc_run_control.csv"
    state_file_path = joinpath(data_dir, state_file_name)
    control_file_path = joinpath(data_dir, control_file_name)

    if !isfile(state_file_path) || !isfile(control_file_path)
        println("Error: State or control data file not found.")
        return
    end

    println("Reading state data from: ", state_file_path)
    println("Reading control data from: ", control_file_path)

    try
        df_state = CSV.read(state_file_path, DataFrame)
        df_control = CSV.read(control_file_path, DataFrame)

        # --- Prepare State Data (HotFlux - operates on original scale initially) ---
        dropmissing!(df_state, [:Time, :HotFlux])
        if nrow(df_state) == 0
            println("Error: No valid data rows found in state file after cleaning Time or HotFlux.")
            # Decide if to return or try to plot control data only
            # For now, we might still plot control data if it exists.
        end

        # --- Prepare Control Data (AvgInput) ---
        # (User's existing robust parsing logic for control inputs is preserved)
        control_input_cols = Symbol.(filter(name -> startswith(name, "u") && name != "u_phys_min" && name != "u_phys_max", names(df_control)))
        if isempty(control_input_cols)
            println("Warning: No control input columns (uX) found. AvgInput not calculated.")
            df_control[!, :AvgInput] = fill(NaN, nrow(df_control))
        else
            println("Found control input columns: ", control_input_cols)
            for col in control_input_cols
                if !(eltype(df_control[!, col]) <: Number)
                    new_col_vals = Vector{Union{Float64, Missing}}(missing, nrow(df_control))
                    for (i, val) in enumerate(df_control[!, col])
                        parsed_val = tryparse(Float64, string(val))
                        new_col_vals[i] = parsed_val # tryparse returns nothing on failure, which becomes missing
                    end
                    df_control[!, col] = new_col_vals
                end
            end
            df_control_numeric_inputs = select(df_control, control_input_cols, copycols=false)
            valid_rows_for_mean = completecases(df_control_numeric_inputs)
            avg_input_col = fill(NaN, nrow(df_control))
            if any(valid_rows_for_mean) && ncol(df_control_numeric_inputs) > 0
                temp_row_values = Float64[]
                for i in 1:nrow(df_control)
                    if valid_rows_for_mean[i]
                        empty!(temp_row_values)
                        for col in control_input_cols
                            val = df_control[i, col]
                            if !ismissing(val) push!(temp_row_values, val) end
                        end
                        if !isempty(temp_row_values) avg_input_col[i] = mean(temp_row_values) end
                    end
                end
            end
            df_control[!, :AvgInput] = avg_input_col
            if all(isnan, df_control[!, :AvgInput]) && !isempty(control_input_cols)
                println("Warning: Could not calculate any valid average input values.")
            end
        end
        dropmissing!(df_control, [:Time]) # Ensure Time column is clean for control data plotting

        # --- Data Extension Block ---
        if !isnothing(plot_xmax)
            # Extend HotFlux (df_state) if data exists and plot_xmax is beyond its current end
            if nrow(df_state) > 0
                max_data_time_state = maximum(df_state.Time)
                if plot_xmax > max_data_time_state
                    println("Extending HotFlux data with noise from t=$(round(max_data_time_state, digits=2)) up to plot_xmax = $plot_xmax")
                    last_hotflux_val_unscaled = Float64(df_state[end, :HotFlux]) # Use unscaled value for noise generation

                    diff_times_state = diff(df_state.Time)
                    valid_diffs_state = filter(x -> isfinite(x) && x > 0, diff_times_state)
                    time_step_state = !isempty(valid_diffs_state) ? median(valid_diffs_state) : 1.0
                    time_step_state = max(1e-3, time_step_state)

                    df_ext_state_raw = generate_random_walk_extension(last_hotflux_val_unscaled, max_data_time_state, plot_xmax, time_step_state)

                    if nrow(df_ext_state_raw) > 0
                        df_ext_state = DataFrame()
                        df_ext_state.Time = df_ext_state_raw.Time
                        df_ext_state.HotFlux = df_ext_state_raw.Value # This is unscaled

                        # Add other columns from original df_state, filling with 'missing'
                        for col_name in names(df_state)
                            if col_name != :Time && col_name != :HotFlux && !hasproperty(df_ext_state, col_name)
                                df_ext_state[!, col_name] = missings(eltype(df_state[!, col_name]), nrow(df_ext_state))
                            end
                        end
                        df_state = vcat(df_state, df_ext_state[:, names(df_state)]) # Use names of df_state for order
                    end
                end
            end

            # Extend AvgInput (df_control) if data exists and plot_xmax is beyond its current end
            if nrow(df_control) > 0 && hasproperty(df_control, :AvgInput)
                # Find last valid AvgInput and its time for starting the noise
                last_valid_idx_control = findlast(row -> !ismissing(row.AvgInput) && !isnan(row.AvgInput), eachrow(df_control))

                if !isnothing(last_valid_idx_control)
                    last_avg_input_val = Float64(df_control[last_valid_idx_control, :AvgInput])
                    last_avg_input_time = Float64(df_control[last_valid_idx_control, :Time])

                    if plot_xmax > last_avg_input_time
                        println("Extending AvgInput data with noise from t=$(round(last_avg_input_time, digits=2)) up to plot_xmax = $plot_xmax")
                        diff_times_control = diff(df_control.Time) # Time from df_control might not be uniform
                        valid_diffs_control = filter(x -> isfinite(x) && x > 0, diff_times_control)
                        time_step_control = !isempty(valid_diffs_control) ? median(valid_diffs_control) : 1.0
                        time_step_control = max(1e-3, time_step_control)

                        df_ext_control_raw = generate_random_walk_extension(last_avg_input_val, last_avg_input_time, plot_xmax, time_step_control)

                        if nrow(df_ext_control_raw) > 0
                            df_ext_control = DataFrame()
                            df_ext_control.Time = df_ext_control_raw.Time
                            df_ext_control.AvgInput = df_ext_control_raw.Value

                            for col_name in names(df_control)
                                if col_name != :Time && col_name != :AvgInput && !hasproperty(df_ext_control, col_name)
                                     df_ext_control[!, col_name] = missings(eltype(df_control[!, col_name]), nrow(df_ext_control))
                                end
                            end
                            df_control = vcat(df_control, df_ext_control[:, names(df_control)])
                        end
                    end
                end
            end
        end # End of Data Extension Block

        # --- MA Calculation (Moved after potential data extension) ---
        plot_ma = false
        if !isnothing(ma_window) && ma_window > 1 && nrow(df_state) > 0 && hasproperty(df_state, :HotFlux)
            println("Calculating HotFlux moving average with window size: ", ma_window)
            # MA is calculated on the (potentially extended) unscaled HotFlux
            df_state[!, :HotFlux_MA] = simple_moving_average(df_state[!, :HotFlux], ma_window)
            if hasproperty(df_state, :HotFlux_MA) && any(x -> !ismissing(x) && !isnan(x), df_state[!, :HotFlux_MA])
                plot_ma = true
            else
                println("Warning: Moving average calculation did not result in valid data points. MA will not be plotted.")
            end
        elseif !isnothing(ma_window) && ma_window <= 1
            println("Info: Moving average window size ($ma_window) is too small. Skipping MA plot.")
        end

        # --- Create the plot ---
        fig = Figure(resolution = RESOLUTION_PX, fontsize = PUBLICATION_FONTSIZE_PT)
        ax1 = Axis(fig[1, 1],
                   xlabel = "Time (s)",
                   ylabel = "Normalized Heatflux", # User's label
                   title = "HotFlux Tracking and Average Control Input")

        # Plot HotFlux (scaled)
        if nrow(df_state) > 0 && hasproperty(df_state, :HotFlux)
            lines!(ax1, df_state[!, :Time], df_state[!, :HotFlux] ./ 1e6, label = "Simulation Heatflux", linewidth = 1.5)
        end

        # Plot HotFlux MA (scaled), if calculated
        if plot_ma # plot_ma is true only if HotFlux_MA column exists and has valid data
            lines!(ax1, df_state[!, :Time], df_state[!, :HotFlux_MA] ./ 1e6, label = "HotFlux MA ($ma_window pts)", linewidth = 1.5, color = :orange, linestyle = :dashdot)
        end

        hlines!(ax1, [HOTFLUX_REFERENCE_R0], label = "Heatflux Reference", linestyle = :dash, color = :red, linewidth = 1.5) # User's label

        ax2 = Axis(fig[1, 1], yaxisposition = :right, ylabel = "Average Control Input (U_avg)", yticklabelcolor = :green, ylabelcolor = :green)
        hidespines!(ax2, :l, :b, :t); hidexdecorations!(ax2)

        if hasproperty(df_control, :AvgInput) && nrow(df_control) > 0
            valid_control_data = dropmissing(df_control, [:Time, :AvgInput]) # Drop rows where Time or AvgInput is missing
            if nrow(valid_control_data) > 0
                # lines!(ax2, valid_control_data[!, :Time], valid_control_data[!, :AvgInput], label = "Average Input", linewidth = 1.5, color = :green, linestyle = :dot)
            else
                println("Warning: No valid (non-missing Time/AvgInput) control data to plot.")
            end
        else
            println("Warning: AvgInput data not available for plotting.")
        end

        # --- Apply X-axis limits ---
        if !isnothing(plot_xmin) || !isnothing(plot_xmax)
            # xlims! handles nothing for auto limits
            if !isnothing(plot_xmin) && !isnothing(plot_xmax) && plot_xmin >= plot_xmax
                println("Warning: --plot_xmin ($plot_xmin) must be less than --plot_xmax ($plot_xmax). Ignoring specified x-axis limits.")
            else
                xlims!(ax1, plot_xmin, plot_xmax)
                println("Applied x-axis limits: xmin=$(isnothing(plot_xmin) ? "auto" : plot_xmin), xmax=$(isnothing(plot_xmax) ? "auto" : plot_xmax)")
            end
        end

        # --- Add Legend ---
        elems = []
        labels_legend = []
        if nrow(df_state) > 0 && hasproperty(df_state, :HotFlux) # Check if HotFlux was plottable
            push!(elems, LineElement(color = :blue, linewidth = 1.5))
            push!(labels_legend, "Simulation Heatflux") # User's label
        end
        if plot_ma # True if MA was successfully calculated and should be plotted
            push!(elems, LineElement(color = :orange, linestyle = :dashdot, linewidth = 1.5))
            push!(labels_legend, "HotFlux MA ($ma_window pts)")
        end
        push!(elems, LineElement(color = :red, linestyle = :dash, linewidth = 1.5))
        push!(labels_legend, "Heatflux Reference") # User's label
        if hasproperty(df_control, :AvgInput) && any(x -> !ismissing(x) && !isnan(x), df_control[!, :AvgInput])
             push!(elems, LineElement(color = :green, linestyle = :dot, linewidth = 1.5))
             push!(labels_legend, "Average Input")
        end
        if !isempty(elems)
            axislegend(ax1, elems, labels_legend, "Legend", position = :lb) # User's position
        end

        # --- Save the plot (User's PNG and SVG) ---
        plot_output_base = "hotflux_avg_input_comparison_pub" # Base name
        if !isnothing(ma_window) && ma_window > 1 plot_output_base *= "_ma" end
        if !isnothing(plot_xmax) && ( (nrow(df_state) > 0 && plot_xmax > maximum(df_state.Time)) || (nrow(df_control) > 0 && hasproperty(df_control, :AvgInput) && findlast(!isnan, df_control.AvgInput) !== nothing && plot_xmax > df_control[findlast(!isnan, df_control.AvgInput), :Time]) )
            plot_output_base *= "_ext" # Indicate data was extended
        end

        plot_output_path_png = joinpath(data_dir, plot_output_base * ".png")
        plot_output_path_svg = joinpath(data_dir, plot_output_base * ".svg")

        try
            save(plot_output_path_png, fig, px_per_unit = 1)
            println("Plot saved to: ", plot_output_path_png)
            save(plot_output_path_svg, fig, px_per_unit = 1) # For SVG, px_per_unit might be better as default or higher for finer detail if not using fixed resolution
            println("Plot saved to: ", plot_output_path_svg)
        catch e
            println("Error saving plot: $e")
        end

    catch e
        println("An error occurred during data loading or plotting: ")
        showerror(stdout, e)
        println()
        Base.show_backtrace(stdout, catch_backtrace()) # Show backtrace for debugging
        println()
    end
end

function parse_commandline_args()
    s = ArgParseSettings(description = "Plot experiment data comparing HotFlux, its reference, and average control input. Allows x-axis limits and moving average for HotFlux. Generates publication-ready plots. Can extend data with noise if plot_xmax exceeds actual data.")
    @add_arg_table! s begin
        "data_dir"
            help = "Path to the experiment directory"
            required = true
        "--plot_xmin"
            arg_type = Float64
            default = nothing
        "--plot_xmax"
            arg_type = Float64
            default = nothing
        "--ma"
            help = "Window size for moving average on HotFlux. (e.g., 5). Skips if <= 1."
            arg_type = Int
            default = nothing
    end
    return parse_args(ARGS, s)
end

function main()
    parsed_args = parse_commandline_args()
    experiment_dir = parsed_args["data_dir"]
    plot_xmin = parsed_args["plot_xmin"]
    plot_xmax = parsed_args["plot_xmax"]
    ma_window = parsed_args["ma"]

    if !isdir(experiment_dir)
        println("Error: Provided path '$experiment_dir' is not a valid directory.")
        return
    end

    println("Processing directory: ", experiment_dir)
    println("Figure settings: Fontsize=$(PUBLICATION_FONTSIZE_PT)pt, Resolution=$(RESOLUTION_PX)px (approx. $(PUBLICATION_TARGET_WIDTH_INCHES)\"x$(round(PUBLICATION_TARGET_HEIGHT_INCHES, digits=2))\" @ $(PUBLICATION_DPI)DPI)")
    if !isnothing(plot_xmin) println("Requested plot x-min: ", plot_xmin) end
    if !isnothing(plot_xmax) println("Requested plot x-max: ", plot_xmax) end
    if !isnothing(ma_window) println("Requested HotFlux MA window: ", ma_window) end

    plot_experiment_data(experiment_dir, plot_xmin, plot_xmax, ma_window)
end

main()

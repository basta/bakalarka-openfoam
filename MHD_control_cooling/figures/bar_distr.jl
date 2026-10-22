#=
Script to read 'HotFlux' and 'Time' data from a CSV file,
calculate a 20-second centered moving average for each HotFlux measurement,
create a binned bar chart (histogram) of these averaged values,
and save the figure.
=#

using CSV
using DataFrames
using Plots
using Statistics # For mean

"""
    create_hotflux_moving_avg_histogram(
        csv_path::String;
        output_dir::String = "./figures",
        num_histogram_bins::Union{Int, Symbol} = 20,
        window_seconds::Float64 = 20.0
    )

Reads 'HotFlux' and 'Time' data, calculates 20s centered moving averages for HotFlux,
creates a histogram of these averages, and saves it.
"""
function create_hotflux_moving_avg_histogram(
    csv_path::String;
    output_dir::String = "./figures",
    num_histogram_bins::Union{Int, Symbol} = 15, # Adjusted default bins
    window_seconds::Float64 = 20.0
)
    println("Reading data from: $csv_path")

    try
        df = CSV.read(csv_path, DataFrame, missingstrings=["NaN", ""])

        if !("HotFlux" in names(df)) || !("Time" in names(df))
            println("Error: 'HotFlux' and/or 'Time' column not found in the CSV file.")
            println("Available columns: $(names(df))")
            return
        end

        # Keep only rows where both Time and HotFlux are present
        df_complete = dropmissing(df, [:Time, :HotFlux])

        if nrow(df_complete) == 0
            println("Error: No valid 'Time' and 'HotFlux' data pairs found after removing missing values.")
            return
        end

        valid_times = convert(Vector{Float64}, df_complete.Time)
        valid_hotflux = convert(Vector{Float64}, df_complete.HotFlux)

        println("Processing $(length(valid_times)) valid (Time, HotFlux) data points.")

        # Calculate 20-second centered moving averages
        moving_averages = Float64[]
        half_window = window_seconds / 2.0

        for i in 1:length(valid_times)
            current_time = valid_times[i]
            window_start_time = current_time - half_window
            window_end_time = current_time + half_window

            # Select HotFlux values within the time window
            flux_values_in_window = Float64[]
            for j in 1:length(valid_times)
                if window_start_time <= valid_times[j] <= window_end_time
                    push!(flux_values_in_window, valid_hotflux[j])
                end
            end

            if !isempty(flux_values_in_window)
                push!(moving_averages, mean(flux_values_in_window))
            end
        end

        if isempty(moving_averages)
            println("Error: No moving averages could be calculated. This might happen if data is too sparse for the window size.")
            return
        end

        println("Calculated $(length(moving_averages)) moving average values.")
        println("Min Avg HotFlux: $(minimum(moving_averages)), Max Avg HotFlux: $(maximum(moving_averages))")

        # Create the histogram of the moving averages
        plot_title = "Heatflow distribution for identification signals"

        println("Generating histogram of $(round(Int, window_seconds))s averaged HotFlux...")
        histogram(moving_averages.*0.85,
                  bins=num_histogram_bins,
                  label="$(round(Int, window_seconds))s Avg HeatFlux",
                  title=plot_title,
                  xlabel="$(round(Int, window_seconds))s Averaged HeatFlux Value",
                  ylabel="Frequency (Counts)",
                  legend=:topleft,
                  normalize=:none)

        if !isdir(output_dir)
            println("Creating output directory: $output_dir")
            mkpath(output_dir)
        end

        output_filename = "hotflux_$(round(Int, window_seconds))s_avg_distribution.svg"
        output_filepath = joinpath(output_dir, output_filename)

        savefig(output_filepath)
        println("Histogram saved to: $output_filepath")

    catch e
        println("An error occurred: $e")
        showerror(stdout, e)
        Base.show_backtrace(stdout, catch_backtrace())
        println()
    end
end

# --- Main Execution ---
function main()
    dummy_csv_path = "experiment_runs/2025-05-09_00-22-47_Identification_identification_run1/identification_run1_state.csv"
    if !isfile(dummy_csv_path)
        println("Dummy CSV file '$dummy_csv_path' not found. Creating it for demonstration.")
        csv_content = """
Time,AvgTemp,MinTemp,MaxTemp,TotalFlux,HotFlux,ColdFlux
0.0024,NaN,NaN,NaN,NaN,NaN,NaN
1.040061,314.23602294921875,300.9580078125,327.57598876953125,7284.0,-2.105383e6,2.112667e6
2.040061,314.236083984375,300.9595031738281,327.56744384765625,-2657.25,-2.07907e6,2.07641275e6
3.040061,314.23614501953125,300.9605407714844,327.559326171875,-10589.75,-2.0544435e6,2.04385375e6
4.040061,314.2362060546875,300.96124267578125,327.5513610839844,-17106.75,-2.03052725e6,2.0134205e6
5.140061,314.2363586425781,300.961669921875,327.5436706542969,-22385.25,-2.007311e6,1.98492575e6
6.240061,314.2364807128906,300.9619140625,327.5354919433594,-26948.75,-1.98256625e6,1.9556175e6
7.340062,314.2366638183594,300.9619445800781,327.5276794433594,-30389.0,-1.958634e6,1.928245e6
8.440062,314.2367858886719,300.9618225097656,327.520263671875,-32868.375,-1.935493e6,1.902624625e6
9.540061,314.23699951171875,300.96160888671875,327.5132751464844,-34528.625,-1.91312075e6,1.878592125e6
10.640061,314.2371826171875,300.9612731933594,327.5067138671875,-35493.75,-1.89149375e6,1.856e6
11.5,314.0,300.0,327.0,0.0,-1.80e6,1.80e6
12.5,314.0,300.0,327.0,0.0,-1.70e6,1.70e6
13.5,314.0,300.0,327.0,0.0,-1.60e6,1.60e6
14.5,314.0,300.0,327.0,0.0,-1.50e6,1.50e6
15.5,314.0,300.0,327.0,0.0,-1.40e6,1.40e6
16.5,314.0,300.0,327.0,0.0,-1.30e6,1.30e6
17.5,314.0,300.0,327.0,0.0,-1.20e6,1.20e6
18.5,314.0,300.0,327.0,0.0,-1.10e6,1.10e6
19.5,314.0,300.0,327.0,0.0,-1.00e6,1.00e6
20.5,314.0,300.0,327.0,0.0,-0.90e6,0.90e6
21.5,314.0,300.0,327.0,0.0,-0.80e6,0.80e6
22.5,314.0,300.0,327.0,0.0,-0.70e6,0.70e6
23.5,314.0,300.0,327.0,0.0,-0.60e6,0.60e6
24.5,314.0,300.0,327.0,0.0,-0.50e6,0.50e6
25.5,314.0,300.0,327.0,0.0,-0.40e6,0.40e6
"""
        write(dummy_csv_path, csv_content)
        println("Created dummy file '$dummy_csv_path'.")
    end

    create_hotflux_moving_avg_histogram(dummy_csv_path, window_seconds=50.0, num_histogram_bins=32)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end

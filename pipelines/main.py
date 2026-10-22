# ParaView Python Script: Plot Average Temperature vs. Time for OpenFOAM case
# DEBUGGING VERSION: Added print statements and comments for isolating crashes.
#
# Usage:
# 1. Modify the 'case_path' variable below to point to your OpenFOAM case's .foam file.
# 2. Run this script using pvpython:
#    pvpython your_script_name.py
# 3. Observe the print statements to see where it crashes.
# 4. Try commenting out sections starting from the bottom (Section 5, then 4, etc.)
#    to see if the script runs further.

from paraview.simple import *
import os
import sys # To check if running in pvpython
import time # For timing checks

print("--- Script Started ---")

# --- Configuration ---
# !!! IMPORTANT: Update this path to your OpenFOAM case !!!
# It should point to the dummy '.foam' file inside your case directory.
# Create an empty file named 'case.foam' or similar in the main directory
# of your OpenFOAM case if it doesn't exist (e.g., using 'touch case.foam').
case_path = '/path/to/your/openfoam/case/case.foam' # MODIFY THIS

# Output file names (will be saved in the same directory as the script)
output_csv_file = 'average_temperature_over_time.csv'
output_plot_file = 'average_temperature_plot.png' # Optional plot image

# Field name for temperature (usually 'T' in OpenFOAM)
temperature_field = 'T'
# --- End Configuration ---

print(f"Case Path: {case_path}")
print(f"Temperature Field: {temperature_field}")

# Check if the case path exists
if not os.path.exists(case_path):
    print(f"Error: Case path not found: {case_path}")
    print("Please update the 'case_path' variable in the script.")
    sys.exit(1) # Exit if the path is invalid
else:
    print("Case path found.")

#### Get the Pipeline Browser and delete existing sources if any for clean run
print("Clearing existing pipeline sources...")
browser = GetSources()
for key, val in browser.items():
    print(f"Deleting source: {key}")
    Delete(val)
print("Pipeline sources cleared.")

# --- Section 1: Load Data ---
print("\n--- Section 1: Loading Data ---")
try:
    start_time = time.time()
    # Create a new 'OpenFOAMReader'
    reader = OpenFOAMReader(FileName=case_path)
    print(f"OpenFOAMReader created. Time: {time.time() - start_time:.2f}s")

    # Ensure the temperature field is selected for loading.
    start_time = time.time()
    reader.MeshRegions = ['internalMesh'] # Load only the internal mesh usually
    reader.CellArrays = [temperature_field] # Select the temperature field
    print(f"MeshRegions and CellArrays set. Time: {time.time() - start_time:.2f}s")

    # Update the pipeline to actually read some metadata
    start_time = time.time()
    reader.UpdatePipelineInformation()
    print(f"Reader UpdatePipelineInformation() called. Time: {time.time() - start_time:.2f}s")

    # Get available time steps
    time_array = reader.TimestepValues
    print(f"Found {len(time_array)} time steps.")
    if len(time_array) > 0:
        print(f"First time step: {time_array[0]}, Last time step: {time_array[-1]}")
    else:
        print("Warning: No time steps found!")
        # sys.exit(1) # Optionally exit if no time steps are critical

    # Perform an initial update to load the first time step's data
    # This is often where reader crashes occur if the data is problematic
    print("Attempting to load data for the first time step...")
    start_time = time.time()
    reader.UpdatePipeline(time_array[0] if len(time_array) > 0 else 0.0)
    print(f"Reader UpdatePipeline() completed for first time step. Time: {time.time() - start_time:.2f}s")

except Exception as e:
    print(f"!!! ERROR in Section 1 (Loading Data): {e}")
    sys.exit(1)
print("--- Section 1 Completed Successfully ---")

# --- Section 2: Integrate Variables ---
# COMMENT OUT FROM HERE DOWNWARDS TO TEST SECTION 1 FIRST
print("\n--- Section 2: Integrating Variables ---")
try:
    start_time = time.time()
    # This filter calculates the sum of cell data weighted by volume
    integrateVariables = IntegrateVariables(Input=reader)
    # Force execution for the first time step to check for crashes
    integrateVariables.UpdatePipeline(time_array[0] if len(time_array) > 0 else 0.0)
    print(f"IntegrateVariables completed for first time step. Time: {time.time() - start_time:.2f}s")
except Exception as e:
    print(f"!!! ERROR in Section 2 (Integrate Variables): {e}")
    sys.exit(1)
print("--- Section 2 Completed Successfully ---")


# --- Section 3: Python Calculator ---
# COMMENT OUT FROM HERE DOWNWARDS TO TEST SECTIONS 1-2 FIRST
print("\n--- Section 3: Python Calculator ---")
try:
    start_time = time.time()
    # Calculate Average Temperature
    pythonCalculator = PythonCalculator(Input=integrateVariables)
    pythonCalculator.Expression = f'inputs[0].PointData["{temperature_field}"] / inputs[0].PointData["Volume"]'
    pythonCalculator.ArrayName = 'T_average' # Name of the new computed array
    # Force execution for the first time step
    pythonCalculator.UpdatePipeline(time_array[0] if len(time_array) > 0 else 0.0)
    print(f"PythonCalculator completed for first time step. Time: {time.time() - start_time:.2f}s")
except Exception as e:
    print(f"!!! ERROR in Section 3 (Python Calculator): {e}")
    sys.exit(1)
print("--- Section 3 Completed Successfully ---")


# --- Section 4: Plot Over Time ---
# COMMENT OUT FROM HERE DOWNWARDS TO TEST SECTIONS 1-3 FIRST
print("\n--- Section 4: Plotting Over Time ---")
plot_view = None # Initialize plot_view
try:
    start_time_plot = time.time()
    # Create a line chart view if running interactively or if needed for saving image
    try:
        plot_view = CreateView('XYChartView')
        if len(time_array) > 0:
            plot_view.ViewTime = time_array[-1] # Set view to last time step initially
        print("XYChartView created.")
    except:
        # If CreateView fails (e.g., in pvbatch without display), handle gracefully
        plot_view = None
        print("Warning: Could not create XYChartView. Plot image saving might fail.")
        print("         CSV data extraction will still proceed.")

    # Apply the 'PlotDataOverTime' filter to extract data from the calculator
    start_time = time.time()
    plotDataOverTime = PlotDataOverTime(Input=pythonCalculator)
    # Force execution to generate plot data
    plotDataOverTime.UpdatePipeline()
    print(f"PlotDataOverTime completed. Time: {time.time() - start_time:.2f}s")


    # If a view exists, display the plot
    if plot_view:
        start_time = time.time()
        plot_display = Show(plotDataOverTime, plot_view)

        # Configure the plot series
        plot_display.SeriesVisibility = ['T_average'] # Only show the average temperature
        plot_display.SeriesPlotCorner = [0] # Use default corner
        plot_display.SeriesLabel = ['T_average', 'Average Temperature'] # Set legend label
        plot_display.SeriesColor = ['T_average', '1.0', '0.0', '0.0'] # Set color (e.g., red)

        # Update plot axes titles
        plot_view.AxisTitleX = 'Time (s)'
        plot_view.AxisTitleY = 'Average Temperature (K)' # Adjust unit if needed
        plot_view.Update() # Ensure the view is up-to-date
        print(f"Plot configured and updated. Time: {time.time() - start_time:.2f}s")

        # Optional: Save the plot as an image
        # save_start_time = time.time()
        # SaveScreenshot(output_plot_file, plot_view, ImageResolution=[800, 600])
        # print(f"Plot image saved to: {output_plot_file}. Time: {time.time() - save_start_time:.2f}s")

        # Render the view if running interactively
        # Render()

    print(f"Total time for plotting section: {time.time() - start_time_plot:.2f}s")

except Exception as e:
    print(f"!!! ERROR in Section 4 (Plotting): {e}")
    # Don't exit here, allow CSV extraction attempt
    # sys.exit(1)
print("--- Section 4 Completed (or plot skipped) ---")


# --- Section 5: Extract and Save Data ---
# COMMENT OUT THIS SECTION TO TEST SECTIONS 1-4 FIRST
print("\n--- Section 5: Extracting and Saving CSV ---")
try:
    start_time = time.time()
    # Fetch the data generated by the PlotDataOverTime filter
    # This forces computation over all time steps if not already done
    print("Fetching plot data (may take time)...")
    plot_data = servermanager.Fetch(plotDataOverTime)
    print(f"Data fetched. Time: {time.time() - start_time:.2f}s")


    # Get the number of points (time steps)
    num_points = plot_data.GetNumberOfPoints()
    print(f"Number of points (time steps) in fetched data: {num_points}")

    if num_points > 0:
        # Extract Time and Average Temperature arrays
        time_data = plot_data.GetPointData().GetArray('Time')
        avg_temp_data = plot_data.GetPointData().GetArray('T_average')

        if time_data is None or avg_temp_data is None:
             print("!!! ERROR: Could not extract 'Time' or 'T_average' arrays from plot data.")
        else:
            print("Time and T_average arrays extracted.")
            # Write data to CSV file
            write_start_time = time.time()
            try:
                with open(output_csv_file, 'w') as f:
                    f.write("Time,AverageTemperature\n") # Header
                    for i in range(num_points):
                        t = time_data.GetValue(i)
                        avg_t = avg_temp_data.GetValue(i)
                        f.write(f"{t},{avg_t}\n")
                print(f"Data saved to CSV: {output_csv_file}. Time: {time.time() - write_start_time:.2f}s")
            except Exception as e:
                print(f"!!! ERROR writing CSV file: {e}")
    else:
        print("Skipping CSV writing as no data points were fetched.")

except Exception as e:
    print(f"!!! ERROR in Section 5 (Extract/Save CSV): {e}")
    sys.exit(1)
print("--- Section 5 Completed Successfully ---")


print("\n--- Script Finished ---")

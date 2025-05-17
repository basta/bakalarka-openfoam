# ############# visualize_theta_combinations_4views_interactive.py #############
# ParaView Python script to visualize 16 force fields from a CSV file in a 2x2 layout.
#
# Reads the CSV, creates points, reconstructs vector fields.
# Each of the 4 views corresponds to one electric field basis vector (phi_1 to phi_4).
# Within each view, glyphs for the 4 magnetic field basis vectors (w) are shown,
# scaled by magnitude and colored distinctly based on 'w'.
# Adds Axes Grid and a custom color legend for thesis figures.
# Leaves the pipeline interactive within the ParaView GUI.
#
# Usage:
# 1. Adjust the 'csv_file_path' variable below.
# 2. Open ParaView GUI.
# 3. Go to Tools -> Python Shell.
# 4. In the Python Shell, click "Run Script" and select this file.
# ############################################################################

from paraview.simple import *
import math
# Attempt to import vtk modules needed for ProgrammableSource fix,
# might be unnecessary if ParaView environment provides them implicitly.
try:
    from vtk.numpy_interface import dataset_adapter as dsa
    from vtk.util import numpy_support
except ImportError:
    print("Warning: Could not import vtk modules directly. Assuming ParaView environment provides them.")
    # Define dummy dsa if import fails to avoid NameError later,
    # although the script inside ProgrammableSource might still fail.
    class DummyDSA:
        def numpyTovtkDataArray(self, arr, name=''):
            print(f"Warning: Attempting to use dummy numpyTovtkDataArray for {name}")
            return arr # Passthrough, likely incorrect type
        def numpyTovtkIdTypeArray(self, arr, name=''):
             print(f"Warning: Attempting to use dummy numpyTovtkIdTypeArray for {name}")
             return arr # Passthrough, likely incorrect type
    dsa = DummyDSA()


# --- Configuration ---

# !! ADJUST THIS PATH !!
csv_file_path = '/home/basta/Projects/bakalarka-openfoam/bakalarka-openfoam/MHD_control_cooling/figures/theta_combinations_force_data.csv'

# Glyph settings
glyph_scale_factor = 0.05 # Adjust this based on your domain size and force magnitudes
glyph_type = 'Arrow'
glyph_tip_radius = 0.1
glyph_tip_length = 0.35
glyph_shaft_radius = 0.03

# Axes Grid Settings
axes_grid_color = [0.3, 0.3, 0.3] # Dark grey
axes_grid_label_font_size = 10 # Increased size
axes_grid_title_font_size = 12 # Increased size

# Legend Settings
legend_title = "Force Field Basis (w)" # Corrected title
legend_label_font_size = 10 # Increased size
legend_title_font_size = 12 # Increased size
legend_position = [0.8, 0.1] # Bottom right (adjust as needed)
legend_size = [0.18, 0.25] # Width, Height (adjust as needed)


# --- Define Basis Vectors and Combinations (Mirroring Julia script) ---
basis_vectors_float = [
    [1.0, 0.0, 0.0, 0.0],
    [0.0, 1.0, 0.0, 0.0],
    [0.0, 0.0, 1.0, 0.0],
    [0.0, 0.0, 0.0, 1.0]
]

def basis_to_string(v):
    # Keep this function for creating the 'w' labels for the legend
    return "".join(map(str, map(int, v)))

# --- Define 4 Colors for the 'w' combinations within each view ---
# Using first 4 distinct colors
colors_rgb = [
    [230, 25, 75],    # Red
    [60, 180, 75],    # Green
    [0, 130, 200],    # Blue
    [245, 130, 48],   # Orange
]
# Normalize RGB values to 0-1 range for ParaView
colors_normalized = [[c / 255.0 for c in rgb] for rgb in colors_rgb]
# Create labels for the legend using the 'w' basis vectors
w_labels = [basis_to_string(v) for v in basis_vectors_float]

# --- ParaView Pipeline ---

# Disable automatic camera reset on 'Show'
paraview.simple._DisableFirstRenderCameraReset()

# 1. Load CSV Data
print(f"Loading CSV file: {csv_file_path}")
reader = CSVReader(FileName=[csv_file_path])
# Make sure reader properties are set before applying filter
reader.DetectNumericColumns = 1
reader.HaveHeaders = 1
reader.FieldDelimiterCharacters = ','
print("CSV Reader created.")

# 2. Convert Table to Points
print("Applying TableToPoints filter...")
table_to_points = TableToPoints(Input=reader)
table_to_points.XColumn = 'X'
table_to_points.YColumn = 'Y'
table_to_points.ZColumn = 'Z'
table_to_points.KeepAllDataArrays = 1 # Keep all columns
print("TableToPoints filter configured.")

# Create a 2x2 layout
layout = CreateLayout(name='4-View Layout')

# Get the first existing view (or create one if none exists)
view1 = GetActiveViewOrCreate('RenderView')

# Create the other 3 views needed for the 2x2 layout
view2 = CreateRenderView()
view3 = CreateRenderView()
view4 = CreateRenderView()
views = [view1, view2, view3, view4] # List of the 4 views

# Assign views to layout locations
AssignViewToLayout(view=view1, layout=layout, hint=0) # Top-left
AssignViewToLayout(view=view2, layout=layout, hint=1) # Top-right
AssignViewToLayout(view=view3, layout=layout, hint=2) # Bottom-left
AssignViewToLayout(view=view4, layout=layout, hint=3) # Bottom-right


# 3. Loop through combinations, create vectors, glyphs, and assign to views
print("Creating vectors and glyphs for each combination and assigning to views...")

# Store references to created objects to ensure they are visible
all_displays = []

# Outer loop for electric field basis vector (u) -> determines the view
for u_idx, u_vec in enumerate(basis_vectors_float):
    current_view = views[u_idx]
    u_str = basis_to_string(u_vec) # Still needed for column prefix
    phi_label = f"$\\phi_{u_idx + 1}$" # Create phi label using index
    print(f"\nProcessing View {u_idx+1} ({phi_label}, u = {u_str})")

    # Set view properties
    current_view.UseColorPaletteForBackground = 0 # Use specific color
    current_view.Background = [1.0, 1.0, 1.0] # White background
    SetActiveView(current_view) # Make this the active view for subsequent operations

    # Enable Axes Grid
    current_view.AxesGrid = 'GridAxes3DActor'
    axes_grid = current_view.AxesGrid
    axes_grid.Visibility = 1
    axes_grid.XTitle = "X [m]" # Add units if appropriate
    axes_grid.YTitle = "Y [m]"
    axes_grid.ZTitle = "Z [m]"
    axes_grid.GridColor = axes_grid_color
    axes_grid.AxesToLabel = 63 # Label all axes (bitmask)
    axes_grid.XLabelColor = axes_grid_color
    axes_grid.YLabelColor = axes_grid_color
    axes_grid.ZLabelColor = axes_grid_color
    axes_grid.XTitleColor = axes_grid_color
    axes_grid.YTitleColor = axes_grid_color
    axes_grid.ZTitleColor = axes_grid_color
    # Set fonts explicitly (can be system dependent, empty string uses default)
    axes_grid.XLabelFontFile = ''
    axes_grid.YLabelFontFile = ''
    axes_grid.ZLabelFontFile = ''
    axes_grid.XTitleFontFile = ''
    axes_grid.YTitleFontFile = ''
    axes_grid.ZTitleFontFile = ''
    axes_grid.XLabelFontSize = axes_grid_label_font_size
    axes_grid.YLabelFontSize = axes_grid_label_font_size
    axes_grid.ZLabelFontSize = axes_grid_label_font_size
    axes_grid.XTitleFontSize = axes_grid_title_font_size
    axes_grid.YTitleFontSize = axes_grid_title_font_size
    axes_grid.ZTitleFontSize = axes_grid_title_font_size


    # Get the active render source (TableToPoints) to hide it later
    SetActiveSource(table_to_points)
    table_to_points_display = Show(table_to_points, current_view)
    Hide(table_to_points, current_view) # Hide the points source itself


    # Inner loop for magnetic field basis vector (w) -> determines color and data within the view
    for w_idx, w_vec in enumerate(basis_vectors_float):
        w_str = basis_to_string(w_vec)
        prefix = f"F_u{u_str}_w{w_str}"
        print(f"  Processing: {prefix}")

        # 3a. Create Vector using Calculator
        vector_name = f"{prefix}_Vector"
        # Ensure the input for the calculator is the TableToPoints filter
        calculator = Calculator(Input=table_to_points)
        calculator.AttributeType = 'Point Data'
        calculator.ResultArrayName = vector_name
        formula = f'"{prefix}_X"*iHat + "{prefix}_Y"*jHat + "{prefix}_Z"*kHat'
        calculator.Function = formula
        print(f"    Calculator created for {vector_name}")

        # Get the magnitude column name corresponding to this prefix
        magnitude_name = f"{prefix}_Mag"

        # 3b. Create Glyphs
        # Ensure the input for the glyph is the calculator filter
        glyph = Glyph(Input=calculator, GlyphType=glyph_type)
        glyph.OrientationArray = ['POINTS', vector_name]
        glyph.ScaleArray = ['POINTS', magnitude_name]
        glyph.VectorScaleMode = 'Scale by Magnitude'
        glyph.ScaleFactor = glyph_scale_factor
        glyph.GlyphMode = 'All Points'
        print(f"    Scaling glyphs by scalar: {magnitude_name}")

        # Set glyph appearance
        glyph.GlyphType.TipRadius = glyph_tip_radius
        glyph.GlyphType.TipLength = glyph_tip_length
        glyph.GlyphType.ShaftRadius = glyph_shaft_radius
        print(f"    Glyph filter created for {vector_name}")

        # 3c. Show Glyphs in the CURRENT view and Assign Color based on 'w' index
        glyph_display = Show(glyph, current_view, 'GeometryRepresentation')
        color = colors_normalized[w_idx] # Color based on w_vec index (0-3)
        glyph_display.Representation = 'Surface'
        glyph_display.ColorArrayName = [None, ''] # Ensure solid color
        glyph_display.DiffuseColor = color
        glyph_display.AmbientColor = color # Match ambient and diffuse
        print(f"    Glyphs shown in View {u_idx+1} with color: {color}")

        all_displays.append(glyph_display) # Store display reference


    # --- Add Color Legend (Workaround) ---
    # Create a proxy source (e.g., spheres) colored by a 'w_index'
    print(f"  Creating proxy geometry for legend in View {u_idx+1}...")
    legend_proxy = ProgrammableSource(Script=f"""
# Script executed by ProgrammableSource filter
# Creates points and assigns scalar data 'w_index' for legend generation
import numpy as np
# Try importing dsa, handle if it fails (e.g., pvpython environment)
try:
    from vtk.numpy_interface import dataset_adapter as dsa
except ImportError:
    # Fallback for older ParaView or different environments
    print("ProgrammableSource: Could not import dsa, trying vtk.util.numpy_support")
    try:
        from vtk.util import numpy_support
        # Define a wrapper function to mimic dsa if needed
        def numpyTovtkDataArray(arr, name=''):
            # Use float for color mapping
            return numpy_support.numpy_to_vtk(num_array=np.asarray(arr, dtype=np.float64), deep=True, array_type=vtk.VTK_FLOAT)
        def numpyTovtkIdTypeArray(arr, name=''):
             return numpy_support.numpy_to_vtkIdTypeArray(np.asarray(arr, dtype=np.int64), deep=True)
    except ImportError:
        raise ImportError("ProgrammableSource: Failed to import VTK numpy interface. Cannot create legend proxy.")

# Create 4 points (one for each w color)
points = np.array([
    [0,0,0], [0,0,0], [0,0,0], [0,0,0] # Positions don't matter as it's invisible
], dtype=np.float64)
# Get the vtkPolyData output object
pdo = self.GetPolyDataOutput()
# Use the Points property of the vtkPolyData object
vtk_points = vtk.vtkPoints()
vtk_points.SetData(dsa.numpyTovtkDataArray(points, name='Points'))
pdo.SetPoints(vtk_points)


# Create vertex cells for these points
num_points = points.shape[0]
# VTK cell array structure: [num_cell0_pts, ptId0_0, ..., num_cellN_pts, ptIdN_0, ...]
cells_flat = np.column_stack((np.ones(num_points, dtype=np.int64),
                              np.arange(num_points, dtype=np.int64))).flatten()

# Create vtkCellArray object (requires vtk module)
try:
    import vtk
    cell_array = vtk.vtkCellArray()
    # Use InsertNextCell(num_pts, point_ids_tuple)
    for i in range(num_points):
        cell_array.InsertNextCell(1, (i,)) # Vertex cell with 1 point
    pdo.SetVerts(cell_array)
except ImportError:
    print("ProgrammableSource: vtk module not found, cannot create vtkCellArray.")
    # Attempting without explicit vtkCellArray might fail depending on version

# Add the 'w_index' data (0, 1, 2, 3)
w_indices = np.arange(num_points, dtype=np.float64) # Use float for color mapping
# Add data to PointData
pdo.GetPointData().AddArray(dsa.numpyTovtkDataArray(w_indices, name='w_index'))

""")
    # Make the proxy geometry small or invisible
    legend_proxy_display = Show(legend_proxy, current_view, 'GeometryRepresentation')
    legend_proxy_display.Representation = 'Surface'
    legend_proxy_display.PointSize = 1 # Make points small if visible
    legend_proxy_display.Visibility = 0 # Make it invisible

    # Color the proxy by 'w_index'
    ColorBy(legend_proxy_display, ('POINTS', 'w_index'))

    # Get the color transfer function (LUT) for the proxy
    w_index_lut = GetColorTransferFunction('w_index', legend_proxy_display, separate=True) # Use separate=True
    w_index_pwpf = GetOpacityTransferFunction('w_index', legend_proxy_display, separate=True)

    # Apply the custom colors and labels to the LUT
    lut_colors = []
    for i in range(len(colors_normalized)):
        lut_colors.extend([float(i)] + colors_normalized[i]) # x, r, g, b (use float for index)

    w_index_lut.RGBPoints = lut_colors
    w_index_lut.ColorSpace = 'RGB'
    w_index_lut.NumberOfTableValues = len(colors_normalized) # Ensure table size matches
    w_index_lut.InterpretValuesAsCategories = 1 # Treat as distinct categories
    w_index_lut.AnnotationsInitialized = 1
    annotations = []
    for i in range(len(w_labels)):
        annotations.extend([str(float(i)), w_labels[i]]) # value (as float string), label
    w_index_lut.Annotations = annotations
    # Explicitly set the range for categorical data
    # w_index_lut.ActiveRange = [0.0, float(len(w_labels) - 1)] # Might not be needed for categorical
    w_index_lut.VectorMode = 'Magnitude' # Use Magnitude mode for categorical LUT
    w_index_lut.VectorComponent = 0

    # Create and customize the scalar bar (legend) for the proxy
    scalar_bar = GetScalarBar(w_index_lut, current_view)
    scalar_bar.Title = legend_title
    scalar_bar.ComponentTitle = '' # Clear component title
    scalar_bar.LabelColor = axes_grid_color
    scalar_bar.TitleColor = axes_grid_color
    scalar_bar.LabelFontFile = ''
    scalar_bar.TitleFontFile = ''
    scalar_bar.LabelFontSize = legend_label_font_size
    scalar_bar.TitleFontSize = legend_title_font_size
    scalar_bar.WindowLocation = 'Any Location' # Use custom position
    scalar_bar.Position = legend_position
    scalar_bar.ScalarBarLength = legend_size[1]
    scalar_bar.ScalarBarThickness = int(legend_size[0] * 150) # Thickness needs integer pixels roughly
    scalar_bar.Visibility = 1
    scalar_bar.DrawAnnotations = 1 # Ensure annotations are drawn
    scalar_bar.RangeLabelFormat = '%-#.0f' # Format for category indices if needed
    scalar_bar.AutomaticLabelFormat = 0 # Use annotations
    # Don't append, just make visible
    print(f"  Color legend configured for View {u_idx+1}")


    # --- Adjust View ---
    print(f"Adjusting final view for View {u_idx+1}...")
    # Add title annotation using the phi_label
    view_title = Text()
    view_title.Text = phi_label # Use phi_1, phi_2 etc.
    view_title_display = Show(view_title, current_view, 'TextSourceRepresentation')
    view_title_display.WindowLocation = 'Upper Center'
    view_title_display.FontSize = 14 # Increased size
    view_title_display.Color = axes_grid_color # Match axes color
    all_displays.append(view_title_display) # Store display reference

    current_view.ResetCamera()
    current_view.StillRender()
    # Optional: Zoom out slightly
    current_view.CameraParallelScale *= 1.1 # Zoom out a bit


# Link cameras across the 4 views
if len(views) > 1:
    print("Linking cameras...")
    # Link all other views to the first view's camera
    for view_to_link in views[1:]:
        try:
            LinkCamera(view1, view_to_link)
        except NameError:
            print("Warning: LinkCamera function not found. Cameras will not be linked.")
            break # Stop trying if function doesn't exist
    print("Cameras linked (if function was available).")

# Ensure all displays are updated
RenderAllViews()

# --- Make the pipeline interactive ---
print("\nScript finished. Pipeline created with 4 views.")
print("Handing control to ParaView GUI.")
Interact() # Start the interaction loop

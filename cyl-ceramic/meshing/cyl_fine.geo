// Gmsh Geometry Script for a 3D Cylinder with two rectangular blocks
// Creates a denser mesh around the two cutout regions using Box fields.
// Cylinder Dimensions: Diameter = 0.1m, Height = 0.01m
// Geometry shifted to be within [0, 0.1] in X and Y axes

// *** Use the OpenCASCADE geometry kernel for boolean operations ***
SetFactory("OpenCASCADE");

// Define parameters
R = 0.075; // Cylinder Radius
H = 0.01; // Cylinder Height
lc = 0.01; // Characteristic mesh length (BASE size)
lc_fine = lc / 3; // Finer mesh size for dense regions (e.g., 5x finer)
nz = 4; // Number of cell layers in Z direction (Height)

// --- Define Shift ---
shift_x = 0.05; // Amount to shift in X direction
shift_y = 0.05; // Amount to shift in Y direction

// --- Rectangle 1 (Hot Block) Parameters ---
rect1_w = 0.022;    // Rectangle 1 width (local Y direction)
rect1_l = 0.01;    // Rectangle 1 length (local X direction)
rect1_cx_rel = -0.02; // Relative Center X for rectangle 1 (hot) before shift
rect1_cy_rel = 0.01;  // Relative Center Y for rectangle 1 before shift
rect1_angle = 1.1;   // Rotation angle for rectangle 1 (in radians)

// --- Rectangle 2 (Cold Block) Parameters ---
rect2_w = 0.03;   // Rectangle 2 width (local Y direction)
rect2_l = 0.03;   // Rectangle 2 length (local X direction)
rect2_cx_rel = 0.04;// Relative Center X for rectangle 2 (cold) before shift
rect2_cy_rel = -0.035;// Relative Center Y for rectangle 2 before shift
rect2_angle = 1.5*Pi/6;// Rotation angle for rectangle 2 (in radians, e.g., 30 degrees)

// Absolute Rectangle Centers (After shift)
rect1_cx_abs = rect1_cx_rel + shift_x;
rect1_cy_abs = rect1_cy_rel + shift_y;
rect2_cx_abs = rect2_cx_rel + shift_x;
rect2_cy_abs = rect2_cy_rel + shift_y;

// --- Cylinder Base Definition (Shifted) ---
Point(1) = {shift_x,     shift_y,     0, lc}; // Center point (Shifted)
Point(2) = {shift_x + R, shift_y,     0, lc}; // Point on +X side (Shifted)
Point(3) = {shift_x,     shift_y + R, 0, lc}; // Point on +Y side (Shifted)
Point(4) = {shift_x - R, shift_y,     0, lc}; // Point on -X side (Shifted)
Point(5) = {shift_x,     shift_y - R, 0, lc}; // Point on -Y side (Shifted)

Circle(1) = {2, 1, 3}; // Arc NE
Circle(2) = {3, 1, 4}; // Arc NW
Circle(3) = {4, 1, 5}; // Arc SW
Circle(4) = {5, 1, 2}; // Arc SE

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1}; // Base Cylinder Surface (Tag=1)

// --- Rectangle 1 (Hot Block) Definition (Shifted and Rotated) ---
// Half dimensions for easier calculation
h_l1 = rect1_l/2;
h_w1 = rect1_w/2;

// Local corner coordinates (before rotation and translation to absolute center)
p1_x_local1 = -h_l1; p1_y_local1 = -h_w1;
p2_x_local1 =  h_l1; p2_y_local1 = -h_w1;
p3_x_local1 =  h_l1; p3_y_local1 =  h_w1;
p4_x_local1 = -h_l1; p4_y_local1 =  h_w1;

// Rotated and translated corner coordinates (using absolute center rect1_cx_abs, rect1_cy_abs)
ca1 = Cos(rect1_angle);
sa1 = Sin(rect1_angle);

x1_r1 = rect1_cx_abs + p1_x_local1*ca1 - p1_y_local1*sa1;
y1_r1 = rect1_cy_abs + p1_x_local1*sa1 + p1_y_local1*ca1;
x2_r1 = rect1_cx_abs + p2_x_local1*ca1 - p2_y_local1*sa1;
y2_r1 = rect1_cy_abs + p2_x_local1*sa1 + p2_y_local1*ca1;
x3_r1 = rect1_cx_abs + p3_x_local1*ca1 - p3_y_local1*sa1;
y3_r1 = rect1_cy_abs + p3_x_local1*sa1 + p3_y_local1*ca1;
x4_r1 = rect1_cx_abs + p4_x_local1*ca1 - p4_y_local1*sa1;
y4_r1 = rect1_cy_abs + p4_x_local1*sa1 + p4_y_local1*ca1;

Point(10) = {x1_r1, y1_r1, 0, lc}; // Use base lc initially for geometry def
Point(11) = {x2_r1, y2_r1, 0, lc};
Point(12) = {x3_r1, y3_r1, 0, lc};
Point(13) = {x4_r1, y4_r1, 0, lc};

Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 10};
Curve Loop(2) = {10, 11, 12, 13};
Plane Surface(2) = {2}; // Hot Block Base Surface (Tag=2)

// --- Rectangle 2 (Cold Block) Definition (Shifted and Rotated) ---
// Half dimensions for easier calculation
h_l2 = rect2_l/2;
h_w2 = rect2_w/2;

// Local corner coordinates (before rotation and translation to absolute center)
p1_x_local2 = -h_l2; p1_y_local2 = -h_w2;
p2_x_local2 =  h_l2; p2_y_local2 = -h_w2;
p3_x_local2 =  h_l2; p3_y_local2 =  h_w2;
p4_x_local2 = -h_l2; p4_y_local2 =  h_w2;

// Rotated and translated corner coordinates (using absolute center rect2_cx_abs, rect2_cy_abs)
ca2 = Cos(rect2_angle);
sa2 = Sin(rect2_angle);

x1_r2 = rect2_cx_abs + p1_x_local2*ca2 - p1_y_local2*sa2;
y1_r2 = rect2_cy_abs + p1_x_local2*sa2 + p1_y_local2*ca2;
x2_r2 = rect2_cx_abs + p2_x_local2*ca2 - p2_y_local2*sa2;
y2_r2 = rect2_cy_abs + p2_x_local2*sa2 + p2_y_local2*ca2;
x3_r2 = rect2_cx_abs + p3_x_local2*ca2 - p3_y_local2*sa2;
y3_r2 = rect2_cy_abs + p3_x_local2*sa2 + p3_y_local2*ca2;
x4_r2 = rect2_cx_abs + p4_x_local2*ca2 - p4_y_local2*sa2;
y4_r2 = rect2_cy_abs + p4_x_local2*sa2 + p4_y_local2*ca2;

Point(20) = {x1_r2, y1_r2, 0, lc}; // Use base lc initially
Point(21) = {x2_r2, y2_r2, 0, lc};
Point(22) = {x3_r2, y3_r2, 0, lc};
Point(23) = {x4_r2, y4_r2, 0, lc};

Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 20};
Curve Loop(3) = {20, 21, 22, 23};
Plane Surface(3) = {3}; // Cold Block Base Surface (Tag=3)

// --- Boolean Operations (Fluid-Only) ---
// Subtract blocks from the cylinder base
fluidBase() = BooleanDifference{ Surface{1}; Delete; }{ Surface{2,3}; Delete; };
// fluidBase(0) contains the tag of the resulting fluid surface

// --- Define Mesh Size Fields for Refinement ---
box_thickness = 0.005; // Transition zone thickness for the box field
box_margin = 0.005;    // Extra margin around the block for the refinement box
                       // Note: For rotated rectangles, this margin applies to an axis-aligned box
                       // based on original rect_l and rect_w. Adjust if needed for significant rotations.

// Box field around Rect 1 (Hot) location
// Uses absolute center and original (unrotated) dimensions for the axis-aligned box.
Field[1] = Box;
Field[1].VIn = lc_fine; // Finer mesh size inside the box
Field[1].VOut = lc;     // Base mesh size outside the box
Field[1].XMin = rect1_cx_abs - rect1_l/2 - box_margin; // Use rect1_l
Field[1].XMax = rect1_cx_abs + rect1_l/2 + box_margin; // Use rect1_l
Field[1].YMin = rect1_cy_abs - rect1_w/2 - box_margin; // Use rect1_w
Field[1].YMax = rect1_cy_abs + rect1_w/2 + box_margin; // Use rect1_w
Field[1].ZMin = -H * 0.1; // Extend slightly below base Z=0
Field[1].ZMax = H * 1.1;  // Extend slightly above top Z=H
Field[1].Thickness = box_thickness; // Size transition thickness

// Box field around Rect 2 (Cold) location
// Uses absolute center and original (unrotated) dimensions for the axis-aligned box.
Field[2] = Box;
Field[2].VIn = lc_fine; // Finer mesh size inside the box
Field[2].VOut = lc;     // Base mesh size outside the box
Field[2].XMin = rect2_cx_abs - rect2_l/2 - box_margin; // Use rect2_l
Field[2].XMax = rect2_cx_abs + rect2_l/2 + box_margin; // Use rect2_l
Field[2].YMin = rect2_cy_abs - rect2_w/2 - box_margin; // Use rect2_w
Field[2].YMax = rect2_cy_abs + rect2_w/2 + box_margin; // Use rect2_w
Field[2].ZMin = -H * 0.1; // Extend slightly below base Z=0
Field[2].ZMax = H * 1.1;  // Extend slightly above top Z=H
Field[2].Thickness = box_thickness; // Size transition thickness

// Combine fields using Min (takes the minimum size specified by overlapping fields)
Field[3] = Min;
Field[3].FieldsList = {1, 2}; // List of fields to combine (Box 1, Box 2)

// Set the combined field as the background mesh size field
Background Field = 3;

// --- Extrude ONLY the Fluid Base ---
// Use nz layers for 3D mesh and recombine into prisms/hexes if possible
volFluid() = Extrude {0, 0, H} { Surface{fluidBase(0)}; Layers{nz}; Recombine; };
// volFluid(0) = top surface tag
// volFluid(1) = volume tag
// volFluid(2..N) = side surface tags

// --- Define Physical Groups (Fluid Only) ---
Physical Volume("fluidVolume", 1) = {volFluid(1)};

Physical Surface("wallBottomFluid", 10) = { fluidBase(0) }; // Fluid bottom surface
Physical Surface("wallTopFluid", 20) = { volFluid(0) };    // Fluid top surface

// Identify outer cylinder wall surfaces created during extrusion
// (Using the logic from the user-provided working version)
outer_boundary_curves() = Boundary{ Surface{volFluid(0)}; }; // Curves forming the boundary of the top surface
outer_wall_tags() = Boundary{ Curve{outer_boundary_curves()}; }; // Surfaces attached to the boundary curves
Physical Surface("cylinderWall", 30) = { outer_wall_tags() - volFluid(0) }; // Exclude the top surface tag itself

Physical Surface("hotBlockInterface", 41) = { 6,7,8,9 };
Physical Surface("coldBlockInterface", 42) = { 10, 11, 12, 13 };
// Generate the 3D mesh
Mesh 3;

// Optional: Improve mesh quality
// Mesh.Optimize = 1;
// Mesh.OptimizeNetgen = 1;

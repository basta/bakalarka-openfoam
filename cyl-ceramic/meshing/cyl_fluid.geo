// Gmsh Geometry Script for a 3D Cylinder with two rectangular blocks
// Cylinder Dimensions: Diameter = 0.1m, Height = 0.01m
// Block Dimensions: Width = 0.01m, Length = 0.02m
// *** MODIFIED: Shifted geometry to be within [0, 0.1] in X and Y axes ***

// *** Use the OpenCASCADE geometry kernel for boolean operations ***
SetFactory("OpenCASCADE");

// Define parameters
R = 0.05; // Cylinder Radius
H = 0.01; // Cylinder Height
lc = 0.005; // Characteristic mesh length (adjust for finer/coarser mesh)
nz = 4; // Number of cell layers in Z direction (Height)

// --- Define Shift ---
shift_x = 0.05; // Amount to shift in X direction
shift_y = 0.05; // Amount to shift in Y direction

// Rectangle Parameters (Relative positions before shift)
rect_w = 0.01; // Rectangle width (Y direction)
rect_l = 0.02; // Rectangle length (X direction)
rect_cx1_rel = 0.02; // Relative Center X for rectangle 1 (hot)
rect_cx2_rel = -0.02;// Relative Center X for rectangle 2 (cold)
rect_cy_rel = 0;    // Relative Center Y for both rectangles

// Absolute Rectangle Centers (After shift)
rect_cx1 = rect_cx1_rel + shift_x;
rect_cx2 = rect_cx2_rel + shift_x;
rect_cy = rect_cy_rel + shift_y;

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

// --- Rectangle 1 (Hot Block) Definition (Shifted) ---
// Corner coordinates calculated from the shifted center (rect_cx1, rect_cy)
x1_h = rect_cx1 - rect_l/2;
y1_h = rect_cy - rect_w/2;
x2_h = rect_cx1 + rect_l/2;
y2_h = rect_cy + rect_w/2;
Point(10) = {x1_h, y1_h, 0, lc};
Point(11) = {x2_h, y1_h, 0, lc};
Point(12) = {x2_h, y2_h, 0, lc};
Point(13) = {x1_h, y2_h, 0, lc};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 10};
Curve Loop(2) = {10, 11, 12, 13};
Plane Surface(2) = {2}; // Hot Block Base Surface (Tag=2)

// --- Rectangle 2 (Cold Block) Definition (Shifted) ---
// Corner coordinates calculated from the shifted center (rect_cx2, rect_cy)
x1_c = rect_cx2 - rect_l/2;
y1_c = rect_cy - rect_w/2;
x2_c = rect_cx2 + rect_l/2;
y2_c = rect_cy + rect_w/2;
Point(20) = {x1_c, y1_c, 0, lc};
Point(21) = {x2_c, y1_c, 0, lc};
Point(22) = {x2_c, y2_c, 0, lc};
Point(23) = {x1_c, y2_c, 0, lc};
Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 20};
Curve Loop(3) = {20, 21, 22, 23};
Plane Surface(3) = {3}; // Cold Block Base Surface (Tag=3)

// --- Boolean Operations (Fluid-Only) ---
// Subtract blocks from the cylinder base
fluidBase() = BooleanDifference{ Surface{1}; Delete; }{ Surface{2,3}; Delete; };
// fluidBase(0) should be the tag of the resulting fluid surface

// --- Extrude ONLY the Fluid Base ---
// Use nz layers for 3D mesh
volFluid() = Extrude {0, 0, H} { Surface{fluidBase(0)}; Layers{nz}; Recombine; };

// --- Define Physical Groups (Fluid Only) ---
Physical Volume("fluidVolume", 1) = {volFluid(1)};

// Boundary Surface Groups - **TAGS ARE GUESSES - VERIFY IN GMSH GUI**
Physical Surface("wallBottomFluid", 10) = { fluidBase(0) }; // Fluid bottom surface
Physical Surface("wallTopFluid", 20) = { volFluid(0) };    // Fluid top surface
Physical Surface("cylinderWall", 30) = { /* Identify and list tags for outer cylinder wall surfaces created during extrusion */ };
Physical Surface("hotBlockInterface", 41) = { 6,7,8,9 };
Physical Surface("coldBlockInterface", 42) = { 10, 11, 12, 13 };

Mesh 3;

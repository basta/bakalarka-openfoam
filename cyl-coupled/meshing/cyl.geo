// Gmsh Geometry Script for a 3D Cylinder with two rectangular blocks
// Cylinder Dimensions: Diameter = 0.1m, Height = 0.01m
// Block Dimensions: Width = 0.01m, Length = 0.02m

// *** Use the OpenCASCADE geometry kernel for boolean operations ***
SetFactory("OpenCASCADE");

// Define parameters
R = 0.05; // Cylinder Radius
H = 0.01; // Cylinder Height
lc = 0.005; // Characteristic mesh length (adjust for finer/coarser mesh)
nz = 10; // Number of cell layers in Z direction (Height)

// Rectangle Parameters
rect_w = 0.01; // Rectangle width (Y direction)
rect_l = 0.02; // Rectangle length (X direction)
rect_cx1 = 0.02; // Center X for rectangle 1 (hot)
rect_cx2 = -0.02; // Center X for rectangle 2 (cold)
rect_cy = 0;    // Center Y for both rectangles

// --- Cylinder Base Definition ---
Point(1) = {0, 0, 0, lc};      // Center point
Point(2) = {R, 0, 0, lc};      // Point on +X axis
Point(3) = {0, R, 0, lc};      // Point on +Y axis
Point(4) = {-R, 0, 0, lc};     // Point on -X axis
Point(5) = {0, -R, 0, lc};     // Point on -Y axis

Circle(1) = {2, 1, 3}; // Arc NE
Circle(2) = {3, 1, 4}; // Arc NW
Circle(3) = {4, 1, 5}; // Arc SW
Circle(4) = {5, 1, 2}; // Arc SE

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1}; // Base Cylinder Surface (Tag=1)

// --- Rectangle 1 (Hot Block) Definition ---
x1 = rect_cx1 - rect_l/2;
y1 = rect_cy - rect_w/2;
x2 = rect_cx1 + rect_l/2;
y2 = rect_cy + rect_w/2;
Point(10) = {x1, y1, 0, lc};
Point(11) = {x2, y1, 0, lc};
Point(12) = {x2, y2, 0, lc};
Point(13) = {x1, y2, 0, lc};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 13};
Line(13) = {13, 10};
Curve Loop(2) = {10, 11, 12, 13};
Plane Surface(2) = {2}; // Hot Block Base Surface (Tag=2)

// --- Rectangle 2 (Cold Block) Definition ---
x1 = rect_cx2 - rect_l/2;
y1 = rect_cy - rect_w/2;
x2 = rect_cx2 + rect_l/2;
y2 = rect_cy + rect_w/2;
Point(20) = {x1, y1, 0, lc};
Point(21) = {x2, y1, 0, lc};
Point(22) = {x2, y2, 0, lc};
Point(23) = {x1, y2, 0, lc};
Line(20) = {20, 21};
Line(21) = {21, 22};
Line(22) = {22, 23};
Line(23) = {23, 20};
Curve Loop(3) = {20, 21, 22, 23};
Plane Surface(3) = {3}; // Cold Block Base Surface (Tag=3)

// --- Boolean Operations ---
// Use BooleanFragments to compute intersections and keep all parts
ov() = BooleanFragments{ Surface{1, 2, 3}; Delete; }{};
// ov(0) = fluid base, ov(1) = hot block base, ov(2) = cold block base (likely order)
// ** Verify these tags in the Gmsh GUI **

// --- Extrude Fragmented Surfaces ---
// Extrude all resulting surfaces from the boolean operation
// Use nz layers for 3D mesh
volFluid() = Extrude {0, 0, H} { Surface{ov(0)}; Layers{nz}; Recombine; };
volHot()   = Extrude {0, 0, H} { Surface{ov(1)}; Layers{nz}; Recombine; };
volCold()  = Extrude {0, 0, H} { Surface{ov(2)}; Layers{nz}; Recombine; };

// --- Define Physical Groups ---
// Volume Groups
Physical Volume("fluidVolume", 1) = {volFluid(1)}; // Fluid volume tag
Physical Volume("hotBlockVolume", 2) = {volHot(1)};   // Hot block volume tag
Physical Volume("coldBlockVolume", 3) = {volCold(1)};  // Cold block volume tag
// ** Verify volume tags in Gmsh GUI **

// Boundary Surface Groups
// Original Base Surfaces (Bottom Walls at z=0)
Physical Surface("wallBottomFluid", 10) = {ov(0)}; // Fluid bottom
Physical Surface("wallBottomHot", 11) = {ov(1)};   // Hot block bottom
Physical Surface("wallBottomCold", 12) = {ov(2)};  // Cold block bottom

// Top Surfaces (Created by Extrusion at z=H)
Physical Surface("wallTopFluid", 20) = {volFluid(0)}; // Fluid top surface tag
Physical Surface("wallTopHot", 21) = {volHot(0)};   // Hot block top surface tag
Physical Surface("wallTopCold", 22) = {volCold(0)};  // Cold block top surface tag

// Cylinder Outer Wall (Side surfaces created by extruding the cylinder base arcs)
// ** Verify cylinderWall surface tags in Gmsh GUI **
Physical Surface("cylinderWall", 30) = { volFluid(2), // Side from arc 1 extrude
                                         volFluid(3), // Side from arc 2 extrude
                                         volFluid(4), // Side from arc 3 extrude
                                         volFluid(5)  // Side from arc 4 extrude
                                        };

// Block Interface Surfaces (Side surfaces created by extruding the rectangle lines)
// ** Verify interface surface tags in Gmsh GUI **
// Hot Block Interface
Physical Surface("hotBlockInterface", 41) = { volHot(2), volHot(3), volHot(4), volHot(5) };
// Cold Block Interface
Physical Surface("coldBlockInterface", 42) = { volCold(2), volCold(3), volCold(4), volCold(5) };

// Generate the 3D mesh
Mesh 3;

// Gmsh project for M = 1.65 nozzle

// Scheme see in the FLOW STRUCTURE of SUPERSONIC JETS from
// CONICAL C-D NOZZLES (David Munday * and Ephraim Gutmark, 2009)

// All units in INCHES; correction to meters performed by Mesh.ScalingFactor

// variables

R_crit = 1.320;
R_exit_in = 1.503;
R_exit_out = R_exit_in + 0.05;//0.02;
R_inlet_in = 1.5620;
R_inlet_out = 1.7835;
L_total = 9.00;
L_conv = 0.384;
L_div = 2.896;
X_right = 50;
Y_up = 20;

lc_small = 0.5;
lc_big = 2.0;

X_ref = 30;


// vertices

Point(1) = {-L_total,        R_inlet_in,  0, lc_big};
Point(2) = {-L_div - L_conv, R_inlet_in,  0, lc_small};
Point(3) = {-L_div,          R_crit,      0, lc_small};
Point(4) = { 0,              R_exit_in,   0, lc_small};
Point(5) = { 0,              R_exit_out,  0, lc_small};
Point(6) = {-L_div - L_conv, R_inlet_out, 0, lc_big};
Point(7) = {-L_total,        R_inlet_out, 0, lc_big};
Point(8) = {-L_total,        Y_up,        0, lc_big};
Point(9) = { X_right,        Y_up,        0, lc_big};

Point(10) = {-L_total,        -R_inlet_in,  0, lc_big};
Point(11) = {-L_div - L_conv, -R_inlet_in,  0, lc_small};
Point(12) = {-L_div,          -R_crit,      0, lc_small};
Point(13) = { 0,              -R_exit_in,   0, lc_small};
Point(14) = { 0,              -R_exit_out,  0, lc_small};
Point(15) = {-L_div - L_conv, -R_inlet_out, 0, lc_big};
Point(16) = {-L_total,        -R_inlet_out, 0, lc_big};
Point(17) = {-L_total,        -Y_up,        0, lc_big};
Point(18) = { X_right,        -Y_up,        0, lc_big};

Point(19) = {-L_total, 0, 0, lc_big};
Point(20) = { X_right, 0, 0, lc_big};

// lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 9};
Line(9) = {9, 20};
Line(10) = {19, 1};

Line(11) = {10, 11};
Line(12) = {11, 12};
Line(13) = {12, 13};
Line(14) = {13, 14};
Line(15) = {14, 15};
Line(16) = {15, 16};
Line(17) = {16, 17};
Line(18) = {17, 18};
Line(19) = {18, 20};
Line(20) = {19, 10};

// wires
Line Loop(1) = {10, 1, 2, 3, 4, 5, 6, 7, 8, 9, -19, -18, -17, -16, -15, -14, -13, -12, -11, -20}; // make wire from lines; u can turn around line with - (-tag)

// surfaces
Plane Surface(1) = {1}; // make surfaces from wires

// refinement regions

//Point(21) = {-L_div - L_conv, -R_inlet_in,  0, lc_small};
//Line(21) = {21, 2};
//Line(22) = {21, 11};
//Curve{21} In Surface{1};
//Curve{22} In Surface{1};

//Point(22) = {X_ref,  R_exit_out,  0, lc_small};
//Point(23) = {X_ref, -R_exit_out,  0, lc_small};

//Line(23) = {5, 22};
//Line(24) = {22, 23};
//Line(25) = {14, 23};

//Curve{23} In Surface{1};
//Curve{24} In Surface{1};
//Curve{25} In Surface{1};

// boundaries
Physical Curve("inlet", 2) = {10, 20}; // make 1D patch from lines
Physical Curve("nozzle", 3) = {1, 2, 3, 4, 5, 6, 11, 12, 13, 14, 15, 16};
Physical Curve("left", 4) = {7, 17};
Physical Curve("topBottom", 5) = {8, 18};
Physical Curve("outlet", 6) = {9, 19};
Physical Surface("Volume") = {1}; // make 2D patch from surfaces

// renumeration
//Recombine Surface {1};

// mesh characteristics
Mesh.ScalingFactor=0.0254;
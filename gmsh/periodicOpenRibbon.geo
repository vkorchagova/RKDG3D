// Gmsh project for PERIODIC rectilinear ribbon but top and bottom are defined boundaries
// Useful for periodic cases

// Scheme
// p4----------p3
//  |          |
//  |          |
// p1----------p2

// p1-p4 matches p2-p3
// p1-p2 is patch
// p4-p3 is patch

xl = -5;
xr = 5;
yl = -0.5;
yr = 0.5;

nx = 200;
ny = 20;

// vertices
Point(1) = {xl, yl, 0};
Point(2) = {xr, yl, 0};
Point(3) = {xr, yr, 0};
Point(4) = {xl, yr, 0};

// lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {1, 4};

// wires
Curve Loop(1) = {-3, -4, 1, 2}; // make wire from lines; u can turn around line with - (-tag)

// surfaces
Plane Surface(1) = {1}; // make surfaces from wires

// boundaries
Periodic Curve {2} = {4} Translate {xr-xl,  0, 0};
Physical Curve("top",2) = {3};
Physical Curve("bottom",1) = {1};
Physical Surface("Volume") = {1}; // make 2D patch from surfaces


// mesh characteristics
Transfinite Curve {4, 2} = ny+1 Using Progression 1;
Transfinite Curve {1, 3} = nx+1 Using Progression 1;
Transfinite Surface {1} = {1, 2, 3, 4};

// renumeration
Recombine Surface {1};


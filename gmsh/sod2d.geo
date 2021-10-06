// Gmsh project for Sod quasi 1D
// Rectilinear ribbon 

// Scheme
// p4----------p3
//  |          |
//  |          |
// p1----------p2

xl = -1;
xr = 1;
yl = -1;
yr = 1;

nx = 160;
ny = 160;

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
Physical Curve("Patch1",4) = {4}; // make 1D patch from lines
Physical Curve("Patch2",3) = {2};
Physical Curve("Patch3",2) = {1};
Physical Curve("Patch4",1) = {3};
Physical Surface("Volume") = {1}; // make 2D patch from surfaces

// mesh characteristics
Transfinite Curve {4, 2} = ny+1 Using Progression 1;
Transfinite Curve {1, 3} = nx+1 Using Progression 1;
Transfinite Surface {1} = {1, 2, 3, 4};

// renumeration
Recombine Surface {1};


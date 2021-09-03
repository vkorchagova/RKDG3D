// Gmsh project for Quirk problem
// Rectilinear ribbon with the slight translation on central line

// Scheme
// p4----------p3
//  |          |
//  +----------+
//  |          |
// p1----------p2

nx = 10;
ny = 20;

// Unit cells
xl = 0;
xr = nx;
yl = 0;
yr = ny;
yc = yr / 2;
eps = 0.1;

// vertices
Point(1) = {xl, yl, 0};
Point(2) = {xr, yl, 0};
Point(3) = {xr, yr, 0};
Point(4) = {xl, yr, 0};
For k In {0:nx:1}
    Point(k+5) = {xl + k, yc + eps*Cos(k*Pi), 0};
EndFor


// lines
Line(1) = {nx+5, 2};
Line(2) = {2, 1};
Line(3) = {1, 5};
Line(4) = {nx+5, 3};
Line(5) = {3, 4};
Line(6) = {4, 5};
For k In {0:nx-1:1}
    Line(k+7) = {k+5,k+5+1};
EndFor

// wires
beam1[] = {1,2,3}; // make wire from lines; u can turn around line with - (-tag)
For k In {0:nx-1:1}
    beam1[] += {k+7};
EndFor
beam2[] = {4,5,6}; // make wire from lines; u can turn around line with - (-tag)
For k In {0:nx-1:1}
    beam2[] += {k+7};
EndFor

Curve Loop(1) = beam1[];
Curve Loop(2) = beam2[];

// surfaces
Plane Surface(1) = {1}; // make surfaces from wires
Plane Surface(2) = {2};

// boundaries
//Physical Curve("Patch1") = {3,5}; // make 1D patch from lines
//Physical Curve("Patch2") = {2};
//Physical Curve("Patch3") = {1,4};
//Physical Curve("Patch4") = {5};
//Physical Surface("Volume") = {1,2}; // make 2D patch from surfaces

// mesh characteristics
Transfinite Curve {2} = nx+1 Using Progression 1;
Transfinite Curve {5} = nx+1 Using Progression 1;
Transfinite Curve {3,6} = 0.5*ny+1 Using Progression 1;
Transfinite Curve {1,4} = 0.5*ny+1 Using Progression 1;

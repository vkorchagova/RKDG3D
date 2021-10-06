// Gmsh project created on Tue Jan 29 12:30:02 2019

d = 5.35e-3;
lenX = 40.0;
lenY = 13.0;

nperd = 20;
nodes1 = lenX*nperd*0.5;
nodes2 = lenY*nperd*0.5;

// Points
Point(1) = {0.0, 0.5*d, 0.0, 1.0};
Point(2) = { -d, 0.5*d, 0.0, 1.0};
Point(3) = { -d, (lenY + 0.5)*d, 0.0, 1.0};
Point(4) = {0.0, (lenY + 0.5)*d, 0.0, 1.0};
Point(5) = {lenX*d, (lenY + 0.5)*d, 0.0, 1.0};
Point(6) = {lenX*d, 0.5*d, 0.0, 1.0};
Point(7) = {0.0, -0.5*d, 0.0, 1.0};
Point(8) = { -d, -0.5*d, 0.0, 1.0};
Point(9) = { -d, -(lenY + 0.5)*d, 0.0, 1.0};
Point(10) = {0.0, -(lenY + 0.5)*d, 0.0, 1.0};
Point(11) = {lenX*d, -(lenY + 0.5)*d, 0.0, 1.0};
Point(12) = {lenX*d, -0.5*d, 0.0, 1.0};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {4, 5};
Line(6) = {5, 6};
Line(7) = {6, 1};

Line(8) = {7, 8};
Line(9) = {8, 9};
Line(10) = {9, 10};
Line(11) = {10, 7};
Line(12) = {10, 11};
Line(13) = {11, 12};
Line(14) = {12, 7};

Line(15) = {1, 7};
Line(16) = {12, 6};

// Wires
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, -4};
Curve Loop(3) = {8, 9, 10, 11};
Curve Loop(4) = {12, 13, 14, -11};
Curve Loop(5) = {15, -14, 16, 7};

// Surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};
Plane Surface(4) = {4};
Plane Surface(5) = {5};

// Set boundaries
Physical Curve("inlet") = {15};
Physical Curve("outlet") = {6, 16, 13, 3, 5, 10, 12};
Physical Curve("wall") = {1, 2, 8, 9};
Physical Surface("Volume") = {1, 2, 3, 4, 5};

// Numbers of segments
Transfinite Curve {5, -7, -14, 12} = nodes1+1 Using Progression 1.002; // here is number of nodes!!!!!!
Transfinite Curve {3, -1, -8, 10} = nperd+1 Using Progression 1;
Transfinite Curve {15, -16} = nperd+1 Using Progression 1;
Transfinite Curve {2, -4, -6} = nodes2+1 Using Progression 1.005;
Transfinite Curve {9, -11, -13} = nodes2+1 Using Progression 1.005;

// Which surfaces are meshed
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};
Transfinite Surface {5};

// Make your mesh quadrangle
Recombine Surface {1, 2, 3, 4, 5};

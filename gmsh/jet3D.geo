// Gmsh project created for Jet 3D

d = 5.35e-3;

nperd = 2;
nodes1 = 10*nperd;
nodes2 = 13.0*nperd*0.5;

// Points
Point(1) = {0.0, 0.5*d, 0.0, 1.0};
Point(2) = {0.0, 13.5*d, 0.0, 1.0};
Point(3) = {20*d, 13.5*d, 0.0, 1.0};
Point(4) = {20*d, 0.5*d, 0.0, 1.0};
Point(5) = {0.0, -0.5*d, 0.0, 1.0};
Point(6) = {0.0, -13.5*d, 0.0, 1.0};
Point(7) = {20*d, -13.5*d, 0.0, 1.0};
Point(8) = {20*d, -0.5*d, 0.0, 1.0};

// Lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {4, 8};

// Wires
Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {5, 6, 7, 8};
Curve Loop(3) = {-9, -4, 10 ,8};

// Surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Set boundaries
Physical Curve("inlet") = {9};
Physical Curve("outlet") = {3, 10, 7, 2, 6};
Physical Curve("wall") = {1, 5};
Physical Surface("Volume") = {1, 2, 3};

// Numbers of segments
Transfinite Curve {2, 4, 8, 6} = nodes1 Using Progression 1; // here is number of nodes!!!!!!
Transfinite Curve {9, 10} = nperd Using Progression 1;
Transfinite Curve {1, -3} = nodes2 Using Progression 1.01;
Transfinite Curve {-5, 7} = nodes2 Using Progression 0.99;

// Which surfaces are meshed
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};

// Make your mesh quadrangle
Recombine Surface {1, 2, 3};

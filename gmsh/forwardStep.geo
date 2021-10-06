// Gmsh project created on Tue Jan 29 12:30:02 2019

cellsStep = 1;

nodes1 = cellsStep + 1;
nodes2 = cellsStep*4 + 1;
nodes3 = cellsStep*3 + 1;
nodes4 = cellsStep*12 + 1;

// Points
Point(1) = {0, 0, 0, 1.0};
Point(2) = {-3, 0, 0, 1.0};
Point(3) = {-3, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Point(5) = {-2.4, 0, 0, 1.0};
Point(6) = {-2.4, 0.2, 0, 1.0};
Point(7) = {0, 0.2, 0, 1.0};
Point(8) = {-2.4, 1, 0, 1.0};
Point(9) = {-3, 0.2, 0, 1.0};

// Lines
Line(1) = {2, 5};
Line(2) = {5, 6};
Line(3) = {6, 7};
Line(4) = {7, 4};
Line(5) = {4, 8};
Line(6) = {8, 3};
Line(7) = {3, 9};
Line(8) = {9, 2};
Line(9) = {9, 6};
Line(10) = {6, 8};

// Wires
Curve Loop(1) = {7, 9, 10, 6};
Curve Loop(2) = {8, 1, 2, -9};
Curve Loop(3) = {10, -5, -4, -3};

// Surfaces
Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// Set boundaries
Physical Curve("inlet") = {7, 8};
Physical Curve("outlet") = {4};
Physical Curve("obstacle") = {2, 3};
Physical Curve("wall") = {1, 6, 5};
Physical Surface("Volume") = {1, 2, 3};

// Numbers of segments
Transfinite Curve {8, 2} = nodes1 Using Progression 1; // here is number of nodes!!!!!!
Transfinite Curve {7, 10, 4} = nodes2 Using Progression 1;
Transfinite Curve {6, 9, 1} = nodes3 Using Progression 1;
Transfinite Curve {5, 3} = nodes4 Using Progression 1;

// Which surfaces are meshed
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};

// Make your mesh quadrangle
Recombine Surface {1, 2, 3};

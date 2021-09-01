// Gmsh project for bump test
// 2D bump 


x1 = 0.0;
x2 = 1.0;
x3 = 2.0;
x4 = 3.0;

y1 = 0.0;
y2 = 1.0;

xc = 1.5;
yc = -1.2;//-3.105;//-1.2;

nx1 = 20;
nx2 = 20;
nx3 = 20;
ny1 = 20;

// vertices
Point(1) = {x1, y1, 0};
Point(2) = {x2, y1, 0};
Point(3) = {x3, y1, 0};
Point(4) = {x4, y1, 0};
Point(5) = {x1, y2, 0};
Point(6) = {x2, y2, 0};
Point(7) = {x3, y2, 0};
Point(8) = {x4, y2, 0};
Point(9) = {xc, yc, 0};

// lines
Line(1) = {1,2};
Circle(2) = {2,9,3};
Line(3) = {3,4};
Line(4) = {5,6};
Line(5) = {6,7};
Line(6) = {7,8};
Line(7) = {1,5};
Line(8) = {2,6};
Line(9) = {3,7};
Line(10) = {4,8};

// wires
Curve Loop(1) = {1,8,-4,-7}; // make wire from lines; u can turn around line with - (-tag)
Curve Loop(2) = {2,9,-5,-8};
Curve Loop(3) = {3,10,-6,-9};


// surfaces

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3};

// boundaries
Physical Curve("inlet") = {7}; // make 1D patch from lines
Physical Curve("outlet") = {10}; // make 1D patch from lines
Physical Curve("bottom") = {1,2,3}; // make 1D patch from lines
Physical Curve("top") = {4,5,6}; // make 1D patch from lines

Physical Surface("Volume") = {1,2,3}; // make 2D patch from surfaces

// mesh characteristics
Transfinite Curve {1,4} = nx1+1 Using Progression 1;
Transfinite Curve {2,5} = nx2+1 Using Progression 1;
Transfinite Curve {3,6} = nx3+1 Using Progression 1;
Transfinite Curve {7,8,9,10} = ny1+1 Using Progression 1.01;
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};

// renumeration
Recombine Surface {1,2,3};


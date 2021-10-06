// Gmsh project for Sod quasi 1D
// Cutted plane for Double Mach reflection

// Scheme
// p5----------p4
//  |          |
//  |          p3
// p1--p2

x1 = 0; 	y1 = 0;
x2 = 0.2;   	y2 = 0;
x3 = 3.2;	y3 = 1.732051;
x4 = 3.2;	y4 = 2.2;
x5 = 0; 	y5 = 2.2;

cellStep = 1;

n1 = cellStep + 1;
n2 = 2*cellStep + 1;
n3 = 10*cellStep + 1;
n4 = 80*cellStep + 1;


// vertices
Point(1) = {x1, y1, 0};
Point(2) = {x2, y2, 0};
Point(3) = {x3, y3, 0};
Point(4) = {x4, y4, 0};
Point(5) = {x5, y5, 0};

// lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 1};

// wires
Curve Loop(1) = {5, 1, 2, 3, 4}; // make wire from lines; u can turn around line with - (-tag)

// surfaces
Plane Surface(1) = {1}; // make surfaces from wires


// boundaries
Physical Curve("inlet") = {5}; // make 1D patch from lines
Physical Curve("outlet") = {3};
Physical Curve("top") = {4};
Physical Curve("bottom") = {1,2};


Physical Surface("Volume") = {1}; // make 2D patch from surfaces

// mesh characteristics
Transfinite Curve {1} = n1 Using Progression 1;
Transfinite Curve {2} = n4 Using Progression 1;
Transfinite Curve {3} = n2 Using Progression 1;
Transfinite Curve {4} = n3 Using Progression 1;
Transfinite Curve {5} = n3 Using Progression 1;

//Transfinite Curve {1, 3} = nx+1 Using Progression 1;
//Transfinite Surface {1} = {1, 2, 3, 4};

// renumeration
//Recombine Surface {1};


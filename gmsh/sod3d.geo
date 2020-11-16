// Gmsh project for Sod 3D
// Hexahedron

xl = -1.0;
xr =  1.0;
yl = -1.0;
yr =  1.0;
zl = -1.0;
zr =  1.0;

// vertices
Point(1) = { xl, yl, zl, 1.0};
Point(2) = { xr, yl, zl, 1.0};
Point(3) = { xr, yr, zl, 1.0};
Point(4) = { xl, yr, zl, 1.0};
Point(5) = { xl, yl, zr, 1.0};
Point(6) = { xr, yl, zr, 1.0};
Point(7) = { xr, yr, zr, 1.0};
Point(8) = { xl, yr, zr, 1.0};

// lines
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {4, 3};
Line(4) = {1, 4};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {8, 7};
Line(8) = {5, 8};
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};

// cube sides
Line Loop(13) = {1, 2, -3, -4};
Plane Surface(14) = {13};
Line Loop(15) = {5, 6, -7, -8};
Plane Surface(16) = {15};
Line Loop(17) = {1, 10, -5, -9};
Plane Surface(18) = {17};
Line Loop(19) = {2, 11, -6, -10};
Plane Surface(20) = {19};
Line Loop(21) = {-3, 12, 7, -11};
Plane Surface(22) = {21};
Line Loop(23) = {-4, 9, 8, -12};
Plane Surface(24) = {23};

Surface Loop(25) = {18, 20, 16, 22, 14, 24};

// cube volume
Volume(26) = {25};


// boundaries
Physical Surface("top") = {16};
Physical Surface("bottom") = {14};
Physical Surface("sides") = {18, 20, 22, 24};
Physical Volume("internal") = {26};


// mesh characteristics
Transfinite Line "*" = 101 Using Progression 1.0;
Transfinite Surface "*";
Recombine Surface "*";
Transfinite Volume "*";
Recombine Volume "*";



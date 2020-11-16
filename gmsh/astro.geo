// Gmsh project for astro test
// 2D thorus 

// Scheme
//          p7
//     *    |    *
//   *      p3     *
//  *    *      *   *
// p8---p4  p1  p2---p6
//  *    *      *   *
//   *      p5     *
//     *    |    *
//          p9

xc = 0.0;
yc = 0.0;

r1 = 1.1;
r2 = 3.1;

nr = 100;
nphi = 50;

// vertices
Point(1) = {xc, yc, 0};
Point(2) = {r1, 0, 0};
Point(3) = {0, r1, 0};
Point(4) = {-r1, 0, 0};
Point(5) = {0, -r1, 0};
Point(6) = {r2, 0, 0};
Point(7) = {0, r2, 0};
Point(8) = {-r2, 0, 0};
Point(9) = { 0, -r2, 0};

// lines
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};
Circle(5) = {6,1,7};
Circle(6) = {7,1,8};
Circle(7) = {8,1,9};
Circle(8) = {9,1,6};
Line(9) = {2,6};
Line(10) = {3,7};
Line(11) = {4,8};
Line(12) = {5,9};

// wires
Curve Loop(1) = {-1,9,5,-10}; // make wire from lines; u can turn around line with - (-tag)
Curve Loop(2) = {-2,10,6,-11};
Curve Loop(3) = {-3,11,7,-12};
Curve Loop(4) = {-4,12,8,-9};

// surfaces

Plane Surface(1) = {1};
Plane Surface(2) = {2};
Plane Surface(3) = {3}; // make surfaces from wires
Plane Surface(4) = {4};

// boundaries
Physical Curve("Patch1") = {1,2,3,4}; // make 1D patch from lines
Physical Curve("Patch2") = {5,6,7,8};
Physical Surface("Volume") = {1, 2, 3, 4};

// mesh characteristics
Transfinite Curve {9,10,11,12} = nr+1 Using Progression 1;
Transfinite Curve {6,2} = nphi+1 Using Progression 1;
Transfinite Curve {7,3} = nphi+1 Using Progression 1;
Transfinite Curve {8,4} = nphi+1 Using Progression 1;
Transfinite Curve {5,1} = nphi+1 Using Progression 1;
Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};

// renumeration
Recombine Surface {1,2,3,4};


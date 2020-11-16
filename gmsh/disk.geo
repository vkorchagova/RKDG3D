// Gmsh project for astro test
// 2D disk 

// Scheme
//         p3
//     *        *
//   *            *
//  *              *
// p4      p1      p2
//  *              *
//   *            *
//     *        *
//         p5

xc = 0.0;
yc = 0.0;

r1 = 2.0;

nr = 100;
nphi = 50;

// vertices
Point(1) = {xc, yc, 0};
Point(2) = {r1, 0, 0};
Point(3) = {0, r1, 0};
Point(4) = {-r1, 0, 0};
Point(5) = {0, -r1, 0};

// lines
Circle(1) = {2,1,3};
Circle(2) = {3,1,4};
Circle(3) = {4,1,5};
Circle(4) = {5,1,2};

// wires
Curve Loop(1) = {1,2,3,4}; // make wire from lines; u can turn around line with - (-tag)


// surfaces

Plane Surface(1) = {1};

// boundaries
Physical Curve("Patch1") = {1}; // make 1D patch from lines
Physical Curve("Patch2") = {2}; // make 1D patch from lines
Physical Curve("Patch3") = {3}; // make 1D patch from lines
Physical Curve("Patch4") = {4}; // make 1D patch from lines
Physical Surface("Volume") = {1}; // make 2D patch from surfaces

// mesh characteristics
Transfinite Curve {1,2,3,4} = nr+1 Using Progression 1;
Transfinite Surface {1};

// renumeration
//Recombine Surface {1};


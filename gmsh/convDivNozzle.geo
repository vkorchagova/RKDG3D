// Gmsh project for conv-div noxxle


xl = -1.0;
xr =  1.0;
xc =  0.0;
yIn  = 0.075;
yOut = 0.225;

nx = 100;
ny = 1;

// vertices
Point(1) = {xl, -yOut, 0};
Point(2) = {xl,  yOut, 0};
Point(3) = {xc, -yIn, 0};
Point(4) = {xc,  yIn, 0};
Point(5) = {xr, -yOut, 0};
Point(6) = {xr,  yOut, 0};


// lines
Line(1) = {1, 2};
Line(2) = {2, 4};
Line(3) = {4, 6};
Line(4) = {6, 5};
Line(5) = {5, 3};
Line(6) = {3, 1};
Line(7) = {4, 3};

// wires
Curve Loop(1) = {1, 2, 7, 6}; // make wire from lines; u can turn around line with - (-tag)
Curve Loop(2) = {3, 4, 5, -7};

// surfaces
Plane Surface(1) = {1}; // make surfaces from wires
Plane Surface(2) = {2};

// boundaries
Physical Curve("inlet",2) = {1}; // make 1D patch from lines
Physical Curve("outlet",3) = {4};
Physical Curve("walls",4) = {2, 3, 5, 6};
Physical Surface("Volume") = {1, 2}; // make 2D patch from surfaces

// mesh characteristics
Transfinite Curve {1, 7, 4} = ny+1 Using Progression 1;
Transfinite Curve {2, 6} = 0.5*nx+1 Using Progression 1;
Transfinite Curve {3, 5} = 0.5*nx+1 Using Progression 1;
Transfinite Surface {1, 2};
//Transfinite Surface {2};


// renumeration
Recombine Surface {1, 2};


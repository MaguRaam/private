// Trapezoid Points 

Point(1) = {0.0, 0.0, 0.0};
Point(2) = {2.5, 1.44337567297, 0.0};
Point(3) = {2.5, 2.0, 0.0};
Point(4) = {0.0, 2-1.44337567297, 0.0};

// Lines Rectangle 

Line(5) = {1,2};
Line(6) = {2,3};
Line(7) = {3,4};
Line(8) = {4,1};



// Line Loops 

Line Loop(9) = {5,6,7,8};
Plane Surface(10)  = {9};

Ny = 30; 
Nx = 5*Ny;


Transfinite Line{5,7} = Nx;  // x-direction // (5/4)*Nx
Transfinite Line{6,8} = Ny; // y-direction 
Transfinite Surface{10};
Recombine Surface{10};

 


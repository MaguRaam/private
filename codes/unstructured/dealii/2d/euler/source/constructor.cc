#include "../include/Euler.h"


// Construtor for the Weno3_2D class 

Euler_2D::Euler_2D (double ft, double cfl_no)
  :
	dt(0.0), 
    finalTime(ft),
	cfl (cfl_no),
    fv (0),
    dof_handler (triangulation)
{}

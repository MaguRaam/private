#include "../include/Weno32.h"


// Construtor for the WENO4 class 

Weno3_2D::Weno3_2D (double ft, double cfl_no)
  :
	dt(0.0), 
    finalTime(ft),
	cfl (cfl_no),
    fv (0),
    dof_handler (triangulation)
{}

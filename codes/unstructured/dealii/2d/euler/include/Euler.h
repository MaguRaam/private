#ifndef WENO32_H_
#define WENO32_H_

#include "Headers.h"
#include "Riemann_Solvers.h"
#include "Exceptions.h"

// Complementary Function Definitions   

Vector<double> initial_condition(Point <2>);

const double M = 10.0; 

// Main Class Declaration 

class Euler_2D {
	
	void make_grid(); 
    void setup_system ();
	void initialize();
	void compute_rhs();
    void compute_time_step_based_on_cfl_number(double);
    void solve(); 
	void output_results (unsigned int);
    void restart(double, unsigned int); 

	Triangulation<2> triangulation;

	Vector<double> RHO;
	Vector<double> RHO_U;
	Vector<double> RHO_V;
	Vector<double> E;

	Vector<double> rhs1;
	Vector<double> rhs2;
	Vector<double> rhs3;
	Vector<double> rhs4;
 
    double dt;
	double finalTime;
	double cfl;
    double h_min; 

    FE_DGQ<2> fv;
	DoFHandler<2> dof_handler;

public:
 	Euler_2D(double, double);
	void run ();
};



#endif /* WENO32_H_ */

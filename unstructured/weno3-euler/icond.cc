/*
 * icond.cc
 *      Author: sunder
 */

#include "hype.hh"

// Some common test problems

void shu_osher_shock_tube(double, double, Array1D<double>&);
void double_mach_reflection(double, double, Array1D<double>&);


//-----------------------------------------------------------------------
// Initial condition function
//-----------------------------------------------------------------------

void icond(double x, double y, Array1D<double>& q0) {

	double_mach_reflection(x,y,q0);

	/*
	const double h = 0.00196969;
	const double rhoL = 1.0;
	const double rhoR = 0.125;
	const double pL = 1.0;
	const double pR = 0.1;
	const double x0 = 0.0;

	v0[0] = tanh_profile(rhoL, rhoR, x0, 3.0*h, x);
	v0[1] = 0.0;
	v0[2] = 0.0;
	v0[3] = tanh_profile(pL, pR, x0, 3.0*h, x);
	*/

	//double h = 0.0177913;
	//double spread = 5.0;
	/*
	v0[0] = tanh_profile(rhoL, rhoR, x0, spread*h, x);
	v0[1] = tanh_profile(uL, uR, x0, spread*h, x);
	v0[2] = 0.0;
	v0[3] = tanh_profile(pL, pR, x0, spread*h, x);
	*/



	/*
	double rho0 = 1.0;
	double vlx0 = 1.0;
	double vly0 = 1.0;
	double prs0 = 1.0;

	double kappa = 5.0; // Strength of the vortex

	double r2 = x*x + y*y;
	double exp_r2 = exp(0.5*(1.0 - r2));

	double tempaa = -exp_r2*exp_r2*kappa*kappa*(GAMMA - 1.0)/(8.0*GAMMA*M_PI*M_PI);
	double tempab = tempaa + prs0/rho0;

	tempab = tempab*std::pow(rho0,GAMMA)/prs0;

	v0[0] = std::pow(tempab, ( 1.0 / (GAMMA - 1.0)));
	v0[1] = vlx0 - y*kappa *exp_r2/(2.0*M_PI);
	v0[2] = vly0 + x*kappa*exp_r2/(2.0*M_PI);
	v0[3] = v0[0]*(tempaa + prs0/rho0);
	*/
	/*
	double M = 3.0;

	v0[0] = 1.4;
	v0[1] = M;
	v0[2] = 0.0;
	v0[3] = 1.0;
	*/

}

//-----------------------------------------------------------------------
// Shu-Osher Oscillatory Shock Tube
//-----------------------------------------------------------------------

void shu_osher_shock_tube(double x, double y, Array1D<double>& q0) {

	Array1D<double> v0(nVar);

	if (x < -4.0) {
		v0[0] = 3.857143;
		v0[1] = 2.629369;
		v0[2] = 0.0;
		v0[3] = 10.33333;
	}

	else {
		v0[0] = 1.0 + 0.2*sin(5.0*x);
		v0[1] = 0.0;
		v0[2] = 0.0;
		v0[3] = 1.0;
	}

	PDEPrim2Cons(v0, q0);
}

//-----------------------------------------------------------------------
// Double Mach Reflection Problem
//-----------------------------------------------------------------------

void double_mach_reflection(double x, double y, Array1D<double>& q0) {

	Array1D<double> v0(nVar);

	if (x < 0.25) {
		v0[0] = 8.0;
		v0[1] = 8.25;
		v0[2] = 0.0;
		v0[3] = 116.5;
	}

	else {
		v0[0] = 1.4;
		v0[1] = 0.0;
		v0[2] = 0.0;
		v0[3] = 1.0;
	}

	PDEPrim2Cons(v0, q0);
}


//-----------------------------------------------------------------------
// Exact solution (for checking convergence)
//-----------------------------------------------------------------------

void exact_solution(double x, double y, double t, double* q) {


	double rho0 = 1.0;
	double vlx0 = 1.0;
	double vly0 = 1.0;
	double prs0 = 1.0;

	double kappa = 5.0; // Strength of the vortex

	double r2 = (x-t)*(x-t) + (y-t)*(y-t);
	double exp_r2 = exp(0.5*(1.0 - r2));

	double tempaa = -exp_r2*exp_r2*kappa*kappa*(GAMMA - 1.0)/(8.0*GAMMA*M_PI*M_PI);
	double tempab = tempaa + prs0/rho0;

	tempab = tempab*std::pow(rho0,GAMMA)/prs0;

	q[0] = std::pow(tempab, ( 1.0 / (GAMMA - 1.0)));
	q[1] = rho0*vlx0;
	q[2] = rho0*vly0;
	q[3] = rho0*vlx0;

}

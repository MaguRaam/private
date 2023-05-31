/*
 * path_conservative_fv.hh
 *      Author: sunder
 */ 

#ifndef PATH_CONSERVATIVE_FV_HH_
#define PATH_CONSERVATIVE_FV_HH_

#include "riemann_solver.hh"
#include "weno.hh"
#include "quadrature.hh"
#include "polynomials.hh"

enum bc_type {
	inflow,
	transmissive,
	reflective,
	periodic
};

// Enumerator for output format

enum data_format {
	vtk,
	tecplot
};

struct AppCtx {
	int N_x;
	int N_y;
	int order;
	double CFL;
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	double initial_time = 0.0;
	double final_time;
	int write_interval = 5;
	double p1 = 0.0;
	double p2 = 0.0;
	double g1 = 1.4;
	double g2 = 1.4; 
	bc_type left_boundary   = periodic;
	bc_type right_boundary  = periodic;
	bc_type top_boundary    = periodic;
	bc_type bottom_boundary = periodic;
	bool reconstruct_primitive_variables = false; 
	data_format output_format = vtk;
};

class Path_Conservative_FV {
	
	typedef boost::multi_array<double, 2> array_type;
	typedef array_type::index index;
	
	static const int dim = 2;
	static const int N_ph = 3; // No. of ghost points in each direction
	
	std::vector<double> (*init_func)(double, double);
	AppCtx Params;
	Legendre_m_2D FV;
	Riemann_Solver riemann;
	
	boost::multi_array<std::vector<double>, 3> U;
	boost::multi_array<std::vector<double>, 3> W;
	boost::multi_array<double, 3> F;
	boost::multi_array<double, 3> G;
	boost::multi_array<double, 3> D;
	boost::multi_array<double, 3> E;
	boost::multi_array<double, 3> RHS;
	boost::multi_array<double, 4> U_L;
	boost::multi_array<double, 4> U_R;
	boost::multi_array<double, 4> U_B;
	boost::multi_array<double, 4> U_T;
	
	boost::multi_array<double, 1> x;
	boost::multi_array<double, 1> y;

	double dx;
	double dy;
	double dt;
	double time;
	int time_step;
	int rk_stage; 
	
	void initialize();
	void compute_primitive_variables(); 
	void apply_boundary_conditions(); 
	void compute_rhs(double);
	void solve_ssprk22();
	void solve_ssprk33();
	void compute_errors(double&, double&);
	void plot_tecplot(int = 0, const int = 4); 
	

	 
	
public:
	Path_Conservative_FV(std::vector<double>(*)(double, double), AppCtx);
	void run();
	
};






#endif /* PATH_CONSERVATIVE_FV_HH_ */

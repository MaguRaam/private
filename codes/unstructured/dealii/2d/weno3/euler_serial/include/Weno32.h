#ifndef WENO43_H_
#define WENO43_H_

#include "Headers.h"
#include "CLS.h"
#include "LU.h"
#include "Riemann_Solvers.h"
#include "Exceptions.h"
#include "Slope_Limiters.h"
#include "stencil_selector.h"
#include "vertex_to_cell.h"
#include "cell_properties.h"

// Complementary Function Definitions   

Vector<double> initial_condition(Point<2>);

double evaluate_weno_polynomial(Vector<double>, Vector<double>,  Point<2>); 

double compute_second_order_smoothness_indicator(Vector<double>, Vector<double>, double);
double compute_third_order_smoothness_indicator(Vector<double>, Vector<double>, double);

DoFHandler<2>::active_cell_iterator return_cell_pointer(const DoFHandler<2>&, unsigned int);

// Main Class Declaration 

class Weno3_2D {
    
    void make_grid();  
    void setup_system ();
    void precompute_matrices();
    void precompute_matrices_veclocity();
    void compute_IS_constants(); 
    void compute_weno_polynomial_constants();
    void reconstruct();
	void initialize();
	void compute_rhs();
    void compute_time_step_based_on_cfl_number(double);
    void solve_ssprk33();
	void compute_cell_properties(); 
	void select_stencils();
	void post_process(unsigned int); 
    void solve_ader(); 
	void output_results (unsigned int);
    void restart(double, unsigned int); 

	Triangulation<2> triangulation;

	Vector<double> RHO;
	Vector<double> RHO_U;
	Vector<double> RHO_V;
	Vector<double> E;

    // Coefficients for WENO polynomials
    std::vector< Vector<double> >  coeffs_RHO;
    std::vector< Vector<double> >  coeffs_RHO_U;
    std::vector< Vector<double> >  coeffs_RHO_V;
    std::vector< Vector<double> >  coeffs_E;
    
    
    // WENO polynomial constants (only depend on mesh)
    std::vector< Vector<double> >  WENO_poly_consts;
    std::vector< Vector<double> >  IS_constants;
    std::vector<bool> is_corner_cell;
	std::vector<cell_properties> Cell; 
	std::vector<stencil_selector> Stencil; 
    

	Vector<double> rhs1;
	Vector<double> rhs2;
	Vector<double> rhs3;
	Vector<double> rhs4;
    
    // Third Order Stencils 
    std::vector< bool > is_admissible_R3;  // Centered third order stencil stencil
    std::vector< Constrained_LS > CLS_R3; 
    
    std::vector< Constrained_LS > CLS_R3_slip;
    std::vector< Constrained_LS > CLS_R2_slip;
    
    std::map <unsigned int, unsigned int> wall_boundary_global_index_map;
    
    std::vector< LUdcmp > LU_R21;
    std::vector< LUdcmp > LU_R22;
    std::vector< LUdcmp > LU_R23;
    std::vector< LUdcmp > LU_R24;
    
 
    double dt;
	double finalTime;
	double cfl;
    double h_min; 

    FE_DGQ<2> fv;
	DoFHandler<2> dof_handler;

public:
 	Weno3_2D(double, double);
	void run ();
};



#endif /* WENO43_H_ */

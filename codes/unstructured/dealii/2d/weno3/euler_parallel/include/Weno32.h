#ifndef WENO43_H_
#define WENO43_H_

#include "Headers.h"
#include "CLS.h"
#include "LU.h"
#include "Riemann_Solvers.h"
#include "Exceptions.h"
#include "Slope_Limiters.h"
#include "cell_properties.h"

const double M = 0.18; 

// Complementary Function Definitions   

Vector<double> initial_condition(Point<2>);

Vector<double> solve_system(FullMatrix<double>, Vector<double>);

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
    void compute_time_step_based_on_cfl_number();
    void solve_ssprk33();
	void compute_cell_properties(); 
    void solve_ader(); 
	void output_results ();
    void restart(); 
    void restart_r();
	void output_grid();
	void plot_tecplot() const; 

    MPI_Comm                                  mpi_communicator;

    parallel::distributed::Triangulation<2> triangulation;

    IndexSet                                  locally_owned_dofs;
    IndexSet                                  locally_relevant_dofs;

    PETScWrappers::MPI::Vector                           RHO;
    PETScWrappers::MPI::Vector                           RHO_U;
    PETScWrappers::MPI::Vector                           RHO_V;
    PETScWrappers::MPI::Vector                           E;

    PETScWrappers::MPI::Vector                           local_RHO;
    PETScWrappers::MPI::Vector                           local_RHO_U;
    PETScWrappers::MPI::Vector                           local_RHO_V;
    PETScWrappers::MPI::Vector                           local_E;

    // Coefficients for WENO polynomials
    PETScWrappers::MPI::Vector	coeffs_x_RHO;
    PETScWrappers::MPI::Vector	coeffs_y_RHO;
    PETScWrappers::MPI::Vector	coeffs_xx_RHO;
    PETScWrappers::MPI::Vector	coeffs_yy_RHO;
    PETScWrappers::MPI::Vector	coeffs_xy_RHO;

    PETScWrappers::MPI::Vector	coeffs_x_RHO_U;
    PETScWrappers::MPI::Vector	coeffs_y_RHO_U;
    PETScWrappers::MPI::Vector	coeffs_xx_RHO_U;
    PETScWrappers::MPI::Vector	coeffs_yy_RHO_U;
    PETScWrappers::MPI::Vector	coeffs_xy_RHO_U;

    PETScWrappers::MPI::Vector	coeffs_x_RHO_V;
    PETScWrappers::MPI::Vector	coeffs_y_RHO_V;
    PETScWrappers::MPI::Vector	coeffs_xx_RHO_V;
    PETScWrappers::MPI::Vector	coeffs_yy_RHO_V;
    PETScWrappers::MPI::Vector	coeffs_xy_RHO_V;

    PETScWrappers::MPI::Vector	coeffs_x_E;
    PETScWrappers::MPI::Vector	coeffs_y_E;
    PETScWrappers::MPI::Vector	coeffs_xx_E;
    PETScWrappers::MPI::Vector	coeffs_yy_E;
    PETScWrappers::MPI::Vector	coeffs_xy_E;

    PETScWrappers::MPI::Vector	local_coeffs_x_RHO;
    PETScWrappers::MPI::Vector	local_coeffs_y_RHO;
    PETScWrappers::MPI::Vector	local_coeffs_xx_RHO;
    PETScWrappers::MPI::Vector	local_coeffs_yy_RHO;
    PETScWrappers::MPI::Vector	local_coeffs_xy_RHO;

    PETScWrappers::MPI::Vector	local_coeffs_x_RHO_U;
    PETScWrappers::MPI::Vector	local_coeffs_y_RHO_U;
    PETScWrappers::MPI::Vector	local_coeffs_xx_RHO_U;
    PETScWrappers::MPI::Vector	local_coeffs_yy_RHO_U;
    PETScWrappers::MPI::Vector	local_coeffs_xy_RHO_U;

    PETScWrappers::MPI::Vector	local_coeffs_x_RHO_V;
    PETScWrappers::MPI::Vector	local_coeffs_y_RHO_V;
    PETScWrappers::MPI::Vector	local_coeffs_xx_RHO_V;
    PETScWrappers::MPI::Vector	local_coeffs_yy_RHO_V;
    PETScWrappers::MPI::Vector	local_coeffs_xy_RHO_V;

    PETScWrappers::MPI::Vector	local_coeffs_x_E;
    PETScWrappers::MPI::Vector	local_coeffs_y_E;
    PETScWrappers::MPI::Vector	local_coeffs_xx_E;
    PETScWrappers::MPI::Vector	local_coeffs_yy_E;
    PETScWrappers::MPI::Vector	local_coeffs_xy_E; 
    
    // WENO polynomial constants (only depend on mesh)
	PETScWrappers::MPI::Vector WENO_poly_consts_x;
	PETScWrappers::MPI::Vector WENO_poly_consts_y;
	PETScWrappers::MPI::Vector WENO_poly_consts_xx;
	PETScWrappers::MPI::Vector WENO_poly_consts_yy;
	PETScWrappers::MPI::Vector WENO_poly_consts_xy;

    std::vector< Vector<double> >  IS_constants;
    std::vector<bool> is_corner_cell;
	std::vector<cell_properties> Cell;   

	Vector<double> rhs1;
	Vector<double> rhs2;
	Vector<double> rhs3;
	Vector<double> rhs4;
    
    // Third Order Stencils 
    std::vector< bool > is_admissible_R3;  // Centered third order stencil stencil
    std::vector< Constrained_LS > CLS_R3; 
    
    std::vector< Constrained_LS > CLS_R3_slip;
    std::vector< Constrained_LS > CLS_R2_slip;

    std::map <unsigned int, unsigned int> global_to_local_index_map;    
    std::map <unsigned int, unsigned int> wall_boundary_global_index_map;
	std::vector< std::set<DoFHandler<2>::active_cell_iterator> > cell_neighbor_iterator;
	std::vector< std::set<DoFHandler<2>::active_cell_iterator> > cell_all_neighbor_iterator;
    
    std::vector< LUdcmp > LU_R21;
    std::vector< LUdcmp > LU_R22;
    std::vector< LUdcmp > LU_R23;
    std::vector< LUdcmp > LU_R24;   
 
	unsigned int dofs_per_cell ;
	unsigned int n_locally_cells;
	unsigned int n_vertices;

	std::vector<types::global_dof_index> local_dof_indices;
	std::vector<types::global_dof_index> local_neighbor_dof_indices;

    ConditionalOStream                        pcout;
    TimerOutput                               computing_timer;

    double dt;
	double finalTime;
	double cfl;
    double h_min; 
	double time;

	bool RESTART;

    FE_DGQ<2> fv;
	DoFHandler<2> dof_handler;

public:
 	Weno3_2D(double, double, bool);
	void run ();
};



#endif /* WENO43_H_ */

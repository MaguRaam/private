#include "../include/Euler.h"


// Initialize the solution 

void Euler_2D::initialize() {

    std::cout << "Initializing the solution" << std::endl;

    unsigned int N_gp = 4;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    Point<2> q_point;

    double V0;

    Vector<double> U(4); Vector<double> W(4);
    
    h_min = std::sqrt(cell->measure()); 

    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
		V0 = cell->measure();
		fv_values.reinit(cell);
		
		if (std::sqrt(cell->measure()) < h_min) {
		
			h_min = std::sqrt(cell->measure()); 
		}

		RHO(c) = 0.0; RHO_U(c) = 0.0; RHO_V(c) = 0.0; E(c) = 0.0; 

		for (unsigned int i = 0; i < N_gp*N_gp; i++) {

			q_point = fv_values.quadrature_point(i);

			W = initial_condition(q_point); 
			U = primitive_to_conserved(W); 

			RHO(c)   += (1./V0)*fv_values.JxW (i)*U(0);
			RHO_U(c) += (1./V0)*fv_values.JxW (i)*U(1);
			RHO_V(c) += (1./V0)*fv_values.JxW (i)*U(2);
			E(c)     += (1./V0)*fv_values.JxW (i)*U(3);
		}
    }
    
    std::cout << "Done!" << std::endl;  
    std::cout << "===========================" << std::endl;
} 

#include "../include/Weno32.h"


//  Evaluate time step using the CFL condition 

void Weno3_2D::compute_time_step_based_on_cfl_number(double time) {

	Vector<double> u(dof_handler.n_dofs());
	Vector<double> v(dof_handler.n_dofs());
	Vector<double> a(dof_handler.n_dofs());

	// Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
    
    Vector<double> U(4); Vector<double> W(4); 

    
	for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        U(0) = RHO(c); U(1) = RHO_U(c); U(2) = RHO_V(c); U(3) = E(c); 
        
        W = conserved_to_primitive(U);  

        a(c) = std::sqrt(GAMMA*W(3)/W(0));
        u(c) = a(c) + std::abs(W(1));
        v(c) = a(c) + std::abs(W(2));
        
	}
	
	double u_max = u.linfty_norm();
	double v_max = v.linfty_norm();

	dt = (cfl*h_min)/(sqrt(2.0)*sqrt(u_max*u_max + v_max*v_max));

	if((time + dt)>finalTime) {
		dt = finalTime- time;
	}
} 

#include "../include/Weno32.h"


//  Evaluate time step using the CFL condition 

void Weno3_2D::compute_time_step_based_on_cfl_number() {

	Vector<double> u(n_locally_cells);
	Vector<double> v(n_locally_cells);
	Vector<double> a(n_locally_cells);

	// Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
    
    Vector<double> U(4); Vector<double> W(4); 
	unsigned int g_i;
    
	for (; cell != endc; ++cell) {

	if (cell->is_locally_owned()){

        cell->get_dof_indices(local_dof_indices);
		g_i = local_dof_indices[0];
		unsigned int c = global_to_local_index_map[local_dof_indices[0] ];
        
        U(0) = RHO(g_i); U(1) = RHO_U(g_i); U(2) = RHO_V(g_i); U(3) = E(g_i); 
        
        W = conserved_to_primitive(U);  

        a(c) = std::sqrt(GAMMA*W(3)/W(0));
        u(c) = a(c) + std::abs(W(1));
        v(c) = a(c) + std::abs(W(2));
	
	}        
	}
	
	double u_max_local = u.linfty_norm();
	double v_max_local = v.linfty_norm();

	double u_max = Utilities::MPI::max (u_max_local, MPI_COMM_WORLD);
	double v_max = Utilities::MPI::max (v_max_local, MPI_COMM_WORLD);

	dt = (cfl*h_min)/(sqrt(2.0)*sqrt(u_max*u_max + v_max*v_max));

	if((time + dt)>finalTime) {
		dt = finalTime- time;
	}
} 

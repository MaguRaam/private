#include "../include/Weno432.h"


//  Evaluate time step using the CFL condition 

void Weno4_2D::compute_time_step_based_on_cfl_number() {

	Vector<double> u(n_locally_cells);
	Vector<double> v(n_locally_cells);
	Vector<double> a(n_locally_cells);
    
    Vector<double> U(4); Vector<double> W(4); 

	unsigned int g_i;
    
	for (unsigned int c = 0; c < n_locally_cells; ++c) {

		g_i = local_to_global_index_map[c];
	        
        U(0) = RHO(g_i); U(1) = RHO_U(g_i); U(2) = RHO_V(g_i); U(3) = E(g_i); 
	        
        W = conserved_to_primitive(U);  

        a(c) = std::sqrt(GAMMA*W(3)/W(0));
        u(c) = a(c) + std::abs(W(1));
        v(c) = a(c) + std::abs(W(2));

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

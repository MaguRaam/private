#include "../include/Weno32.h"

// Solve using SSPRK(3,3)

void Weno3_2D::solve_ssprk33() {

    double time = 0.0;
    unsigned int count = 0;
    
    Vector<double> RHO_old(dof_handler.n_dofs());
    Vector<double> RHO_U_old(dof_handler.n_dofs());
	Vector<double> RHO_V_old(dof_handler.n_dofs());
	Vector<double> E_old(dof_handler.n_dofs()); 

	while (time < finalTime) {

		compute_time_step_based_on_cfl_number(time);
		time += dt;
        
        if (count%500== 0) {
			post_process(count);
            output_results(count);
            restart(time, count); 
        }
        
        
        std::cout << "time = " << time << ", Final time = " << finalTime<< std::endl;
        count++;

        RHO_old = RHO;
        RHO_U_old = RHO_U;
        RHO_V_old = RHO_V;
        E_old = E;

		// SSPRK Step 1

        compute_rhs();

		for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
            RHO(c) = RHO(c) + dt*rhs1(c);
            RHO_U(c) = RHO_U(c) + dt*rhs2(c);
            RHO_V(c) = RHO_V(c) + dt*rhs3(c);
            E(c) = E(c) + dt*rhs4(c);
        
		}

		// SSPRK Step 2
		
        
		compute_rhs();

        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
            RHO(c)   = (3./4.)*RHO_old(c)   + (1./4.)*RHO(c)   + (1./4.)*dt*rhs1(c);
            RHO_U(c) = (3./4.)*RHO_U_old(c) + (1./4.)*RHO_U(c) + (1./4.)*dt*rhs2(c);
            RHO_V(c) = (3./4.)*RHO_V_old(c) + (1./4.)*RHO_V(c) + (1./4.)*dt*rhs3(c);
            E(c)     = (3./4.)*E_old(c)     + (1./4.)*E(c)     + (1./4.)*dt*rhs4(c);
        
        }

		// SSPRK Step 3

		compute_rhs();

        for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
            RHO(c)   = (1./3.)*RHO_old(c)   + (2./3.)*RHO(c)   + (2./3.)*dt*rhs1(c);
            RHO_U(c) = (1./3.)*RHO_U_old(c) + (2./3.)*RHO_U(c) + (2./3.)*dt*rhs2(c);
            RHO_V(c) = (1./3.)*RHO_V_old(c) + (2./3.)*RHO_V(c) + (2./3.)*dt*rhs3(c);
            E(c)     = (1./3.)*E_old(c)     + (2./3.)*E(c)     + (2./3.)*dt*rhs4(c);
            
		}
		
    }
    
    post_process(count); 
    output_results(count);
    restart(time, count); 

}
 

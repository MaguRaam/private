#include "../include/Euler.h"


void Euler_2D::solve() {

    // Using SSPRK(2,2)

	double time = 0.0;
    unsigned int count = 0;

	while (time  < finalTime) {

		compute_time_step_based_on_cfl_number(time);
		time += dt;
        
        if (count%100 == 0) {
            output_results(count);
        }
        
        std::cout << "time = " << time << ", Final time = " << finalTime<< std::endl;
        count++;

        compute_rhs();

		for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
		
			RHO(c) = RHO(c) + dt*rhs1(c);
			RHO_U(c) = RHO_U(c) + dt*rhs2(c);
			RHO_V(c) = RHO_V(c) + dt*rhs3(c);
			E(c) = E(c) + dt*rhs4(c);
            
		}
    
    }
    
    output_results(count);
}
 

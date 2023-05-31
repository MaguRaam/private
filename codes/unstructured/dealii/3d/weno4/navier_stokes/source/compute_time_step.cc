#include "../include/Weno432.h"


//  Evaluate time step using the CFL condition 

void Weno4_3D::compute_time_step_based_on_cfl_number() {

	{
    
	    Vector<double> U(5); Vector<double> W(5);
		double lambda_v_max_local = 0.0;
		double lambda_c_max_local = 0.0;
		double lambda_v, lambda_c;
		unsigned int g_i;
		double a, mu, T;	// 10 is a random number
//    	pcout<<"dt start: \n";
		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];
		        
    	    U(0) = RHO(g_i); U(1) = RHO_U(g_i); U(2) = RHO_V(g_i); U(3) = RHO_W(g_i); U(4) = E(g_i); 
	
//			std::cout<<"U: "<<U<<std::endl;
		        
    	    claw.conserved_to_primitive(U,W);
	
    	    a = std::sqrt(claw.gamma()*W(4)/W(0));
			T = W(4)/W(0); 
			mu = claw.viscosity(T);
			lambda_c = std::sqrt(W[1]*W[1] + W[2]*W[2] + W[3]*W[3]) + a;
			lambda_v = std::max(1.33*mu, claw.gamma()*mu/claw.Pr())/U(0); // 10 is a random number
/*
			if(	U(0) < 1e-10)
				std::cout<<"g_i: "<<g_i<<"\trho: "<<U(0)<<std::endl;

			if(	lambda_v > 10)
				std::cout<<"rho: "<<U(0)<<"\t"<<claw.Pr()<<"\t"<<claw.gamma()*mu<<std::endl;
*/
			lambda_c_max_local = std::max(lambda_c, lambda_c_max_local);
			lambda_v_max_local = std::max(lambda_v, lambda_v_max_local);	
		}
	
		double lambda_c_max = Utilities::MPI::max (lambda_c_max_local, MPI_COMM_WORLD);
		double lambda_v_max = Utilities::MPI::max (lambda_v_max_local, MPI_COMM_WORLD);
	
		dt = (cfl*h_min)/(lambda_c_max + (lambda_v_max/h_min)*2.0);
//		pcout<<"lambda_c_max: "<<lambda_c_max<<"\tlambda_v_max: "<<lambda_v_max<<std::endl;
		if((time + dt)>finalTime) {
			dt = finalTime- time;
		}
	}
} 

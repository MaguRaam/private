#include "../include/Weno432.h"

// Solve using SSPRK(3,3)

void Weno4_2D::solve_ssprk33() {

	pcout<<"solve by ssprk33: "<<std::endl;

	Use_ader = false;

    auto start = std::chrono::system_clock::now();

    unsigned int count = 0;
    
    Vector<double> RHO_old(n_locally_cells);
    Vector<double> RHO_U_old(n_locally_cells);
	Vector<double> RHO_V_old(n_locally_cells);
	Vector<double> E_old(n_locally_cells); 

	compute_time_step_based_on_cfl_number();

	unsigned int output_count;
	output_count = 0.2 / dt ;

	unsigned int g_i ;

	while (time < finalTime) {

	    auto start_ssprk33 = std::chrono::system_clock::now();
        
		compute_time_step_based_on_cfl_number();

		if (count%output_count == 0 ) {
            output_results();
        }

		if (count%100 == 0 ) {
            restart();
        }

		if ((count-5)%100 == 0 ) {
            restart_r();
        }

		time += dt;
        
    	if ( (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) && ( ( count%50 == 0 ) || std::fabs(time - finalTime) < 1e-8) ){
		    std::ofstream fout_convergence ;
    	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
    	    fout_convergence.precision(7) ;
	
    	    const std::string filename = "log.dat";
		    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    		fout_convergence <<"time = " << time <<",\t dt: "<<dt<< ",\t Final time = " << finalTime<< std::endl;
    	    fout_convergence.close();
		}

		pcout<<"time = " << time <<",\t dt: "<<dt<< ",\t Final time = " << finalTime<< std::endl;
	
        count++;

		for (unsigned int c = 0; c < n_locally_cells; ++c) {
			g_i = local_to_global_index_map[c];
			RHO_old(c) = RHO(g_i);
			RHO_U_old(c) = RHO_U(g_i);
			RHO_V_old(c) = RHO_V(g_i);
			E_old(c) = E(g_i);	
		}

		// SSPRK Step 1

        compute_rhs();

//		pcout<<"1st step"<<std::endl;

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

            local_RHO(g_i) = local_RHO(g_i) + dt*rhs1(c);
   	        local_RHO_U(g_i) = local_RHO_U(g_i) + dt*rhs2(c);
   	        local_RHO_V(g_i) = local_RHO_V(g_i) + dt*rhs3(c);
   	        local_E(g_i) = local_E(g_i) + dt*rhs4(c);
   
		}
    
		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;


		// SSPRK Step 2
		
        
		compute_rhs();

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

            local_RHO(g_i)   = (3./4.)*RHO_old(c)   + (1./4.)*local_RHO(g_i)   + (1./4.)*dt*rhs1(c);
            local_RHO_U(g_i) = (3./4.)*RHO_U_old(c) + (1./4.)*local_RHO_U(g_i) + (1./4.)*dt*rhs2(c);
            local_RHO_V(g_i) = (3./4.)*RHO_V_old(c) + (1./4.)*local_RHO_V(g_i) + (1./4.)*dt*rhs3(c);
            local_E(g_i)     = (3./4.)*E_old(c)     + (1./4.)*local_E(g_i)     + (1./4.)*dt*rhs4(c);        

        }

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;
        
		// SSPRK Step 3

		compute_rhs();

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];
 
            local_RHO(g_i)   = (1./3.)*RHO_old(c)   + (2./3.)*local_RHO(g_i)   + (2./3.)*dt*rhs1(c);
   	        local_RHO_U(g_i) = (1./3.)*RHO_U_old(c) + (2./3.)*local_RHO_U(g_i) + (2./3.)*dt*rhs2(c);
   	        local_RHO_V(g_i) = (1./3.)*RHO_V_old(c) + (2./3.)*local_RHO_V(g_i) + (2./3.)*dt*rhs3(c);
   	        local_E(g_i)     = (1./3.)*E_old(c)     + (2./3.)*local_E(g_i)     + (2./3.)*dt*rhs4(c);        

		}

	    auto start_com = std::chrono::system_clock::now();

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;

	    auto end_com = std::chrono::system_clock::now();
	    std::chrono::duration<double> elapsed_seconds_com = end_com - start_com;
	    std::chrono::duration<double> elapsed_seconds_ssprk33 = end_com - start_ssprk33;

    	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0 && count == 2){
		    std::ofstream fout_convergence ;
    	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
    	    fout_convergence.precision(7) ;
	
    	    const std::string filename = "timer.dat";
		    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);
	
    		fout_convergence << "time taken by data transfer of ssprk33 = " << elapsed_seconds_com.count() << std::endl;
    		fout_convergence << "time taken by 1 step of ssprk33 = " << elapsed_seconds_ssprk33.count() << std::endl;
    	    fout_convergence.close();
		}

    }

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "timer.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "time taken by ssprk33 = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}
    
	output_results();
    restart(); 
}
 

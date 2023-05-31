#include "../include/Weno432.h"


// Update the solution using SSPRK(5,4) method

void Weno4_3D::solve_rk6() {

	Use_ader = false;

    // Using SSPRK(5,4)

    auto start = std::chrono::system_clock::now();

    unsigned int count = 0;

    Vector<double> RHO_old(n_locally_cells);
    Vector<double> RHO_U_old(n_locally_cells);
    Vector<double> RHO_V_old(n_locally_cells);
    Vector<double> RHO_W_old(n_locally_cells);
    Vector<double> E_old(n_locally_cells);

    std::vector<Vector<double> > K1(7); 
    std::vector<Vector<double> > K2(7);
    std::vector<Vector<double> > K3(7);
    std::vector<Vector<double> > K4(7);
    std::vector<Vector<double> > K5(7);

    for (unsigned int i = 0; i < 7; ++i) {
        K1[i].reinit(n_locally_cells);
        K2[i].reinit(n_locally_cells);
        K3[i].reinit(n_locally_cells);
        K4[i].reinit(n_locally_cells);
        K5[i].reinit(n_locally_cells);
    }

	Vector<double> B_RK(7), C_RK(7);
	FullMatrix<double> A_RK(7,7); 
	A_RK = 0.0;

	// Stage - 1 ;
	C_RK[0] = 0.0 ; B_RK[0] = 9.0/180.0 ;
	// Stage - 2 ;
	C_RK[1] = 1.0 ; B_RK[1] = 0.0 ; A_RK[1][0] = 1.0 ;
	// Stage - 3 ;
	C_RK[2] = 0.5 ; B_RK[2] = 64.0/180.0 ; A_RK[2][0] = 3.0/8.0 ; A_RK[2][1] = 1.0/8.0 ;
	// Stage - 4 ;
	C_RK[3] = 2.0/3.0 ; B_RK[3] = 0.0 ; A_RK[3][0] = 8.0/27.0 ; A_RK[3][1] = 2.0/27.0 ; A_RK[3][2] = 8.0/27.0 ;
	// Stage - 5 ;
	C_RK[4] = (7.0 - sqrt(21.0) )/14.0 ; B_RK[4] = 49.0/180.0 ; A_RK[4][0] = (9.0*sqrt(21.0) - 21.0)/392.0 ; A_RK[4][1] = -8.0*(7.0 - sqrt(21.0))/392.0 ; A_RK[4][2] = 48.0*(7.0 - sqrt(21.0))/392.0 ; A_RK[4][3] = -3.0*(21.0 - sqrt(21.0))/392.0 ;
	// Stage - 6 ;
	C_RK[5] = (7.0 + sqrt(21.0) )/14.0 ; B_RK[5] = 49.0/180.0 ; A_RK[5][0] = -5.0*(231.0 + 51.0*sqrt(21.0))/1960.0 ; A_RK[5][1] = -40.0*(7.0 + sqrt(21.0))/1960.0 ; A_RK[5][2] = -320.0*sqrt(21.0)/1960.0 ; A_RK[5][3] = 3.0*(21.0 + 121.0*sqrt(21.0))/1960.0 ; A_RK[5][4] = 392.0*(6.0 + sqrt(21.0))/1960.0 ;
	// Stage - 7 ;
	C_RK[6] = 1.0 ; B_RK[6] = 9.0/180.0 ; A_RK[6][0] = 15.0*(22.0 + 7.0*sqrt(21.0))/180.0 ; A_RK[6][1] = 120.0/180.0 ; A_RK[6][2] = 40.0*(7.0*sqrt(21.0) - 5.0)/180.0 ; A_RK[6][3] = -63.0*(3.0*sqrt(21.0) - 2.0)/180.0 ; A_RK[6][4] = -14.0*(49.0 + 9.0*sqrt(21.0))/180.0 ; A_RK[6][5] = 70.0*(7.0 - sqrt(21.0))/180.0 ;


	compute_time_step_based_on_cfl_number();
	unsigned int output_count;
	output_count = 0.05 / dt ;
	if(output_count < 10) output_count = 100;
	pcout<<"dt: "<<dt<<"\toutput_count: "<<output_count<<std::endl;
	reconstruct();
	unsigned int g_i ;
	
	while (time < finalTime) {

	    auto start_rk6 = std::chrono::system_clock::now();

		compute_time_step_based_on_cfl_number();

		if (count%output_count == 0 ) {
            output_results();
			post_process_force();
			reconstruct_quad();
			reconstruct_check();
        }



//		pcout<<count <<" output_results"<<std::endl;
		if (count%100 == 0 ) {
            restart();
			reconstruct_quad();
			reconstruct_check();
        }
//		pcout<<count <<" reconstruct_quad"<<std::endl;
		if ((count - 5)%500 == 0 ) {
            restart_r();
        }
//		pcout<<count <<" restart_r"<<std::endl;
		if (count%10 == 0 ) {
//			reconstruct_quad();
            post_process_force();
        }
        
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
		time += dt;
		for (unsigned int c = 0; c < n_locally_cells; ++c) {
			g_i = local_to_global_index_map[c];
			RHO_old(c) = RHO(g_i);
			RHO_U_old(c) = RHO_U(g_i);
			RHO_V_old(c) = RHO_V(g_i);
			RHO_W_old(c) = RHO_W(g_i);
			E_old(c) = E(g_i);
		}
       
		unsigned int NUM_STAGE = 7;
		double rho, rho_u, rho_v, rho_w, e;

		for(unsigned int stage = 1 ; stage < NUM_STAGE ; stage++) {
			compute_rhs();
			K1[stage-1] = rhs1;
			K2[stage-1] = rhs2;
			K3[stage-1] = rhs3;
			K4[stage-1] = rhs4;
			K5[stage-1] = rhs5;
			for (unsigned int c = 0; c < n_locally_cells; ++c) {
	   			g_i = local_to_global_index_map[c];
				rho  = RHO_old(c);	rho_u  = RHO_U_old(c);	rho_v  = RHO_V_old(c);
				rho_w  = RHO_W_old(c);	e  = E_old(c);
				for(unsigned int inner = 0 ; inner < stage ; inner++) {
					rho  += A_RK[stage][inner]*K1[inner](c)*dt;
					rho_u  += A_RK[stage][inner]*K2[inner](c)*dt;
					rho_v  += A_RK[stage][inner]*K3[inner](c)*dt;
					rho_w  += A_RK[stage][inner]*K4[inner](c)*dt;
					e  += A_RK[stage][inner]*K5[inner](c)*dt;
				}
				local_RHO(g_i) = rho;
				local_RHO_U(g_i) = rho_u;
				local_RHO_V(g_i) = rho_v;
				local_RHO_W(g_i) = rho_w;
				local_E(g_i) = e;
			}
			local_RHO.compress(VectorOperation::insert);
			local_RHO_U.compress(VectorOperation::insert);
			local_RHO_V.compress(VectorOperation::insert);
			local_RHO_W.compress(VectorOperation::insert);
			local_E.compress(VectorOperation::insert);

			RHO = local_RHO;
			RHO_U = local_RHO_U;
			RHO_V = local_RHO_V;
			RHO_W = local_RHO_W;
			E = local_E;
		}
      
		compute_rhs();
		K1[6] = rhs1;
		K2[6] = rhs2;
		K3[6] = rhs3;
		K4[6] = rhs4;
		K5[6] = rhs5;

		for (unsigned int c = 0; c < n_locally_cells; ++c) {
			g_i = local_to_global_index_map[c];
			rho  = RHO_old(c);	rho_u  = RHO_U_old(c);	rho_v  = RHO_V_old(c);
			rho_w  = RHO_W_old(c);	e  = E_old(c);
			for(unsigned int stage = 0 ; stage < NUM_STAGE ; stage++) {
				rho  += B_RK[stage]*K1[stage](c)*dt;
				rho_u  += B_RK[stage]*K2[stage](c)*dt;
				rho_v  += B_RK[stage]*K3[stage](c)*dt;
				rho_w  += B_RK[stage]*K4[stage](c)*dt;
				e  += B_RK[stage]*K5[stage](c)*dt;
			}
			local_RHO(g_i) = rho;
			local_RHO_U(g_i) = rho_u;
			local_RHO_V(g_i) = rho_v;
			local_RHO_W(g_i) = rho_w;
			local_E(g_i) = e;
        }

	    auto start_com = std::chrono::system_clock::now();

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_RHO_W.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		RHO_W = local_RHO_W;
		E = local_E;

	    auto end_com = std::chrono::system_clock::now();
	    std::chrono::duration<double> elapsed_seconds_com = end_com - start_com;
	    std::chrono::duration<double> elapsed_seconds_rk6 = end_com - start_rk6;

    	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0 && count == 2){
		    std::ofstream fout_convergence ;
    	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
    	    fout_convergence.precision(7) ;
	
    	    const std::string filename = "timer.dat";
		    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);
	
    		fout_convergence << "time taken by data transfer of rk6 = " << elapsed_seconds_com.count() << std::endl;
    		fout_convergence << "time taken by 1 step of rk6 = " << elapsed_seconds_rk6.count() << std::endl;
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

    	fout_convergence << "time taken by rk6 = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}

    output_results();
    restart(); 
} 

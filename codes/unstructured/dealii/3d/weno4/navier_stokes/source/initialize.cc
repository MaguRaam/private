#include "../include/Weno432.h"


// Initialize the solution 

void Weno4_3D::initialize() {

    unsigned int N_gp = 3;               // No. of quadrature points
    QGauss<3> quadrature_formula(N_gp);
//    FEValues<3> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);
    FEValues<3> fv_values (mapping, fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<3>::active_cell_iterator endc = dof_handler.end();

    Point<3> q_point;

    double V0;

	unsigned int g_i;

    Vector<double> U_(6); Vector<double> W_(6);
	double rho, rho_u, rho_v, rho_w, e, phi;
    
	if( RESTART ) {

		unsigned int status;

		std::ifstream status_S("restart/Status.dat", std::ios::in);
		if (!status_S) {
			std::cerr << "Error: Status of Restart files cannot be verified.\n";
			std::exit(EXIT_FAILURE);
		}
		status_S >> time;
		status_S >> status;
		status_S.close();

		if(status == 1 ) {

			double rho ;
   
		    std::ifstream fin_RHO("restart/RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 5) + ".dat");

	    	if ( !(fin_RHO.is_open()) ) {
	    	    std::cerr << "Restart files could not be opened" << std::endl; 
	    	    std::exit(EXIT_FAILURE);
		    }

			fin_RHO >> time;

			pcout << "Restarting from:" << time << std::endl;

			for (unsigned int c = 0; c < n_locally_cells; ++c) {

				g_i = local_to_global_index_map[c];

				fin_RHO >> rho; 		
				fin_RHO >> rho_u;
				fin_RHO >> rho_v;
				fin_RHO >> rho_w;
				fin_RHO >> e;
				
			
				local_RHO(g_i) = rho;
			    local_RHO_U(g_i) = rho_u;
			    local_RHO_V(g_i) = rho_v;
			    local_RHO_W(g_i) = rho_w;
				local_E(g_i) = e;
				

			}
		}
		else {

   			pcout << " The Restart files are corrupted; Checking Redundant files " << std::endl;

			std::ifstream status_S("restart_r/Status.dat", std::ios::in);
			if (!status_S) {
				std::cerr << "Error: Status of Restart_r files cannot be verified.\n";
				std::exit(EXIT_FAILURE);
			}
			status_S >> time;
			status_S >> status;
			status_S.close();

			if(status == 1 ) {
	
				double rho;
	   
			    std::ifstream fin_RHO("restart_r/RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 5) + ".dat");

		    	if ( !(fin_RHO.is_open()) ) {
		    	    std::cerr << "Unable to open restart_r folder" << std::endl; 
		    	    std::exit(EXIT_FAILURE);
			    }

				fin_RHO >> time;

				pcout << "Restarting from:" << time << std::endl;
	
				for (unsigned int c = 0; c < n_locally_cells; ++c) {

					g_i = local_to_global_index_map[c];

					fin_RHO >> rho; 		
					fin_RHO >> rho_u;
					fin_RHO >> rho_v;
					fin_RHO >> rho_w;
					fin_RHO >> e;
								

					local_RHO(g_i) = rho;	
				    local_RHO_U(g_i) = rho_u;
				    local_RHO_V(g_i) = rho_v;
				    local_RHO_W(g_i) = rho_w;
					local_E(g_i) = e;
					
				}
//				copy_data_initialize();
			}
			else {
        			std::cerr << "Error: The Restart Files and Redundant files are corrupted; " << "Simulation Cannot be Started."<<std::endl;
		    	    std::exit(EXIT_FAILURE);
			}

		}

	}
	else{

	    Point<3> q_point;
		cell = dof_handler.begin_active();
	    for (unsigned int c = 0; cell != endc; ++cell, ++c) {

		if (cell->is_locally_owned()){	

	        cell->get_dof_indices(local_dof_indices);						
			g_i = local_dof_indices[0] ;

	        V0 = 0.0; 
			fv_values.reinit(cell);
			for (unsigned int i=0; i<fv_values.n_quadrature_points; ++i)
				V0 += fv_values.JxW (i);
				
			rho = 0.0, rho_u = 0.0, rho_v = 0.0, rho_w = 0.0, e = 0.0, phi = 0.0;

			for (unsigned int i = 0; i < fv_values.n_quadrature_points; i++) {

				q_point = fv_values.quadrature_point(i);

				W_ = initial_condition(q_point, h_min); 
				claw.primitive_to_conserved(W_, U_); 

				rho   += (1./V0)*fv_values.JxW (i)*U_(0);
				rho_u += (1./V0)*fv_values.JxW (i)*U_(1);
				rho_v += (1./V0)*fv_values.JxW (i)*U_(2);
				rho_w += (1./V0)*fv_values.JxW (i)*U_(3);
				e     += (1./V0)*fv_values.JxW (i)*U_(4);
			}


			local_RHO(g_i) = rho;; 		
		    local_RHO_U(g_i) = rho_u;;
		    local_RHO_V(g_i) = rho_v;
		    local_RHO_W(g_i) = rho_w;
			local_E(g_i) = e;
			

   	    }
		}
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

//	std::cout<<"rank: "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<" RHO: "<<RHO(465653)<<std::endl;

    	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
		    std::ofstream fout_convergence ;
    	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
    	    fout_convergence.precision(7) ;
	
    	    const std::string filename = "log.dat";
		    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    		fout_convergence << "h_min: " <<h_min<< std::endl;  
    	    fout_convergence.close();
		}

} 

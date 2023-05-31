#include "../include/Weno432.h"


// Initialize the solution 

void Weno4_2D::initialize() {

    unsigned int N_gp = 4;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<2>::active_cell_iterator cell;

    Point<2> q_point;

    double V0;

	double h_old,h_min_local; 
	h_min_local = 1e6; 

	unsigned int g_i;

    Vector<double> U(4); Vector<double> W(4);
    
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

			double rho, rho_u, rho_v, e ;
   
		    std::ifstream fin_RHO("restart/RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 5) + ".dat");

	    	if ( !(fin_RHO.is_open()) ) {
	    	    std::cerr << "Restart files could not be opened" << std::endl; 
	    	    std::exit(EXIT_FAILURE);
		    }

			fin_RHO >> time;

			pcout << "Restarting from:" << time << std::endl;

			for (unsigned int c = 0; c < n_locally_cells; ++c) {

				g_i = local_to_global_index_map[c];
				cell = local_index_to_iterator[c];			

				V0 = cell->measure(); 
		
				h_old = std::sqrt(V0);

				if (h_old < h_min_local) {
					h_min_local = h_old; 
				}
		
				fin_RHO >> rho; 	
    		    fin_RHO >> rho_u;
    		    fin_RHO >> rho_v;
    		    fin_RHO >> e;
	
				local_RHO(g_i) = rho;; 	
    		    local_RHO_U(g_i) = rho_u;;
    		    local_RHO_V(g_i) = rho_v;
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
	
				double rho, rho_u, rho_v, e ;
	   
			    std::ifstream fin_RHO("restart_r/RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 5) + ".dat");

		    	if ( !(fin_RHO.is_open()) ) {
		    	    std::cerr << "Unable to open restart_r folder" << std::endl; 
		    	    std::exit(EXIT_FAILURE);
			    }

				fin_RHO >> time;

				pcout << "Restarting from:" << time << std::endl;
	
				for (unsigned int c = 0; c < n_locally_cells; ++c) {

					g_i = local_to_global_index_map[c];
					cell = local_index_to_iterator[c];			
	
					V0 = cell->measure(); 
		
					h_old = std::sqrt(V0);
		
					if (h_old < h_min_local) {
						h_min_local = h_old; 
					}
		
					fin_RHO >> rho; 	
	    		    fin_RHO >> rho_u;
	    		    fin_RHO >> rho_v;
	    		    fin_RHO >> e;
		
					local_RHO(g_i) = rho;; 	
	    		    local_RHO_U(g_i) = rho_u;;
	    		    local_RHO_V(g_i) = rho_v;
	    		    local_E(g_i) = e;

				}

			}
			else {
        			std::cerr << "Error: The Restart Files and Redundant files are corrupted; " << "Simulation Cannot be Started."<<std::endl;
		    	    std::exit(EXIT_FAILURE);
			}

		}

	}
	else{

    	Vector<double> U(4); Vector<double> W(4);
	
	    Point<2> q_point;

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];
			cell = local_index_to_iterator[c]; 

			fv_values.reinit(cell);

			V0 = cell->measure(); 
		
			h_old = std::sqrt(V0);
	
			if (h_old < h_min_local) {
				h_min_local = h_old; 
			}
	
			double rho = 0.0, rho_u = 0.0, rho_v = 0.0, e = 0.0;

			for (unsigned int i = 0; i < N_gp*N_gp; i++) {

				q_point = fv_values.quadrature_point(i);

				W = initial_condition(q_point); 
				U = primitive_to_conserved(W); 
	
				rho   += (1./V0)*fv_values.JxW (i)*U(0);
				rho_u += (1./V0)*fv_values.JxW (i)*U(1);
				rho_v += (1./V0)*fv_values.JxW (i)*U(2);
				e     += (1./V0)*fv_values.JxW (i)*U(3);
			}

			local_RHO(g_i)   = rho;
			local_RHO_U(g_i) = rho_u;
			local_RHO_V(g_i) = rho_v;
			local_E(g_i)     = e;

   	    }

	}

	local_RHO.compress(VectorOperation::insert);
	local_RHO_U.compress(VectorOperation::insert);
	local_RHO_V.compress(VectorOperation::insert);
	local_E.compress(VectorOperation::insert);

	RHO = local_RHO;
	RHO_U = local_RHO_U;
	RHO_V = local_RHO_V;
	E = local_E;

//	std::cout<<"rank: "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<" RHO: "<<RHO(465653)<<std::endl;

	h_min = Utilities::MPI::min (h_min_local, MPI_COMM_WORLD);
    pcout << "h_min: " <<h_min<< std::endl;  

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

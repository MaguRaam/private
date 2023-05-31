#include "../include/Weno32.h"

void Weno3_2D::restart() {
   
    std::ofstream fout_RHO("restart/RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
	std::ofstream fout_RHO_U("restart/RHO_U_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
	std::ofstream fout_RHO_V("restart/RHO_V_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
	std::ofstream fout_E("restart/E_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
	std::ofstream fout_status("restart/status_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
    
    if ( !(fout_RHO.is_open()) ) {
        std::cerr << "Unable to open restart folder" << std::endl; 
        std::exit(EXIT_FAILURE);
    }
    
    fout_RHO.flags( std::ios::dec | std::ios::scientific );
	fout_RHO.precision(16);
    
    fout_RHO_U.flags( std::ios::dec | std::ios::scientific );
	fout_RHO_U.precision(16);
    
    fout_RHO_V.flags( std::ios::dec | std::ios::scientific );
	fout_RHO_V.precision(16);
    
    fout_E.flags( std::ios::dec | std::ios::scientific );
	fout_E.precision(16);

    fout_status.flags( std::ios::dec | std::ios::scientific );
	fout_status.precision(16);
    
    fout_status << time << std::endl; 
    
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

	for (; cell != endc; ++cell) {
	    if (cell->is_locally_owned()){
	        cell->get_dof_indices(local_dof_indices);

	        fout_RHO << local_RHO(local_dof_indices[0] ) << std::endl; 
	        fout_RHO_U << local_RHO_U(local_dof_indices[0] ) << std::endl;
	        fout_RHO_V << local_RHO_V(local_dof_indices[0] ) << std::endl;
	        fout_E << local_E(local_dof_indices[0] ) << std::endl;

		}
	}
    
    fout_RHO.close(); 
    fout_RHO_U.close();
    fout_RHO_V.close();
    fout_E.close(); 
    
}

void Weno3_2D::restart_r() {
   
    std::ofstream fout_RHO("restart_r/RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
	std::ofstream fout_RHO_U("restart_r/RHO_U_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
	std::ofstream fout_RHO_V("restart_r/RHO_V_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
	std::ofstream fout_E("restart_r/E_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
	std::ofstream fout_status("restart_r/status_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
    
    if ( !(fout_RHO.is_open()) ) {
        std::cerr << "Unable to open restart_r folder" << std::endl; 
        std::exit(EXIT_FAILURE);
    }
    
    fout_RHO.flags( std::ios::dec | std::ios::scientific );
	fout_RHO.precision(16);
    
    fout_RHO_U.flags( std::ios::dec | std::ios::scientific );
	fout_RHO_U.precision(16);
    
    fout_RHO_V.flags( std::ios::dec | std::ios::scientific );
	fout_RHO_V.precision(16);
    
    fout_E.flags( std::ios::dec | std::ios::scientific );
	fout_E.precision(16);

    fout_status.flags( std::ios::dec | std::ios::scientific );
	fout_status.precision(16);
    
    fout_status << time << std::endl; 
    
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

	for (; cell != endc; ++cell) {
	    if (cell->is_locally_owned()){
	        cell->get_dof_indices(local_dof_indices);

	        fout_RHO << local_RHO(local_dof_indices[0] ) << std::endl; 
	        fout_RHO_U << local_RHO_U(local_dof_indices[0] ) << std::endl;
	        fout_RHO_V << local_RHO_V(local_dof_indices[0] ) << std::endl;
	        fout_E << local_E(local_dof_indices[0] ) << std::endl;

		}
	}
    
    fout_RHO.close(); 
    fout_RHO_U.close();
    fout_RHO_V.close();
    fout_E.close(); 
    
}


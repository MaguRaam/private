#include "../include/Weno32.h"


// Put everything together - solve the actual problem  

void Weno3_2D::run() {
    
    pcout << "Running with "
#ifdef USE_PETSC_LA
          << "PETSc"
#else
          << "Trilinos"
#endif
          << " on "
          << Utilities::MPI::n_mpi_processes(mpi_communicator)
          << " MPI rank(s)..." << std::endl;      

    make_grid(); 

	pcout << "Number of active cells: "
						<< triangulation.n_global_active_cells()
                        << std::endl;    
    pcout << "============================" << std::endl;

	setup_system();
	compute_cell_properties(); 
    compute_weno_polynomial_constants(); 
	initialize();
	precompute_matrices();
    compute_IS_constants(); 
    solve_ssprk33();   
	
	MPI_Barrier(MPI_COMM_WORLD); 
	
	plot_tecplot(); 
	
	MPI_Barrier(MPI_COMM_WORLD);
	
    computing_timer.print_summary ();

} 

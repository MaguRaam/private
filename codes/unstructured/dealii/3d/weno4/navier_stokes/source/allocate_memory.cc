#include "../include/Weno432.h"

// Setup the system - allocate the necessary memory 

void Weno4_3D::allocate_memory() {

//    auto start = std::chrono::system_clock::now();

    pcout << "allocate memory" << std::endl;

	RHO.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	RHO_U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	RHO_V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	RHO_W.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	E.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);


	local_RHO.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	local_RHO_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	local_RHO_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	local_RHO_W.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	local_E.reinit(locally_owned_dofs, MPI_COMM_WORLD);

	error_RHO.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	local_error_RHO.reinit(locally_owned_dofs, MPI_COMM_WORLD);

	error_RHO_U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	local_error_RHO_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);

	error_RHO_V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	local_error_RHO_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);

	error_RHO_W.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	local_error_RHO_W.reinit(locally_owned_dofs, MPI_COMM_WORLD);

	error_E.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	local_error_E.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    coeffs_RHO.resize(n_relevant_cells);
    coeffs_RHO_U.resize(n_relevant_cells);
    coeffs_RHO_V.resize(n_relevant_cells);
    coeffs_RHO_W.resize(n_relevant_cells);
    coeffs_E.resize(n_relevant_cells);

	upwind_coeffs_RHO.resize(n_relevant_cells);
    upwind_coeffs_RHO_U.resize(n_relevant_cells);
    upwind_coeffs_RHO_V.resize(n_relevant_cells);
    upwind_coeffs_RHO_W.resize(n_relevant_cells);
    upwind_coeffs_E.resize(n_relevant_cells);
	
	Cell.resize(n_store_cell);
	is_1st_order.resize(n_relevant_cells);    
    WENO_poly_consts.resize(n_relevant_cells); 

    CLS_R4.resize(n_relevant_cells);
    CLS_R3.resize(n_relevant_cells);
    
    is_admissible_R3_d.resize(n_relevant_cells);
    CLS_R3_d.resize(n_relevant_cells);

    LU_R2_d.resize(n_relevant_cells);

    CLS_R4_slip.resize(n_relevant_cells);
    CLS_R3_slip.resize(n_relevant_cells);
    CLS_R3_d_slip.resize(n_relevant_cells);
    LU_R2_d_slip.resize(n_relevant_cells);

	rhs1.reinit(n_locally_cells);
	rhs2.reinit(n_locally_cells);
	rhs3.reinit(n_locally_cells);
	rhs4.reinit(n_locally_cells);
	rhs5.reinit(n_locally_cells);

    for (unsigned int i = 0; i < n_relevant_cells; i++) {
        is_admissible_R3_d[i].resize(6);
        CLS_R3_d[i].resize(6);
        CLS_R3_d_slip[i].resize(6);
        LU_R2_d[i].resize(8);
        LU_R2_d_slip[i].resize(8);
        coeffs_RHO[i].reinit(20);
        coeffs_RHO_U[i].reinit(20);
        coeffs_RHO_V[i].reinit(20);
        coeffs_RHO_W[i].reinit(20);
        coeffs_E[i].reinit(20);
        upwind_coeffs_RHO[i].reinit(20);
        upwind_coeffs_RHO_U[i].reinit(20);
        upwind_coeffs_RHO_V[i].reinit(20);
        upwind_coeffs_RHO_W[i].reinit(20);
        upwind_coeffs_E[i].reinit(20);
        WENO_poly_consts[i].reinit(19); 
    }

	Utilities::System::MemoryStats stat;
	Utilities::System::get_memory_stats(stat);
	pcout<<stat.VmRSS/std::pow(2,20)<<std::endl;
	pcout<<"total (on 24 core) : "<<24.0*stat.VmRSS/std::pow(2,20)<<std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "log.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Memory consumption in Allocate memory per node in GB = " << stat.VmRSS/std::pow(2,20) << std::endl;
        fout_convergence.close();
	}
/*
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "timer.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Time taken to allocate_memory = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}
*/
}

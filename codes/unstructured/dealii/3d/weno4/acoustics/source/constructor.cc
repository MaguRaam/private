#include "../include/Weno432.h"


// Construtor for the WENO4 class 

Weno4_3D::Weno4_3D (double ft, double cfl_no, bool restart, unsigned int refinement, double weight)
  :
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (MPI_COMM_WORLD, typename Triangulation<3>::MeshSmoothing(Triangulation<3>::none), false, parallel::shared::Triangulation<3>::Settings::partition_zorder),
	//triangulation (MPI_COMM_WORLD),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times),
	dt(0.0), 
    finalTime(ft),
	cfl (cfl_no),
	weight_ls(weight),
	RESTART(restart),
	n_refinement(refinement),
    fv (0),
    dof_handler (triangulation),
    claw(1.4, 1.0, 1.0),
	mapping(1)
//	true_solution()
{}

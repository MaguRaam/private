#include "../include/Weno432.h"


// Construtor for the WENO4 class 

Weno4_2D::Weno4_2D (double ft, double cfl_no, bool restart, unsigned int refinement)
  :
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (MPI_COMM_WORLD, typename Triangulation<2>::MeshSmoothing(Triangulation<2>::none), false, parallel::shared::Triangulation<2>::Settings::partition_zorder),
	//triangulation (MPI_COMM_WORLD),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times),
	dt(0.0), 
    finalTime(ft),
	cfl (cfl_no),
	RESTART(restart),
	n_refinement(refinement),
    fv (0),
    dof_handler (triangulation)
{}

#include "../include/Weno32.h"


// Construtor for the WENO4 class 

Weno3_2D::Weno3_2D (double ft, double cfl_no, bool restart)
  :
    mpi_communicator (MPI_COMM_WORLD),
    triangulation (mpi_communicator,
                   typename Triangulation<2>::MeshSmoothing
                   (Triangulation<2>::smoothing_on_refinement |
                    Triangulation<2>::smoothing_on_coarsening)),
    pcout (std::cout, (Utilities::MPI::this_mpi_process(mpi_communicator) == 1)),
    computing_timer (mpi_communicator, pcout, TimerOutput::summary, TimerOutput::wall_times),
	dt(0.0), 
    finalTime(ft),
	cfl (cfl_no),
	RESTART(restart),
    fv (0),
    dof_handler (triangulation)
{}

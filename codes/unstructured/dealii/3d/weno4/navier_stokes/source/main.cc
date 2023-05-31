#include "../include/Weno432.h"


// Main function for the problem 

int main(int argc, char *argv[]){

	try
    {
		std::cout.flags( std::ios::dec | std::ios::scientific ) ; 
		std::cout.precision(16) ;

		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
		double finalTime = 100.0;
		double cfl = 0.25;
		bool restart = false;
		unsigned int refinement = 1;	// 0 for 512x512 and so on ...
		double weight = 1.0e7;

		Weno4_3D test_problem(finalTime, cfl, restart, refinement, weight);
		test_problem.run();
    }
  
  catch (std::exception &exc)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Exception on processing: " << std::endl
                << exc.what() << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;

      return 1;
    }
  catch (...)
    {
      std::cerr << std::endl << std::endl
                << "----------------------------------------------------"
                << std::endl;
      std::cerr << "Unknown exception!" << std::endl
                << "Aborting!" << std::endl
                << "----------------------------------------------------"
                << std::endl;
      return 1;
    }

    return 0;
} 

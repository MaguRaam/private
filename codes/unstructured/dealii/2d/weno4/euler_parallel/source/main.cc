#include "../include/Weno432.h"


// Main function for the problem 

int main(int argc, char *argv[]){

	try
    {
		std::cout.flags( std::ios::dec | std::ios::scientific ) ; 
		std::cout.precision(7) ;

		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
		
		double finalTime = 0.3;
		double cfl = 0.4;
		bool restart = false;
		unsigned int refinement = 1;	// 0 for 262k (512 x 512) , 1 for 1048k, 2 for 4194k amnd so on ..

		Weno4_2D test_problem(finalTime, cfl, restart, refinement);
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

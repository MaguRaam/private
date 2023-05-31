#include "../include/Weno32.h"


// Main function for the problem 
int main(int argc, char *argv[]){

	{

		std::cout.flags( std::ios::dec | std::ios::scientific ) ; 
		std::cout.precision(6) ;

		Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1); 
		
		double finalTime = 1.95;
		double cfl = 0.4;
		bool restart = true;
		Weno3_2D test_problem(finalTime, cfl, restart);
		test_problem.run();

	}
     
    return 0; 
}

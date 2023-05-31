#include "../include/Weno32.h"


// Main function for the problem 

int main() {

	double finalTime = 0.3;
	double cfl = 0.4;
	Weno3_2D test_problem(finalTime, cfl);
	test_problem.run();

    return 0;
} 

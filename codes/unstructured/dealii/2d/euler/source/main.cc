#include "../include/Euler.h"


// Main function for the problem 

int main() {

	double finalTime = 0.2;
	double cfl = 0.3;
	Euler_2D test_problem(finalTime, cfl);
	test_problem.run();

    return 0;
} 

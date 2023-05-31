#include "../include/Euler.h"


// Put everything together - solve the actual problem  

void Euler_2D::run() {
    
    auto start = std::chrono::system_clock::now();
	
	make_grid(); 
	setup_system();
	initialize();
	solve();   
    
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "Time taken = " << elapsed_seconds.count() << std::endl; 
} 

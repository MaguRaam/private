#include "../include/Weno32.h"


// Put everything together - solve the actual problem  

void Weno3_2D::run() {
    
    auto start = std::chrono::system_clock::now();
        

        
    make_grid(); 
	setup_system();
	select_stencils(); 
	compute_cell_properties(); 
    compute_weno_polynomial_constants(); 
	initialize();
	precompute_matrices();
    //precompute_matrices_veclocity();
    compute_IS_constants(); 
    solve_ader();   
	
	
    auto end = std::chrono::system_clock::now();
	std::chrono::duration<double> elapsed_seconds = end-start;

    std::cout << "Time taken = " << elapsed_seconds.count() << std::endl;
} 

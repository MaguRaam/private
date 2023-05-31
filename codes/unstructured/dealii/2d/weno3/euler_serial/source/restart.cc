#include "../include/Weno32.h"

void Weno3_2D::restart(double time, unsigned int count) {
    
    std::ofstream fout_RHO("restart/RHO.dat");
	std::ofstream fout_U("restart/RHO_U.dat");
	std::ofstream fout_V("restart/RHO_V.dat");
	std::ofstream fout_P("restart/E.dat");
    
    if ( !(fout_RHO.is_open()) ) {
        std::cerr << "Unable to open restart folder" << std::endl; 
        std::exit(EXIT_FAILURE);
    }
    
    fout_RHO.flags( std::ios::dec | std::ios::scientific );
	fout_RHO.precision(10);
    
    fout_U.flags( std::ios::dec | std::ios::scientific );
	fout_U.precision(10);
    
    fout_V.flags( std::ios::dec | std::ios::scientific );
	fout_V.precision(10);
    
    fout_P.flags( std::ios::dec | std::ios::scientific );
	fout_P.precision(10);
    
    fout_RHO << "time  = " << time << std::endl; 
    fout_RHO << "count  = " << count << std::endl;
    fout_U << "time  = " << time << std::endl;
    fout_U << "count  = " << count << std::endl;
    fout_V << "time  = " << time << std::endl;
    fout_V << "count  = " << count << std::endl;
    fout_P << "time  = " << time << std::endl;
    fout_P << "count  = " << count << std::endl;
    
    
    for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
        fout_RHO << RHO(c) << std::endl; 
        fout_U << RHO_U(c) << std::endl;
        fout_V << RHO_V(c) << std::endl;
        fout_P << E(c) << std::endl;
    }
    
    fout_RHO.close(); 
    fout_U.close();
    fout_V.close();
    fout_P.close(); 
    
}

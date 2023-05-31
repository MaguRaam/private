#include "../include/LU.h"
#include "../include/CLS.h"
#include "../include/Headers.h"
#include "../include/Exceptions.h"
#include "../include/Riemann_Solvers.h" 
#include "../include/Weno32.h"


// Make the grid  

void Weno3_2D::make_grid() {

    std::cout << "===========================" << std::endl;
    std::cout << "Making grid" << std::endl;

    std::vector<unsigned int> repetions(2); // No. of cells in x and y directions 
    repetions[0] = 100;
    repetions[1] = 100; 
    
    // Diagonal points of the domain 
    Point<2> P1(0.0,  0.0);
    Point<2> P2(1.0,  1.0);
    
    bool colorize = true;  // Set boundary ids for the four boundaries 
    
    // Left face: 0 
    // Right face: 1
    // Bottom face: 2
    // Top face: 3
    
    GridGenerator::subdivided_hyper_rectangle(triangulation, repetions, P1, P2, colorize);
    GridTools::distort_random (0.15, triangulation, true);
    
    
    std::cout << "Number of cells: "
                << triangulation.n_active_cells()
                << std::endl;
    std::cout << "Done!" << std::endl;  
    std::cout << "===========================" << std::endl; 
    
    
}

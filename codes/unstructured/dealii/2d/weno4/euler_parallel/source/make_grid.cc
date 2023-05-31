#include "../include/Weno432.h"

void Weno4_2D::make_grid() {

	std::vector<unsigned int> repetions(2); // No. of cells in x and y directions 
    repetions[0] = 64; 
    repetions[1] = 64; 
    
    // Diagonal points of the domain 
    Point<2> P1( 0.0, 0.0);
    Point<2> P2( 1.0, 1.0);
    
    bool colorize = true;  // Set boundary ids for the four boundaries 
    
    // Left face: 0 
    // Right face: 1
    // Bottom face: 2
    // Top face: 3
    
    GridGenerator::subdivided_hyper_rectangle(triangulation, repetions, P1, P2, colorize);
    GridTools::distort_random (0.15, triangulation, true);

}

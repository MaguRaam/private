#include "../include/Weno32.h"


void Weno3_2D::post_process(unsigned int i) {
    
    // Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
    
    const std::string filename = "contours/RHO_contour_" + Utilities::int_to_string (i, 4) + ".dat";
    
    std::ofstream fout_RHO(filename);
	
	if ( !(fout_RHO.is_open()) ) {
        std::cerr << "Unable to open contours folder" << std::endl; 
        std::exit(EXIT_FAILURE);
    } 
    
    Point<2> V; 
	Point<2> C;
	
	unsigned int V_index, V_neighbor_index; 
	double value; 
	
	Vert2Cell<2> V2C(triangulation);
    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
		
		for (unsigned int v = 0; v < 4; v++) {
			
			value = 0.0;
			V = cell->vertex(v);
			V_index = cell->vertex_index(v); 
			
			for (unsigned int i = 0; i < V2C.n_cells_sharing_vertex(V_index) ; ++i) {
				V_neighbor_index = V2C.cell_sharing_vertex(V_index, i);
				value += RHO(V_neighbor_index); 
			}
			
			value = value/(V2C.n_cells_sharing_vertex(V_index)); 
		
			fout_RHO << value << "\n"; 
			
		}
        
    }
    
    fout_RHO.close(); 
}
 

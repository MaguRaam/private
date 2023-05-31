#include "../include/Weno32.h"

// select stencil for each cell 

void Weno3_2D::select_stencils() {
	
	std::cout << "Selecting stencils for each cell" << "\n"; 
	
	// --------------- Now loop over all the cells to fill the stencil selector structure --------------- //
	
	Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();
	Triangulation<2>::active_cell_iterator endc = triangulation.end();
	
	unsigned int v_NE, v_SW, v_NW, v_SE;
	int vertex_neighbor_index; 
	
	Vert2Cell<2> V2C(triangulation);
	
	for (int c = 0; cell != endc; ++cell, ++c) {
		
		// First fill all the face neighbor information 
	
		if (cell->neighbor_index(0) != -1) {
			Stencil[c].is_W_present = true; 
			Stencil[c].W_index = cell->neighbor_index(0); 
		}
		
		else {
			Stencil[c].is_W_present = false; 
			Stencil[c].W_index = -1;   
		}
	
		if (cell->neighbor_index(1) != -1) {
			Stencil[c].is_E_present = true; 
			Stencil[c].E_index = cell->neighbor_index(1); 
		}
		
		else {
			Stencil[c].is_E_present = false; 
			Stencil[c].E_index = -1;  
		}
	
		if (cell->neighbor_index(2) != -1) {
			Stencil[c].is_S_present = true; 
			Stencil[c].S_index = cell->neighbor_index(2); 
		}
		
		else {
			Stencil[c].is_S_present = false; 
			Stencil[c].S_index = -1;  
		}
	
		if (cell->neighbor_index(3) != -1) {
			Stencil[c].is_N_present = true; 
			Stencil[c].N_index = cell->neighbor_index(3); 
		}
		
		else {
			Stencil[c].is_N_present = false; 
			Stencil[c].N_index = -1;   
		}
		
		// Now come to the vertex neighbors 
		
		// SW Stencil 
		
		v_SW = cell->vertex_index(0); 
		
		// get all the cells sharing the given vertex 
		
		for (unsigned int i = 0; i < V2C.n_cells_sharing_vertex(v_SW) ; ++i) {
			
			vertex_neighbor_index = V2C.cell_sharing_vertex(v_SW, i); 
			
			if (vertex_neighbor_index != Stencil[c].S_index && vertex_neighbor_index != Stencil[c].W_index &&  vertex_neighbor_index != c) {
				Stencil[c].SW_index.push_back(vertex_neighbor_index); 
			}
		}
		
		// SE Stencil 
		
		v_SE = cell->vertex_index(1); 
		
		// get all the cells sharing the given vertex 
		
		for (unsigned int i = 0; i < V2C.n_cells_sharing_vertex(v_SE) ; ++i) {
			
			vertex_neighbor_index = V2C.cell_sharing_vertex(v_SE, i); 
			
			if (vertex_neighbor_index != Stencil[c].S_index && vertex_neighbor_index != Stencil[c].E_index &&  vertex_neighbor_index != c) {
				Stencil[c].SE_index.push_back(vertex_neighbor_index); 
			}
		}
		
		// NW Stencil 
		
		v_NW = cell->vertex_index(2); 
		
		// get all the cells sharing the given vertex 
		
		for (unsigned int i = 0; i < V2C.n_cells_sharing_vertex(v_NW) ; ++i) {
			
			vertex_neighbor_index = V2C.cell_sharing_vertex(v_NW, i); 
			
			if (vertex_neighbor_index != Stencil[c].N_index && vertex_neighbor_index != Stencil[c].W_index &&  vertex_neighbor_index != c) {
				Stencil[c].NW_index.push_back(vertex_neighbor_index); 
			}
		}
		
		// NE Stencil 
		
		v_NE = cell->vertex_index(3); 
		
		// get all the cells sharing the given vertex 
		
		for (unsigned int i = 0; i < V2C.n_cells_sharing_vertex(v_NE) ; ++i) {
			
			vertex_neighbor_index = V2C.cell_sharing_vertex(v_NE, i); 
			
			if (vertex_neighbor_index != Stencil[c].N_index && vertex_neighbor_index != Stencil[c].E_index &&  vertex_neighbor_index != c) {
				Stencil[c].NE_index.push_back(vertex_neighbor_index); 
			}
		}
		
	}
	
	std::cout << "Done!" << std::endl;
    std::cout << "===========================" << std::endl;
}

/*
 * vertex_to_cell.h
 *      Author: sunder
 */

#ifndef VERTEX_TO_CELL_H_
#define VERTEX_TO_CELL_H_ 


#include "Headers.h"

// Class for getting all the cells surrounding a vertex
// Ref. : Applied CFD Techniques: An Introduction based on Finite Element Methods, Rainald Lohner (Wiley 2008) 

// Forward declaration 

template <unsigned int>
class Vert2Cell; 

// Class function definitions 

template <unsigned int dim>
class Vert2Cell {
	unsigned int n_vertices; 
	unsigned int n_cells; 
	int* vert_to_cell1; // These two are linked lists 
	int* vert_to_cell2;
public:
	Vert2Cell(); 
	Vert2Cell(const Triangulation<dim>&); 
	~Vert2Cell(); 
	
	void initialzie(const Triangulation<dim>&); 
	
	unsigned int n_cells_sharing_vertex(unsigned int) const;
	unsigned int cell_sharing_vertex(unsigned int, unsigned int) const;  
}; 

#endif /* VERTEX_TO_CELL_H_ */

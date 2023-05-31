#include "../include/Weno32.h"


// Setup the system - allocate the necessary memory 

void Weno3_2D::setup_system() {

    TimerOutput::Scope t(computing_timer, "Setting up system");

    pcout << "Setting up system" << std::endl;

	time = 0.0;

	dof_handler.distribute_dofs (fv);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
    DoFTools::extract_locally_relevant_dofs (dof_handler,
                                             locally_relevant_dofs);

	dofs_per_cell = fv.dofs_per_cell;
	local_dof_indices.resize(dofs_per_cell);
	local_neighbor_dof_indices.resize(dofs_per_cell);

	n_locally_cells = dof_handler.n_locally_owned_dofs();

	std::cout<<"no of cells on processor: "<<n_locally_cells<<std::endl;

	RHO.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	RHO_U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	RHO_V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
	E.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

	local_RHO.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	local_RHO_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	local_RHO_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);
	local_E.reinit(locally_owned_dofs, MPI_COMM_WORLD);

	rhs1.reinit(n_locally_cells);
	rhs2.reinit(n_locally_cells);
	rhs3.reinit(n_locally_cells);
	rhs4.reinit(n_locally_cells);

    coeffs_x_RHO.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_y_RHO.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_xx_RHO.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_yy_RHO.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_xy_RHO.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

    coeffs_x_RHO_U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_y_RHO_U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_xx_RHO_U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_yy_RHO_U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_xy_RHO_U.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

    coeffs_x_RHO_V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_y_RHO_V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_xx_RHO_V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_yy_RHO_V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_xy_RHO_V.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

    coeffs_x_E.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_y_E.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_xx_E.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_yy_E.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    coeffs_xy_E.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

    local_coeffs_x_RHO.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_y_RHO.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_xx_RHO.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_yy_RHO.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_xy_RHO.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    local_coeffs_x_RHO_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_y_RHO_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_xx_RHO_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_yy_RHO_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_xy_RHO_U.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    local_coeffs_x_RHO_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_y_RHO_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_xx_RHO_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_yy_RHO_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_xy_RHO_V.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    local_coeffs_x_E.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_y_E.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_xx_E.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_yy_E.reinit(locally_owned_dofs, MPI_COMM_WORLD);
    local_coeffs_xy_E.reinit(locally_owned_dofs, MPI_COMM_WORLD);

    WENO_poly_consts_x.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    WENO_poly_consts_y.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    WENO_poly_consts_xx.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    WENO_poly_consts_yy.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);
    WENO_poly_consts_xy.reinit(locally_owned_dofs, locally_relevant_dofs, MPI_COMM_WORLD);

    IS_constants.resize(n_locally_cells);
	Cell.resize(n_locally_cells);
    
    CLS_R3.resize(n_locally_cells);
    LU_R21.resize(n_locally_cells);
    LU_R22.resize(n_locally_cells);
    LU_R23.resize(n_locally_cells);
    LU_R24.resize(n_locally_cells);
    
    is_corner_cell.resize(n_locally_cells);
    
    
    for (unsigned int i = 0; i < n_locally_cells; i++) {
        IS_constants[i].reinit(6);
    }
    
	// Initialize the wall boundary matrices  
    
    // Get the number of cells at the wall boundary 
 /*   
    unsigned int n_wall_boundary_cells = 0; 
    
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        if (cell->at_boundary()) {
            
            for (unsigned int f=0; f < GeometryInfo<2>::faces_per_cell; ++f) {
                
                if (cell->face(f)->at_boundary()) {
                    
                    if (cell->face(f)->boundary_id() == 2) { // boundary_id 2 corresponds to the reflecting wall 
                        
                        wall_boundary_global_index_map[c] = n_wall_boundary_cells; 
                        
                        n_wall_boundary_cells ++; 
                        
                    }
                    
                }
            }
        }
    }
            
    std::cout << "Number of wall boundary cells = " <<  n_wall_boundary_cells << std::endl; 

    
    CLS_R3_slip.resize(n_wall_boundary_cells);
    CLS_R2_slip.resize(n_wall_boundary_cells);
*/
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();


    unsigned int local_index = 0; 

	std::set<unsigned int> vertex_index;	
    std::map <unsigned int, unsigned int> vertex_map;
    
    for (; cell != endc; ++cell) {
        
		if (cell->is_locally_owned()){

	        cell->get_dof_indices(local_dof_indices);
                          
			global_to_local_index_map[local_dof_indices[0] ] = local_index; 
//			local_to_global_index_map[local_index] = local_dof_indices[0]; 

	        for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i){
				vertex_index.insert(cell->vertex_index(i));		
		}

		if (cell->is_ghost()){

	        cell->get_dof_indices(local_dof_indices);
                          
			global_to_local_index_map[local_dof_indices[0] ] = local_index; 
//			local_to_global_index_map[local_index] = local_dof_indices[0]; 

	        for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i){
				vertex_index.insert(cell->vertex_index(i));		
			}
		}
		local_index++;                                             
        }
	}

	n_vertices = vertex_index.size();

	std::vector< std::set<DoFHandler<2>::active_cell_iterator> > vertex_to_cell_iterator(n_vertices);

	cell_neighbor_iterator.resize(n_locally_cells);
	cell_all_neighbor_iterator.resize(n_locally_cells);

	typename std::set<unsigned int>::iterator ver_iter = vertex_index.begin(), ver_end = vertex_index.end();

	unsigned int local_vertex = 0;
	for (; ver_iter != ver_end; ++ver_iter) {
		vertex_map[*ver_iter] = local_vertex;
		local_vertex++;
	}

	cell = dof_handler.begin_active();

    for (; cell!=endc; ++cell) {

		if(cell->is_locally_owned()) {

	        for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i){

				vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].insert(cell);
			
			}

		}

		if(cell->is_ghost() ) {

	        for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i){

				vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].insert(cell);

			}

		}
    }

	cell = dof_handler.begin_active();

	typedef typename std::set<DoFHandler<2>::active_cell_iterator>::iterator ver_cell_iter;

    for (; cell!=endc; ++cell) {

		if(cell->is_locally_owned()) {

	        cell->get_dof_indices(local_dof_indices);
			unsigned int local_i = global_to_local_index_map[local_dof_indices[0] ];		

	        for (unsigned int i=0; i<GeometryInfo<2>::vertices_per_cell; ++i){

				ver_cell_iter adjacent_cell = vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].begin();
				ver_cell_iter end_cell = vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].end();

				for(;adjacent_cell != end_cell; ++adjacent_cell)
					cell_neighbor_iterator[local_i].insert(*adjacent_cell);
			
			}

			cell_neighbor_iterator[local_i].erase(cell);
			cell_all_neighbor_iterator[local_i] = cell_neighbor_iterator[local_i] ;

			for (unsigned int f=0; f < 4; ++f) {

                if (!cell->face(f)->at_boundary()) {
					cell_neighbor_iterator[local_i].erase(cell->neighbor(f));
				}
			}
//			if(cell_neighbor_iterator[local_i].size() == 0 ) std::cout<<local_dof_indices[0]<<"\tcenter: "<<cell->center()<<std::endl;
		}
    }	

    pcout << "Done!" << std::endl;  
    pcout << "===========================" << std::endl;
}

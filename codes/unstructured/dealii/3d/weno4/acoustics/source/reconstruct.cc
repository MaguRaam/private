#include "../include/Weno432.h"


// Perform the actual reconstruction 

void Weno4_3D::reconstruct() {

    unsigned int n_faces = GeometryInfo<3>::faces_per_cell;

    double h; // Measure of cell size 
    
    // Variables for reconstruction of RHO_U

	double u0; 

	// Fourth order stencil
	Vector<double> d_u_4; 
	Vector<double> b_u_4; upwind_coeffs_RHO_V[c](i+1) = rho_v_coeff_4(i);
	Vector<double> u_coeff_4(19);


    // Iterate over all the cells 

    DoFHandler<3>::active_cell_iterator cell, neighbor;

	std::vector<bool> n_n_face(6, false);
	std::vector<unsigned int> n_n_face_index(6), neighbor_index(6, numbers::invalid_unsigned_int);  
	unsigned int index_least;
    unsigned int ROWS, ROWS_slip, index, g_i, fff, n_boundary_quad_points = 2, n_boundary_quad_points_edge = 1;
	bool negative = false;
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];
		is_1st_order[c] = false;

      	for (unsigned int f = 0; f < n_faces; f++){	
			if (cell_neighbor_neighbor_index[c][f].size() > 0) {
				n_n_face[f] = true;
				n_n_face_index[f] = cell_neighbor_neighbor_index[c][f][0];
			}
			else {
				n_n_face[f] = false;
				n_n_face_index[f] = numbers::invalid_unsigned_int;
			}

			if(!cell->face(f)->at_boundary()) {
				neighbor = cell->neighbor(f);
			    neighbor->get_dof_indices(local_neighbor_dof_indices);
				neighbor_index[f] = local_neighbor_dof_indices[0];
			}
			else
				neighbor_index[f] = numbers::invalid_unsigned_int;
		}

    	u0 =   U(g_i);

		h = Cell[c].h();
        
        if ( !(cell->at_boundary()) ) {
//			pcout<<"interior g_i: "<<g_i<<"\tcenter: "<<cell->center()<<std::endl;
            // =====================================================================
            // r = 4 stencil upwind_coeffs_RHO_V[c](i+1) = rho_v_coeff_4(i);
            	
				d_u_4(f)  = (U(neighbor_index[f]) - rho_u0);
			} 

			ROWS = cell_diagonal_neighbor_index[c].size(); 
			
			// neighbor of neighbors 

	      	for (unsigned int f = 0; f < n_faces; f++)
				if (n_n_face[f])
					ROWS++; 
				
			index = 0; 
            
            // Least Squares Part  
            
			b_u_4.reinit(ROWS);	
			
			// vertex neighbors
			
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];
				b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0;  
				index++; 
			}
			
			// neighbors of neighbors 

	      	for (unsigned int f = 0; f < n_faces; f++) {
				if (n_n_face[f]){
					local_neighbor_dof_indices[0] = n_n_face_index[f];
					b_u_4(index)   = U(local_neighbor_dof_indices[0] ) - u0;
					index++;
				}
            }

			CLS_R4[c].solve(b_u_4, d_u_4, u_coeff_4); 
		} // End of interior cells loop

		else {

				unsigned int n_boundary_faces = 0, n_no_slip_boundary = 0;
			    for (unsigned int f=0; f<n_faces; ++f) {
					if(cell->face(f)->at_boundary()) {
						if(cell->face(f)->boundary_id() == 2) {	
							n_no_slip_boundary += 1;
						}
						n_boundary_faces += 1;
					}
				}

				// =====================================================================
	            // r = 4 stencil 
	            // =====================================================================

				index_least = 0;

				ROWS = cell_diagonal_neighbor_index[c].size(); 
			
				// neighbor of neighbors 

		      	for (unsigned int f = 0; f < n_faces; f++)
					if (n_n_face[f])
						ROWS++; 

				ROWS_slip = ROWS + (n_boundary_faces - n_no_slip_boundary)*4;
				ROWS += n_boundary_faces*4;
    	        // Least Squares Part  
				b_u_4.reinit(ROWS_slip);	
			
				// vertex neighbors
			
				for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

					local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];
					b_u_4(index_least)   = U(local_neighbor_dof_indices[0] ) - u0;  
					index_least++; 
				}
			
				// neighbors of neighbors 

		      	for (unsigned int f = 0; f < n_faces; f++) {
					if (n_n_face[f]){
						local_neighbor_dof_indices[0] = n_n_face_index[f];
						b_u_4(index_least)   = U(local_neighbor_dof_indices[0] ) - u0;
						index_least++;
					}
    	        }

				index = 0; 
				ROWS_slip = 6 - n_boundary_faces + n_no_slip_boundary*4;
				ROWS = 6 - n_boundary_faces;
				d_u_4.reinit(ROWS_slip);
				unsigned int index_least_slip = index_least;

				index = 0; 
				ROWS_slip = 6 - n_boundary_faces + n_no_slip_boundary*4;
				d_u_4.reinit(ROWS_slip);
				
	           	for (unsigned int f = 0; f < n_faces; f++) {			

					if(!cell->face(f)->at_boundary()) {	
						neighbor = cell->neighbor(f);
					    neighbor->get_dof_indices(local_neighbor_dof_indices);
						neighbor_index[f] = local_neighbor_dof_indices[0];
						d_u_4(index)  = (U(neighbor_index[f]) - u0);
						index++;
					}
					else {
						if (cell->face(f)->boundary_id() == 2) {
			                for (unsigned int i = 0; i < 4; i++) {
								d_u_4(index) = -u0;
								index++;
							}
						}
						else {
			                for (unsigned int i = 0; i < 4; i++) {
								b_u_4(index_least_slip) = 0.0;
								index_least_slip++;
							}
						}
					}
				}

                CLS_R4_slip[c].solve(b_u_4, d_u_4, u_coeff_4);
                

		coeffs_U[c](0) = u0;
	
		upwind_coeffs_U[c](0) = u0;
		
		for (unsigned int i = 0; i < 19; ++i) {

			coeffs_U[c](i+1) = u_coeff_4(i);
			
			upwind_coeffs_U[c](i+1) = u_coeff_4(i);
		}

    } // End of cell loop 
} // End of function 

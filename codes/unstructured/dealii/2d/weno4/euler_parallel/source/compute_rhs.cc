#include "../include/Weno432.h"

void Weno4_2D::compute_rhs() {

    reconstruct();
    
    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;

	// Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell, neighbor ;

    Point<2> face_quadrature_point_1;
    Point<2> face_quadrature_point_2;

    double V_c, S_f, h; // Volume of the cell and surface area of the face
	
	double nx1, ny1;   // Face normal vectors
	double nx2, ny2; 
	
    Vector<double> UL1(4); Vector<double> UR1(4); // Solving the Riemann Problem
    Vector<double> UL2(4); Vector<double> UR2(4); // Solving the Riemann Problem

	Vector<double> WL1(4); Vector<double> WR1(4); // Solving the Riemann Problem
    Vector<double> WL2(4); Vector<double> WR2(4); // Solving the Riemann Problem
	
    Vector<double> F(4);
    
    bool boundary;

	unsigned int local_face_index;  

    std::vector< Vector<double> > Flux1(n_faces); 
    std::vector< Vector<double> > Flux2(n_faces);
    std::vector< bool > did_not_compute_flux_for_the_face(n_faces); 
    
    for (unsigned int f = 0; f < n_faces; f++) {
        Flux1[f].reinit(4);
        Flux2[f].reinit(4);
        did_not_compute_flux_for_the_face[f] = true; 
    }
    
    
	unsigned int neighbor_c, g_i ;

	for (unsigned int c = 0; c < n_locally_cells; ++c) {

		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];    

        V_c = Cell[c].measure();

        rhs1(c) = 0.0;
        rhs2(c) = 0.0;
        rhs3(c) = 0.0;
        rhs4(c) = 0.0;
        
        
        for (unsigned int f = 0; f < faces_per_cell; ++f) {
            
			local_face_index = face_index_map[ cell->face_index(f) ];
            S_f = Cell[c].S_f(f);
            
            if(did_not_compute_flux_for_the_face[local_face_index]) {
               
                // Get some geometry info

				h = Cell[c].h(); 
    
                nx1 = Cell[c].nx1(f); ny1 = Cell[c].ny1(f);
                nx2 = Cell[c].nx2(f); ny2 = Cell[c].ny2(f);
                
                face_quadrature_point_1 = Cell[c].face_quadrature_point1(f); 
                face_quadrature_point_2 = Cell[c].face_quadrature_point2(f); 
                
                // Left face
                UL1(0) = evaluate_weno_polynomial(coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point_1, h);
                UL1(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point_1, h);
                UL1(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point_1, h);
                UL1(3) = evaluate_weno_polynomial(coeffs_E[c], WENO_poly_consts[c], face_quadrature_point_1, h);
                
                UL2(0) = evaluate_weno_polynomial(coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point_2, h);
                UL2(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point_2, h);
                UL2(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point_2, h);
                UL2(3) = evaluate_weno_polynomial(coeffs_E[c], WENO_poly_consts[c], face_quadrature_point_2, h);
                
                
                WL1 = conserved_to_primitive(UL1); WL2 = conserved_to_primitive(UL2);
	            
				if (WL1(0) < 0.0 || WL1(3) < 0.0 || WL2(0) < 0.0 || WL2(3) < 0.0) {
    	            std::cout<<"in compute rhs "<<std::endl<<"global index: "<<g_i<<"\t local index: "<<c
                    <<"\t rank: "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<"\tface: "<<f<<"\tquad point 1: "<<face_quadrature_point_1<<std::endl
                    <<WL1<<std::endl<<WL2<<std::endl
					<<coeffs_E[c]<<std::endl<<WENO_poly_consts[c]<<std::endl;
                }
/*
				if (c == 0) {
    	            std::cout<<"in compute rhs "<<std::endl<<"global index: "<<c
                    <<"\tface: "<<f<<"\tquad point 1: "<<face_quadrature_point_1<<"\th: "<<h<<std::endl
                    <<WL1<<std::endl<<WL2<<std::endl
					<<coeffs_E[c]<<std::endl<<WENO_poly_consts[c]<<std::endl;
                }
*/                
                if (cell->face(f)->at_boundary()) {
                                          
					WR1 = WL1; WR2 = WL2; 

					boundary = true;                  
                    
                }
                
                else {
                    
                    // Get the right state values

	    	        neighbor = cell->neighbor(f);
		            neighbor->get_dof_indices(local_neighbor_dof_indices);
					neighbor_c = global_to_local_index_map[local_neighbor_dof_indices[0] ];
					
					h = Cell[neighbor_c].h();
                    
                    UR1(0) = evaluate_weno_polynomial(coeffs_RHO[neighbor_c ], WENO_poly_consts[neighbor_c ],  face_quadrature_point_1, h);
                    UR1(1) = evaluate_weno_polynomial(coeffs_RHO_U[neighbor_c ], WENO_poly_consts[neighbor_c ],  face_quadrature_point_1, h);
                    UR1(2) = evaluate_weno_polynomial(coeffs_RHO_V[neighbor_c ], WENO_poly_consts[neighbor_c ],  face_quadrature_point_1, h);
                    UR1(3) = evaluate_weno_polynomial(coeffs_E[neighbor_c ], WENO_poly_consts[neighbor_c ],  face_quadrature_point_1, h);
                    
                    UR2(0) = evaluate_weno_polynomial(coeffs_RHO[neighbor_c ], WENO_poly_consts[neighbor_c ],  face_quadrature_point_2, h);
                    UR2(1) = evaluate_weno_polynomial(coeffs_RHO_U[neighbor_c ], WENO_poly_consts[neighbor_c ],  face_quadrature_point_2, h);
                    UR2(2) = evaluate_weno_polynomial(coeffs_RHO_V[neighbor_c ], WENO_poly_consts[neighbor_c ],  face_quadrature_point_2, h);
                    UR2(3) = evaluate_weno_polynomial(coeffs_E[neighbor_c ], WENO_poly_consts[neighbor_c ],  face_quadrature_point_2, h);
                
                    WR1 = conserved_to_primitive(UR1); WR2 = conserved_to_primitive(UR2);

		            if (WR1(0) < 0.0 || WR1(3) < 0.0 || WR2(0) < 0.0 || WR2(3) < 0.0) {
						std::cout<<"in compute rhs "<<std::endl<<"global index: "<<g_i<<"\t local index: "<<c<<std::endl<<"neighbor global index: "<<local_neighbor_dof_indices[0] <<"\t local index: "<<neighbor_c
						<<"\t rank: "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<"\tface: "<<f<<"\tcenter: "<<cell->face(f)->center()<<std::endl
						<<WR1<<std::endl<<WR2<<std::endl<<coeffs_E[neighbor_c]<<std::endl<<WENO_poly_consts[neighbor_c]<<std::endl;
					}   
/*
					if(c == 1) {
						std::cout<<"in compute rhs "<<std::endl<<"neighbor global index: "<<neighbor_c
						<<"\tface: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\th: "<<h<<std::endl
						<<WR1<<std::endl<<WR2<<std::endl<<coeffs_E[neighbor_c]<<std::endl<<WENO_poly_consts[neighbor_c]<<std::endl;
					}                 
					if( c == 1 ) {
						std::cout<<"coeffs"<<std::endl;
						for(unsigned int i = 0; i < 10; ++i )
							std::cout<<coeffs_RHO[neighbor_c](i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 10; ++i )
							std::cout<<coeffs_RHO_U[neighbor_c](i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 10; ++i )
							std::cout<<coeffs_RHO_V[neighbor_c](i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 10; ++i )
							std::cout<<coeffs_E[neighbor_c](i)<<"\t";
						std::cout<<std::endl;
						std::cout<<"Poly"<<std::endl;
						for(unsigned int i = 0; i < 9; ++i )
							std::cout<<WENO_poly_consts[neighbor_c](i)<<"\t";
						std::cout<<std::endl;							
					}		
*/                    
                    boundary = false; 
                }
                    
                
                Flux1[local_face_index] = rotated_HLLC_riemann_solver(WL1, WR1, nx1, ny1, face_quadrature_point_1, boundary);
                Flux2[local_face_index] = rotated_HLLC_riemann_solver(WL2, WR2, nx2, ny2, face_quadrature_point_2, boundary);
                
                
                did_not_compute_flux_for_the_face[local_face_index] = false; 
            }
            
            else {
                
                Flux1[local_face_index] *= -1.0; 
                Flux2[local_face_index] *= -1.0;
            }
            
            
            F(0) = 0.5*(Flux1[local_face_index](0) + Flux2[local_face_index](0)); 
            F(1) = 0.5*(Flux1[local_face_index](1) + Flux2[local_face_index](1)); 
            F(2) = 0.5*(Flux1[local_face_index](2) + Flux2[local_face_index](2));
            F(3) = 0.5*(Flux1[local_face_index](3) + Flux2[local_face_index](3));

            // Add it to the rhs vectors

            rhs1(c) += (-1.0/V_c)*(F(0)*S_f);
            rhs2(c) += (-1.0/V_c)*(F(1)*S_f);
            rhs3(c) += (-1.0/V_c)*(F(2)*S_f);
            rhs4(c) += (-1.0/V_c)*(F(3)*S_f);
        }
	}
}


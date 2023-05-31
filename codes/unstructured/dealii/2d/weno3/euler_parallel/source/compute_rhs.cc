#include "../include/Weno32.h"


// Compute the rhs vectors

void Weno3_2D::compute_rhs() {

    reconstruct();
   
    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;

	// Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

    Point<2> face_quadrature_point_1;
    Point<2> face_quadrature_point_2;

    double V_c, j_w1, j_w2; // Volume of the cell and jacobian value for surface area of the face
	
	double nx1, ny1;   // Face normal vectors
	double nx2, ny2; 
	
    Vector<double> UL1(4); Vector<double> UR1(4); // Solving the Riemann Problem
    Vector<double> UL2(4); Vector<double> UR2(4); // Solving the Riemann Problem

	Vector<double> WL1(4); Vector<double> WR1(4); // Solving the Riemann Problem
    Vector<double> WL2(4); Vector<double> WR2(4); // Solving the Riemann Problem
    
    bool boundary;

    Vector<double>  Flux1(4); 
    Vector<double>  Flux2(4);
    Vector<double> F(4);

	unsigned int c, g_i ;
    
	for (; cell != endc; ++cell) {

   	  if (cell->is_locally_owned()){

	        cell->get_dof_indices(local_dof_indices);
			c = global_to_local_index_map[local_dof_indices[0] ];
        	g_i = local_dof_indices[0] ;
			V_c = Cell[c].measure();

			double flux_rho = 0.0;	double flux_rho_u = 0.0; 
			double flux_rho_v = 0.0;	double flux_e = 0.0; 

			Vector<double> local_RHO_coeffs(6);
			Vector<double> local_RHO_U_coeffs(6);
			Vector<double> local_RHO_V_coeffs(6);
			Vector<double> local_E_coeffs(6);     

			local_RHO_coeffs(0) = RHO(g_i);
			local_RHO_U_coeffs(0) = RHO_U(g_i);
			local_RHO_V_coeffs(0) = RHO_V(g_i);
			local_E_coeffs(0) = E(g_i);

			local_RHO_coeffs(1) = coeffs_x_RHO(g_i);
			local_RHO_U_coeffs(1) = coeffs_x_RHO_U(g_i);
			local_RHO_V_coeffs(1) = coeffs_x_RHO_V(g_i);
			local_E_coeffs(1) = coeffs_x_E(g_i);

			local_RHO_coeffs(2) = coeffs_y_RHO(g_i);
			local_RHO_U_coeffs(2) = coeffs_y_RHO_U(g_i);
			local_RHO_V_coeffs(2) = coeffs_y_RHO_V(g_i);
			local_E_coeffs(2) = coeffs_y_E(g_i); 

			local_RHO_coeffs(3) = coeffs_xx_RHO(g_i);
			local_RHO_U_coeffs(3) = coeffs_xx_RHO_U(g_i);
			local_RHO_V_coeffs(3) = coeffs_xx_RHO_V(g_i);
			local_E_coeffs(3) = coeffs_xx_E(g_i);

			local_RHO_coeffs(4) = coeffs_yy_RHO(g_i);
			local_RHO_U_coeffs(4) = coeffs_yy_RHO_U(g_i);
			local_RHO_V_coeffs(4) = coeffs_yy_RHO_V(g_i);
			local_E_coeffs(4) = coeffs_yy_E(g_i);

			local_RHO_coeffs(5) = coeffs_xy_RHO(g_i);
			local_RHO_U_coeffs(5) = coeffs_xy_RHO_U(g_i);
			local_RHO_V_coeffs(5) = coeffs_xy_RHO_V(g_i);
			local_E_coeffs(5) = coeffs_xy_E(g_i);  

			Vector<double> local_WENO_poly_consts(5);

			local_WENO_poly_consts(0) = WENO_poly_consts_x(g_i);
			local_WENO_poly_consts(1) = WENO_poly_consts_y(g_i);
			local_WENO_poly_consts(2) = WENO_poly_consts_xx(g_i);
			local_WENO_poly_consts(3) = WENO_poly_consts_yy(g_i);
			local_WENO_poly_consts(4) = WENO_poly_consts_xy(g_i);
        
			for (unsigned int f=0; f < faces_per_cell; ++f) {
            
				j_w1 = Cell[c].face_jxw_1(f);
				j_w2 = Cell[c].face_jxw_2(f);
                               
                nx1 = Cell[c].nx1(f); ny1 = Cell[c].ny1(f);
                nx2 = Cell[c].nx2(f); ny2 = Cell[c].ny2(f);
                
                face_quadrature_point_1 = Cell[c].face_quadrature_point1(f);
                face_quadrature_point_2 = Cell[c].face_quadrature_point2(f);
                
                // Left face
                UL1(0) = evaluate_weno_polynomial(local_RHO_coeffs, local_WENO_poly_consts, face_quadrature_point_1);
                UL1(1) = evaluate_weno_polynomial(local_RHO_U_coeffs, local_WENO_poly_consts, face_quadrature_point_1);
                UL1(2) = evaluate_weno_polynomial(local_RHO_V_coeffs, local_WENO_poly_consts, face_quadrature_point_1);
                UL1(3) = evaluate_weno_polynomial(local_E_coeffs, local_WENO_poly_consts, face_quadrature_point_1);
                
                UL2(0) = evaluate_weno_polynomial(local_RHO_coeffs, local_WENO_poly_consts, face_quadrature_point_2);
                UL2(1) = evaluate_weno_polynomial(local_RHO_U_coeffs, local_WENO_poly_consts, face_quadrature_point_2);
                UL2(2) = evaluate_weno_polynomial(local_RHO_V_coeffs, local_WENO_poly_consts, face_quadrature_point_2);
                UL2(3) = evaluate_weno_polynomial(local_E_coeffs, local_WENO_poly_consts, face_quadrature_point_2);
                
                WL1 = conserved_to_primitive(UL1); WL2 = conserved_to_primitive(UL2);

              
                if (cell->face(f)->at_boundary()) {
                    if (cell->face(f)->boundary_id() == 0 || cell->face(f)->boundary_id() == 1) {
                    
                        // Left and Right (Reflecting)
                        
                        WR1(0) = WL1(0);  
                        WR1(1) = WL1(1) - 2.0*WL1(1)*nx1*nx1 - 2.0*WL1(2)*nx1*ny1; 
                        WR1(2) = WL1(2) - 2.0*WL1(1)*nx1*ny1 - 2.0*WL1(2)*ny1*ny1;
                        WR1(3) = WL1(3);  
                    
						WR2(0) = WL2(0);  
                        WR2(1) = WL2(1) - 2.0*WL2(1)*nx2*nx2 - 2.0*WL2(2)*nx2*ny2; 
                        WR2(2) = WL2(2) - 2.0*WL2(1)*nx2*ny2 - 2.0*WL2(2)*ny2*ny2;
                        WR2(3) = WL2(3);  
						
                    }
                    
                    else if (cell->face(f)->boundary_id() == 2) {
                        
                        // Bottom boundary 
                        
                        WR1(0) =  2.0;
                        WR1(1) =  0.0;
                        WR1(2) =  0.0;
                        WR1(3) =  1.0;
						
						WR2 = WR1;  
                    }
                    
                    
                    else if (cell->face(f)->boundary_id() == 3) {
                        
                        // Top boundary 
                        
                        WR1(0) =  1.0;
                        WR1(1) =  0.0;
                        WR1(2) =  0.0;
                        WR1(3) =  2.5;
						
						WR2 = WR1; 
                    }
                    
                    boundary = true; 
    	        }
                else {

	    	        typename DoFHandler<2>::cell_iterator neighbor = cell->neighbor(f);
		            neighbor->get_dof_indices(local_neighbor_dof_indices);
	
					Vector<double> neighbor_RHO_coeffs(6);
					Vector<double> neighbor_RHO_U_coeffs(6);
					Vector<double> neighbor_RHO_V_coeffs(6);
					Vector<double> neighbor_E_coeffs(6);

					neighbor_RHO_coeffs(0) = RHO(local_neighbor_dof_indices[0]);
					neighbor_RHO_U_coeffs(0) = RHO_U(local_neighbor_dof_indices[0]);
					neighbor_RHO_V_coeffs(0) = RHO_V(local_neighbor_dof_indices[0]);
					neighbor_E_coeffs(0) = E(local_neighbor_dof_indices[0]);

					neighbor_RHO_coeffs(1) = coeffs_x_RHO(local_neighbor_dof_indices[0]);
					neighbor_RHO_U_coeffs(1) = coeffs_x_RHO_U(local_neighbor_dof_indices[0]);
					neighbor_RHO_V_coeffs(1) = coeffs_x_RHO_V(local_neighbor_dof_indices[0]);
					neighbor_E_coeffs(1) = coeffs_x_E(local_neighbor_dof_indices[0]);

					neighbor_RHO_coeffs(2) = coeffs_y_RHO(local_neighbor_dof_indices[0]);
					neighbor_RHO_U_coeffs(2) = coeffs_y_RHO_U(local_neighbor_dof_indices[0]);
					neighbor_RHO_V_coeffs(2) = coeffs_y_RHO_V(local_neighbor_dof_indices[0]);
					neighbor_E_coeffs(2) = coeffs_y_E(local_neighbor_dof_indices[0]); 

					neighbor_RHO_coeffs(3) = coeffs_xx_RHO(local_neighbor_dof_indices[0]);
					neighbor_RHO_U_coeffs(3) = coeffs_xx_RHO_U(local_neighbor_dof_indices[0]);
					neighbor_RHO_V_coeffs(3) = coeffs_xx_RHO_V(local_neighbor_dof_indices[0]);
					neighbor_E_coeffs(3) = coeffs_xx_E(local_neighbor_dof_indices[0]);

					neighbor_RHO_coeffs(4) = coeffs_yy_RHO(local_neighbor_dof_indices[0]);
					neighbor_RHO_U_coeffs(4) = coeffs_yy_RHO_U(local_neighbor_dof_indices[0]);
					neighbor_RHO_V_coeffs(4) = coeffs_yy_RHO_V(local_neighbor_dof_indices[0]);
					neighbor_E_coeffs(4) = coeffs_yy_E(local_neighbor_dof_indices[0]);

					neighbor_RHO_coeffs(5) = coeffs_xy_RHO(local_neighbor_dof_indices[0]);
					neighbor_RHO_U_coeffs(5) = coeffs_xy_RHO_U(local_neighbor_dof_indices[0]);
					neighbor_RHO_V_coeffs(5) = coeffs_xy_RHO_V(local_neighbor_dof_indices[0]);
					neighbor_E_coeffs(5) = coeffs_xy_E(local_neighbor_dof_indices[0]);  

					Vector<double> neighbor_WENO_poly_consts(5);

					neighbor_WENO_poly_consts(0) = WENO_poly_consts_x(local_neighbor_dof_indices[0]);
					neighbor_WENO_poly_consts(1) = WENO_poly_consts_y(local_neighbor_dof_indices[0]);
					neighbor_WENO_poly_consts(2) = WENO_poly_consts_xx(local_neighbor_dof_indices[0]);
					neighbor_WENO_poly_consts(3) = WENO_poly_consts_yy(local_neighbor_dof_indices[0]);
					neighbor_WENO_poly_consts(4) = WENO_poly_consts_xy(local_neighbor_dof_indices[0]);   	                
    	            // Get the right state values
                    
	                UR1(0) = evaluate_weno_polynomial(neighbor_RHO_coeffs, neighbor_WENO_poly_consts,  face_quadrature_point_1);
                    UR1(1) = evaluate_weno_polynomial(neighbor_RHO_U_coeffs, neighbor_WENO_poly_consts,  face_quadrature_point_1);
                    UR1(2) = evaluate_weno_polynomial(neighbor_RHO_V_coeffs, neighbor_WENO_poly_consts,  face_quadrature_point_1);
                    UR1(3) = evaluate_weno_polynomial(neighbor_E_coeffs, neighbor_WENO_poly_consts,  face_quadrature_point_1);
                    
                    UR2(0) = evaluate_weno_polynomial(neighbor_RHO_coeffs, neighbor_WENO_poly_consts,  face_quadrature_point_2);
                    UR2(1) = evaluate_weno_polynomial(neighbor_RHO_U_coeffs, neighbor_WENO_poly_consts,  face_quadrature_point_2);
                    UR2(2) = evaluate_weno_polynomial(neighbor_RHO_V_coeffs, neighbor_WENO_poly_consts,  face_quadrature_point_2);
                    UR2(3) = evaluate_weno_polynomial(neighbor_E_coeffs, neighbor_WENO_poly_consts,  face_quadrature_point_2);
                       
                    WR1 = conserved_to_primitive(UR1); WR2 = conserved_to_primitive(UR2);
                   
                    boundary = false; 
                }       
                
                Flux1 = local_Lax_Friedrichs_riemann_solver(WL1, WR1, nx1, ny1, face_quadrature_point_1, boundary);
                Flux2 = local_Lax_Friedrichs_riemann_solver(WL2, WR2, nx2, ny2, face_quadrature_point_2, boundary);
                                 
	            F(0) = j_w1 * Flux1(0) + j_w2 * Flux2(0); 
    	        F(1) = j_w1 * Flux1(1) + j_w2 * Flux2(1); 
    	        F(2) = j_w1 * Flux1(2) + j_w2 * Flux2(2);
    	        F(3) = j_w1 * Flux1(3) + j_w2 * Flux2(3);

    	        // Add it to the rhs vectors
	
    	        flux_rho += (-1.0/V_c)*F(0);
    	        flux_rho_u += (-1.0/V_c)*F(1);
    	        flux_rho_v += (-1.0/V_c)*F(2);
    	        flux_e += (-1.0/V_c)*F(3);

			}

			rhs1(c) = flux_rho;
			rhs2(c) = flux_rho_u;
			rhs3(c) = flux_rho_v;
			rhs4(c) = flux_e;
        }
    }
}

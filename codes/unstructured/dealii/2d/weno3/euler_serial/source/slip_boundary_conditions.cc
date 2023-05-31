#include "../include/Weno32.h"

// Impement the slip boundary conditions 

/*
void Weno4_2D::precompute_matrices_veclocity() {
    
    std::cout << "Computing reconstruction matrices for slip boundaries" << std::endl;
    
    unsigned int N_gp = 2;  // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);

    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    Point<2> q_point;
    Point<2> C; 
    Point<2> P; 
    Point<2> face_center; 
    
    double V_neighbor; 

    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    
    QGauss<2-1> face_quadrature_formula(2);
    FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);

    Tensor<1,2> face_normal_vector1; // Face normal vector
    Tensor<1,2> face_normal_vector2; // Face normal vector

    double nx1, ny1;   // Face normal vectors
    double nx2, ny2; 
    
    double x_g1, y_g1, x_g2, y_g2, x_c, y_c;


    FullMatrix<double> A_R3(6,10); // Least Squares Matrix for r=3 stencil 
    FullMatrix<double> C_R3(5,10); // Constraint Matrix for r=3 stencil
    
    FullMatrix<double> C_R4(7,18);
    FullMatrix<double> A_R4; 
    
    unsigned int ROWS, index; 
    

    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        if (cell->at_boundary()) {
            
            for (unsigned int f=0; f < GeometryInfo<2>::faces_per_cell; ++f) {
                
                if (cell->face(f)->at_boundary()) {
                    
                    if (cell->face(f)->boundary_id() == 2) { // boundary_id 2 corresponds to the reflecting wall 
                        
                        if (!(is_corner_cell[c])) {
                            
                            C = cell->center();
                            
                            // Get the normals 
                            
                            fv_face_values.reinit(cell, f); 

                            face_normal_vector1 = fv_face_values.normal_vector(0); 
                            nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                
                            face_normal_vector2 = fv_face_values.normal_vector(1);
                            nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
                
                            x_g1 = fv_face_values.quadrature_point(0)(0); 
                            y_g1 = fv_face_values.quadrature_point(0)(1);
                            x_g2 = fv_face_values.quadrature_point(1)(0);
                            y_g2 = fv_face_values.quadrature_point(1)(1);
                            face_center = cell->face(f)->center(); 
                            x_c = face_center(0); y_c = face_center(1);
                            
                            // =====================================================================
                            // r = 4 stencil  
                            // =====================================================================
                            
                            C_R4 = 0.0; 
                        
                            C_R4(0,0)  = nx1*(x_g1 - WENO_poly_consts[c](0)); 
                            C_R4(0,1)  = nx1*(y_g1 - WENO_poly_consts[c](1)); 
                            C_R4(0,2)  = nx1*(x_g1*x_g1 - WENO_poly_consts[c](2)); 
                            C_R4(0,3)  = nx1*(y_g1*y_g1 - WENO_poly_consts[c](3));
                            C_R4(0,4)  = nx1*(x_g1*y_g1 - WENO_poly_consts[c](4));
                            C_R4(0,5)  = nx1*(x_g1*x_g1*x_g1 - WENO_poly_consts[c](5));
                            C_R4(0,6)  = nx1*(y_g1*y_g1*y_g1 - WENO_poly_consts[c](6));
                            C_R4(0,7)  = nx1*(x_g1*x_g1*y_g1 - WENO_poly_consts[c](7));
                            C_R4(0,8)  = nx1*(x_g1*y_g1*y_g1 - WENO_poly_consts[c](8));
                            C_R4(0,9)  = ny1*(x_g1 - WENO_poly_consts[c](0)); 
                            C_R4(0,10) = ny1*(y_g1 - WENO_poly_consts[c](1)); 
                            C_R4(0,11) = ny1*(x_g1*x_g1 - WENO_poly_consts[c](2)); 
                            C_R4(0,12) = ny1*(y_g1*y_g1 - WENO_poly_consts[c](3));
                            C_R4(0,13) = ny1*(x_g1*y_g1 - WENO_poly_consts[c](4));
                            C_R4(0,14) = ny1*(x_g1*x_g1*x_g1 - WENO_poly_consts[c](5));
                            C_R4(0,15) = ny1*(y_g1*y_g1*y_g1 - WENO_poly_consts[c](6));
                            C_R4(0,16) = ny1*(x_g1*x_g1*y_g1 - WENO_poly_consts[c](7));
                            C_R4(0,17) = ny1*(x_g1*y_g1*y_g1 - WENO_poly_consts[c](8));
                            
                            C_R4(1,0)  = nx2*(x_g2 - WENO_poly_consts[c](0)); 
                            C_R4(1,1)  = nx2*(y_g2 - WENO_poly_consts[c](1)); 
                            C_R4(1,2)  = nx2*(x_g2*x_g2 - WENO_poly_consts[c](2)); 
                            C_R4(1,3)  = nx2*(y_g2*y_g2 - WENO_poly_consts[c](3));
                            C_R4(1,4)  = nx2*(x_g2*y_g2 - WENO_poly_consts[c](4));
                            C_R4(1,5)  = nx2*(x_g2*x_g2*x_g2 - WENO_poly_consts[c](5));
                            C_R4(1,6)  = nx2*(y_g2*y_g2*y_g2 - WENO_poly_consts[c](6));
                            C_R4(1,7)  = nx2*(x_g2*x_g2*y_g2 - WENO_poly_consts[c](7));
                            C_R4(1,8)  = nx2*(x_g2*y_g2*y_g2 - WENO_poly_consts[c](8));
                            C_R4(1,9)  = ny2*(x_g2 - WENO_poly_consts[c](0)); 
                            C_R4(1,10) = ny2*(y_g2 - WENO_poly_consts[c](1)); 
                            C_R4(1,11) = ny2*(x_g2*x_g2 - WENO_poly_consts[c](2)); 
                            C_R4(1,12) = ny2*(y_g2*y_g2 - WENO_poly_consts[c](3));
                            C_R4(1,13) = ny2*(x_g2*y_g2 - WENO_poly_consts[c](4));
                            C_R4(1,14) = ny2*(x_g2*x_g2*x_g2 - WENO_poly_consts[c](5));
                            C_R4(1,15) = ny2*(y_g2*y_g2*y_g2 - WENO_poly_consts[c](6));
                            C_R4(1,16) = ny2*(x_g2*x_g2*y_g2 - WENO_poly_consts[c](7));
                            C_R4(1,17) = ny2*(x_g2*y_g2*y_g2 - WENO_poly_consts[c](8));
                            
                            C_R4(2,0)  = nx2*(x_c - WENO_poly_consts[c](0)); 
                            C_R4(2,1)  = nx2*(y_c - WENO_poly_consts[c](1)); 
                            C_R4(2,2)  = nx2*(x_c*x_c - WENO_poly_consts[c](2)); 
                            C_R4(2,3)  = nx2*(y_c*y_c - WENO_poly_consts[c](3));
                            C_R4(2,4)  = nx2*(x_c*y_c - WENO_poly_consts[c](4));
                            C_R4(2,5)  = nx2*(x_c*x_c*x_c - WENO_poly_consts[c](5));
                            C_R4(2,6)  = nx2*(y_c*y_c*y_c - WENO_poly_consts[c](6));
                            C_R4(2,7)  = nx2*(x_c*x_c*y_c - WENO_poly_consts[c](7));
                            C_R4(2,8)  = nx2*(x_c*y_c*y_c - WENO_poly_consts[c](8));
                            C_R4(2,9)  = ny2*(x_c - WENO_poly_consts[c](0)); 
                            C_R4(2,10) = ny2*(y_c - WENO_poly_consts[c](1)); 
                            C_R4(2,11) = ny2*(x_c*x_c - WENO_poly_consts[c](2)); 
                            C_R4(2,12) = ny2*(y_c*y_c - WENO_poly_consts[c](3));
                            C_R4(2,13) = ny2*(x_c*y_c - WENO_poly_consts[c](4));
                            C_R4(2,14) = ny2*(x_c*x_c*x_c - WENO_poly_consts[c](5));
                            C_R4(2,15) = ny2*(y_c*y_c*y_c - WENO_poly_consts[c](6));
                            C_R4(2,16) = ny2*(x_c*x_c*y_c - WENO_poly_consts[c](7));
                            C_R4(2,17) = ny2*(x_c*y_c*y_c - WENO_poly_consts[c](8));
                            
                            // Row 4, 5 and 6 - Zero derivative for tangential velocity
                            
                            C_R4(3,0) = -ny1*nx1; 
                            C_R4(3,1) = -ny1*ny1;
                            C_R4(3,2) = -ny1*2.0*x_g1*nx1;
                            C_R4(3,3) = -ny1*2.0*y_g1*ny1; 
                            C_R4(3,4) = -ny1*(nx1*y_g1 + ny1*x_g1); 
                            C_R4(3,5) = -ny1*3.0*nx1*x_g1*x_g1; 
                            C_R4(3,6) = -ny1*3.0*ny1*y_g1*y_g1; 
                            C_R4(3,7) = -ny1*(2.0*nx1*x_g1*y_g1 + ny1*x_g1*x_g1); 
                            C_R4(3,8) = -ny1*(2.0*ny1*x_g1*y_g1 + nx1*y_g1*y_g1);
                            C_R4(3,9)  = nx1*nx1; 
                            C_R4(3,10) = nx1*ny1;
                            C_R4(3,11) = nx1*2.0*x_g1*nx1;
                            C_R4(3,12) = nx1*2.0*y_g1*ny1; 
                            C_R4(3,13) = nx1*(nx1*y_g1 + ny1*x_g1); 
                            C_R4(3,14) = nx1*3.0*nx1*x_g1*x_g1; 
                            C_R4(3,15) = nx1*3.0*ny1*y_g1*y_g1; 
                            C_R4(3,16) = nx1*(2.0*nx1*x_g1*y_g1 + ny1*x_g1*x_g1); 
                            C_R4(3,17) = nx1*(2.0*ny1*x_g1*y_g1 + nx1*y_g1*y_g1);
                            
                            C_R4(4,0) = -ny2*nx2; 
                            C_R4(4,1) = -ny2*ny2;
                            C_R4(4,2) = -ny2*2.0*x_g2*nx2;
                            C_R4(4,3) = -ny2*2.0*y_g2*ny2; 
                            C_R4(4,4) = -ny2*(nx2*y_g2 + ny2*x_g2); 
                            C_R4(4,5) = -ny2*3.0*nx2*x_g2*x_g2; 
                            C_R4(4,6) = -ny2*3.0*ny2*y_g2*y_g2; 
                            C_R4(4,7) = -ny2*(2.0*nx2*x_g2*y_g2 + ny2*x_g2*x_g2); 
                            C_R4(4,8) = -ny2*(2.0*ny2*x_g2*y_g2 + nx2*y_g2*y_g2);
                            C_R4(4,9)  = nx2*nx2; 
                            C_R4(4,10) = nx2*ny2;
                            C_R4(4,11) = nx2*2.0*x_g2*nx2;
                            C_R4(4,12) = nx2*2.0*y_g2*ny2; 
                            C_R4(4,13) = nx2*(nx2*y_g2 + ny2*x_g2); 
                            C_R4(4,14) = nx2*3.0*nx2*x_g2*x_g2; 
                            C_R4(4,15) = nx2*3.0*ny2*y_g2*y_g2; 
                            C_R4(4,16) = nx2*(2.0*nx2*x_g2*y_g2 + ny2*x_g2*x_g2); 
                            C_R4(4,17) = nx2*(2.0*ny2*x_g2*y_g2 + nx2*y_g2*y_g2);
                            
                            C_R4(5,0) = -ny1*nx1; 
                            C_R4(5,1) = -ny1*ny1;
                            C_R4(5,2) = -ny1*2.0*x_c*nx1;
                            C_R4(5,3) = -ny1*2.0*y_c*ny1; 
                            C_R4(5,4) = -ny1*(nx1*y_c + ny1*x_c); 
                            C_R4(5,5) = -ny1*3.0*nx1*x_c*x_c; 
                            C_R4(5,6) = -ny1*3.0*ny1*y_c*y_c; 
                            C_R4(5,7) = -ny1*(2.0*nx1*x_c*y_c + ny1*x_c*x_c); 
                            C_R4(5,8) = -ny1*(2.0*ny1*x_c*y_c + nx1*y_c*y_c);
                            C_R4(5,9)  = nx1*nx1; 
                            C_R4(5,10) = nx1*ny1;
                            C_R4(5,11) = nx1*2.0*x_c*nx1;
                            C_R4(5,12) = nx1*2.0*y_c*ny1; 
                            C_R4(5,13) = nx1*(nx1*y_c + ny1*x_c); 
                            C_R4(5,14) = nx1*3.0*nx1*x_c*x_c; 
                            C_R4(5,15) = nx1*3.0*ny1*y_c*y_c; 
                            C_R4(5,16) = nx1*(2.0*nx1*x_c*y_c + ny1*x_c*x_c); 
                            C_R4(5,17) = nx1*(2.0*ny1*x_c*y_c + nx1*y_c*y_c);
                            
                            // Row 7 - Zero third derivative for tangential velocity
                            
                            C_R4(6,0) = 0.0; 
                            C_R4(6,1) = 0.0; 
                            C_R4(6,2) = 0.0; 
                            C_R4(6,3) = 0.0;
                            C_R4(6,4) = 0.0; 
                            C_R4(6,5) = -ny1*nx1; 
                            C_R4(6,6) = -ny1*ny1; 
                            C_R4(6,7) = 0.0;
                            C_R4(6,8) = 0.0;
                            C_R4(6,9) = 0.0; 
                            C_R4(6,10) = 0.0; 
                            C_R4(6,11) = 0.0; 
                            C_R4(6,12) = 0.0;
                            C_R4(6,13) = 0.0; 
                            C_R4(6,14) = nx1*nx1; 
                            C_R4(6,15) = nx1*ny1; 
                            C_R4(6,16) = 0.0;
                            C_R4(6,17) = 0.0;
                            
                            ROWS = 0; index = 0; 

                            if (cell->neighbor_index(0) != -1) {
                                ROWS++; // P1 cell
                            } 

                            if (cell->neighbor_index(3) != -1) {
                                ROWS++; // P2 cell 
                            }

                            if (cell->neighbor_index(1) != -1) {
                                ROWS++; // P3 cell
                            }

                            if (cell->neighbor_index(2) != -1) {
                                ROWS++; // P4 cell 
                            }

                            if (cell->neighbor_index(0) != -1) {
                                if (cell->neighbor(0)->neighbor_index(0) != -1) {
                                    ROWS++;   // S1 cell 
                                }
                            }

                            if (cell->neighbor_index(0) != -1) {
                                if (cell->neighbor(0)->neighbor_index(3) != -1) {
                                    ROWS++; // S2 cell 
                                }
                            }

                            if (cell->neighbor_index(3) != -1) {
                                if (cell->neighbor(3)->neighbor_index(3) != -1) {
                                    ROWS++; // S3 cell
                                }
                            }

                            if (cell->neighbor_index(3) != -1) {
                                if (cell->neighbor(3)->neighbor_index(1) != -1) {
                                    ROWS++; // S4 cell 
                                }
                            }

                            if (cell->neighbor_index(1) != -1) {
                                if (cell->neighbor(1)->neighbor_index(1) != -1) {
                                    ROWS++; // S5 cell 
                                }
                            }

                            if (cell->neighbor_index(1) != -1) {
                                if (cell->neighbor(1)->neighbor_index(2) != -1) {
                                    ROWS++; // S6 cell 
                                }
                            }

                            if (cell->neighbor_index(2) != -1) {
                                if (cell->neighbor(2)->neighbor_index(2) != -1) {
                                    ROWS++; // S7 cell
                                }
                            }

                            if (cell->neighbor_index(2) != -1) {
                                if (cell->neighbor(2)->neighbor_index(0) != -1) {
                                    ROWS++; // S8 cell
                                }
                            }

                            A_R4.reinit(2*ROWS, 18); 
                            
                            A_R4 = 0.0;


                            if (cell->neighbor_index(0) != -1) {
                                
                                // (P1 cell)

                                V_neighbor = cell->neighbor(0)->measure(); 
                                fv_values.reinit (cell->neighbor(0));
                        
                                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                    q_point = fv_values.quadrature_point(i);
                                    A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                    A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                    A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                    A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                    A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                    A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                    A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                    A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                    A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                }

                                A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                A_R4(index+ROWS, 11) = A_R4(index, 2);
                                A_R4(index+ROWS, 12) = A_R4(index, 3);
                                A_R4(index+ROWS, 13) = A_R4(index, 4);
                                A_R4(index+ROWS, 14) = A_R4(index, 5);
                                A_R4(index+ROWS, 15) = A_R4(index, 6);
                                A_R4(index+ROWS, 16) = A_R4(index, 7);
                                A_R4(index+ROWS, 17) = A_R4(index, 8);
                                
                                index++;

                            } 

                            if (cell->neighbor_index(3) != -1) {
                                
                                // (P2 cell)

                                V_neighbor = cell->neighbor(3)->measure(); 
                                fv_values.reinit (cell->neighbor(3));
                        
                                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                    q_point = fv_values.quadrature_point(i);
                                    A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                    A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                    A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                    A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                    A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                    A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                    A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                    A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                    A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                }

                                A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                A_R4(index+ROWS, 11) = A_R4(index, 2);
                                A_R4(index+ROWS, 12) = A_R4(index, 3);
                                A_R4(index+ROWS, 13) = A_R4(index, 4);
                                A_R4(index+ROWS, 14) = A_R4(index, 5);
                                A_R4(index+ROWS, 15) = A_R4(index, 6);
                                A_R4(index+ROWS, 16) = A_R4(index, 7);
                                A_R4(index+ROWS, 17) = A_R4(index, 8);

                                index++;

                            }

                            if (cell->neighbor_index(1) != -1) {
                                
                                // (P3 cell)

                                V_neighbor = cell->neighbor(1)->measure(); 
                                fv_values.reinit (cell->neighbor(1));
                        
                                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                    q_point = fv_values.quadrature_point(i);
                                    A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                    A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                    A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                    A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                    A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                    A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                    A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                    A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                    A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                }

                                A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                A_R4(index+ROWS, 11) = A_R4(index, 2);
                                A_R4(index+ROWS, 12) = A_R4(index, 3);
                                A_R4(index+ROWS, 13) = A_R4(index, 4);
                                A_R4(index+ROWS, 14) = A_R4(index, 5);
                                A_R4(index+ROWS, 15) = A_R4(index, 6);
                                A_R4(index+ROWS, 16) = A_R4(index, 7);
                                A_R4(index+ROWS, 17) = A_R4(index, 8);

                                index++;

                            }

                            if (cell->neighbor_index(2) != -1) {
                                
                                // (P4 cell)

                                V_neighbor = cell->neighbor(2)->measure(); 
                                fv_values.reinit (cell->neighbor(2));
                        
                                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                    q_point = fv_values.quadrature_point(i);
                                    A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                    A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                    A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                    A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                    A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                    A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                    A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                    A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                    A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                }

                                A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                A_R4(index+ROWS, 11) = A_R4(index, 2);
                                A_R4(index+ROWS, 12) = A_R4(index, 3);
                                A_R4(index+ROWS, 13) = A_R4(index, 4);
                                A_R4(index+ROWS, 14) = A_R4(index, 5);
                                A_R4(index+ROWS, 15) = A_R4(index, 6);
                                A_R4(index+ROWS, 16) = A_R4(index, 7);
                                A_R4(index+ROWS, 17) = A_R4(index, 8);

                                index++; 
                            }


                            // S1 cell 

                            if (cell->neighbor_index(0) != -1) {
                                if (cell->neighbor(0)->neighbor_index(0) != -1) {

                                    V_neighbor = cell->neighbor(0)->neighbor(0)->measure(); 
                                    fv_values.reinit (cell->neighbor(0)->neighbor(0));
                                
                                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                        q_point = fv_values.quadrature_point(i);
                                        A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                        A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                        A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                        A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                        A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                        A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                        A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                        A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                        A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                    }

                                    A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                    A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                    A_R4(index+ROWS, 11) = A_R4(index, 2);
                                    A_R4(index+ROWS, 12) = A_R4(index, 3);
                                    A_R4(index+ROWS, 13) = A_R4(index, 4);
                                    A_R4(index+ROWS, 14) = A_R4(index, 5);
                                    A_R4(index+ROWS, 15) = A_R4(index, 6);
                                    A_R4(index+ROWS, 16) = A_R4(index, 7);
                                    A_R4(index+ROWS, 17) = A_R4(index, 8);

                                    index++; 
                                }
                            }

                            // S2 cell 

                            if (cell->neighbor_index(0) != -1) {
                                if (cell->neighbor(0)->neighbor_index(3) != -1) {

                                    V_neighbor = cell->neighbor(0)->neighbor(3)->measure(); 
                                    fv_values.reinit (cell->neighbor(0)->neighbor(3));
                                
                                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                        q_point = fv_values.quadrature_point(i);
                                        A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                        A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                        A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                        A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                        A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                        A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                        A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                        A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                        A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                    }

                                    A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                    A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                    A_R4(index+ROWS, 11) = A_R4(index, 2);
                                    A_R4(index+ROWS, 12) = A_R4(index, 3);
                                    A_R4(index+ROWS, 13) = A_R4(index, 4);
                                    A_R4(index+ROWS, 14) = A_R4(index, 5);
                                    A_R4(index+ROWS, 15) = A_R4(index, 6);
                                    A_R4(index+ROWS, 16) = A_R4(index, 7);
                                    A_R4(index+ROWS, 17) = A_R4(index, 8);

                                    index++; 
                                }
                            }

                            // S3 cell

                            if (cell->neighbor_index(3) != -1) {
                                if (cell->neighbor(3)->neighbor_index(3) != -1) {

                                    V_neighbor = cell->neighbor(3)->neighbor(3)->measure(); 
                                    fv_values.reinit (cell->neighbor(3)->neighbor(3));
                                
                                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                        q_point = fv_values.quadrature_point(i);
                                        A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                        A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                        A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                        A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                        A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                        A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                        A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                        A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                        A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                    }

                                    A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                    A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                    A_R4(index+ROWS, 11) = A_R4(index, 2);
                                    A_R4(index+ROWS, 12) = A_R4(index, 3);
                                    A_R4(index+ROWS, 13) = A_R4(index, 4);
                                    A_R4(index+ROWS, 14) = A_R4(index, 5);
                                    A_R4(index+ROWS, 15) = A_R4(index, 6);
                                    A_R4(index+ROWS, 16) = A_R4(index, 7);
                                    A_R4(index+ROWS, 17) = A_R4(index, 8);

                                    index++; 
                                }
                            }

                            // S4 cell 

                            if (cell->neighbor_index(3) != -1) {
                                if (cell->neighbor(3)->neighbor_index(1) != -1) {

                                    V_neighbor = cell->neighbor(3)->neighbor(1)->measure(); 
                                    fv_values.reinit (cell->neighbor(3)->neighbor(1));
                                
                                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                        q_point = fv_values.quadrature_point(i);
                                        A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                        A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                        A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                        A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                        A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                        A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                        A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                        A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                        A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                    }

                                    A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                    A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                    A_R4(index+ROWS, 11) = A_R4(index, 2);
                                    A_R4(index+ROWS, 12) = A_R4(index, 3);
                                    A_R4(index+ROWS, 13) = A_R4(index, 4);
                                    A_R4(index+ROWS, 14) = A_R4(index, 5);
                                    A_R4(index+ROWS, 15) = A_R4(index, 6);
                                    A_R4(index+ROWS, 16) = A_R4(index, 7);
                                    A_R4(index+ROWS, 17) = A_R4(index, 8);

                                    index++; 
                                }
                            }

                            // S5 cell

                            if (cell->neighbor_index(1) != -1) {
                                if (cell->neighbor(1)->neighbor_index(1) != -1) {

                                    V_neighbor = cell->neighbor(1)->neighbor(1)->measure(); 
                                    fv_values.reinit (cell->neighbor(1)->neighbor(1));
                                
                                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                        q_point = fv_values.quadrature_point(i);
                                        A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                        A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                        A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                        A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                        A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                        A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                        A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                        A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                        A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                    }

                                    A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                    A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                    A_R4(index+ROWS, 11) = A_R4(index, 2);
                                    A_R4(index+ROWS, 12) = A_R4(index, 3);
                                    A_R4(index+ROWS, 13) = A_R4(index, 4);
                                    A_R4(index+ROWS, 14) = A_R4(index, 5);
                                    A_R4(index+ROWS, 15) = A_R4(index, 6);
                                    A_R4(index+ROWS, 16) = A_R4(index, 7);
                                    A_R4(index+ROWS, 17) = A_R4(index, 8);

                                    index++; 
                                }
                            }

                            // S6 cell

                            if (cell->neighbor_index(1) != -1) {
                                if (cell->neighbor(1)->neighbor_index(2) != -1) {

                                    V_neighbor = cell->neighbor(1)->neighbor(2)->measure(); 
                                    fv_values.reinit (cell->neighbor(1)->neighbor(2));
                                
                                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                        q_point = fv_values.quadrature_point(i);
                                        A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                        A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                        A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                        A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                        A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                        A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                        A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                        A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                        A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                    }

                                    A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                    A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                    A_R4(index+ROWS, 11) = A_R4(index, 2);
                                    A_R4(index+ROWS, 12) = A_R4(index, 3);
                                    A_R4(index+ROWS, 13) = A_R4(index, 4);
                                    A_R4(index+ROWS, 14) = A_R4(index, 5);
                                    A_R4(index+ROWS, 15) = A_R4(index, 6);
                                    A_R4(index+ROWS, 16) = A_R4(index, 7);
                                    A_R4(index+ROWS, 17) = A_R4(index, 8);

                                    index++; 
                                }
                            }

                            // S7 cell

                            if (cell->neighbor_index(2) != -1) {
                                if (cell->neighbor(2)->neighbor_index(2) != -1) {

                                    V_neighbor = cell->neighbor(2)->neighbor(2)->measure(); 
                                    fv_values.reinit (cell->neighbor(2)->neighbor(2));
                                
                                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                        q_point = fv_values.quadrature_point(i);
                                        A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                        A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                        A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                        A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                        A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                        A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                        A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                        A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                        A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                    }

                                    A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                    A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                    A_R4(index+ROWS, 11) = A_R4(index, 2);
                                    A_R4(index+ROWS, 12) = A_R4(index, 3);
                                    A_R4(index+ROWS, 13) = A_R4(index, 4);
                                    A_R4(index+ROWS, 14) = A_R4(index, 5);
                                    A_R4(index+ROWS, 15) = A_R4(index, 6);
                                    A_R4(index+ROWS, 16) = A_R4(index, 7);
                                    A_R4(index+ROWS, 17) = A_R4(index, 8);

                                    index++; 
                                }
                            }
                            
                            // S8 cell

                            if (cell->neighbor_index(2) != -1) {
                                if (cell->neighbor(2)->neighbor_index(0) != -1) {

                                    V_neighbor = cell->neighbor(2)->neighbor(0)->measure(); 
                                    fv_values.reinit (cell->neighbor(2)->neighbor(0));
                                
                                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                        q_point = fv_values.quadrature_point(i);
                                        A_R4(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                        A_R4(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                        A_R4(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                        A_R4(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                        A_R4(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                        A_R4(index,5) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(0) - WENO_poly_consts[c](5));
                                        A_R4(index,6) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1)*q_point(1) - WENO_poly_consts[c](6));
                                        A_R4(index,7) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0)*q_point(1) - WENO_poly_consts[c](7));
                                        A_R4(index,8) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1)*q_point(1) - WENO_poly_consts[c](8));
                                    }

                                    A_R4(index+ROWS, 9)  = A_R4(index, 0); 
                                    A_R4(index+ROWS, 10) = A_R4(index, 1); 
                                    A_R4(index+ROWS, 11) = A_R4(index, 2);
                                    A_R4(index+ROWS, 12) = A_R4(index, 3);
                                    A_R4(index+ROWS, 13) = A_R4(index, 4);
                                    A_R4(index+ROWS, 14) = A_R4(index, 5);
                                    A_R4(index+ROWS, 15) = A_R4(index, 6);
                                    A_R4(index+ROWS, 16) = A_R4(index, 7);
                                    A_R4(index+ROWS, 17) = A_R4(index, 8);

                                    index++; 
                                }
                            }
                            
                            CLS_R4_slip[wall_boundary_global_index_map[c]].initialize(A_R4, C_R4);
                            
                            // =====================================================================
                            // r = 3 stencil  
                            // =====================================================================
                            
                            // Form the least squares matrix 
                            
                            A_R3 = 0.0; index = 0; 
                            
                            if (cell->neighbor_index(0) != -1 ) {
                            
                                // (P1 cell)
                                
                                V_neighbor = cell->neighbor(0)->measure(); 
                                fv_values.reinit (cell->neighbor(0));
                                 
                                
                                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                    q_point = fv_values.quadrature_point(i);
                                    A_R3(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                    A_R3(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                    A_R3(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                    A_R3(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                    A_R3(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                }
                                
                                A_R3(index+3,5) = A_R3(index,0); A_R3(index+3, 6) = A_R3(index,1); A_R3(index+3,7) = A_R3(index,2); 
                                A_R3(index+3,8) = A_R3(index,3); A_R3(index+3, 9) = A_R3(index,4);
                                
                                index++; 
                            }
                            
                            if (cell->neighbor_index(3) != -1 ) {
                            
                                // (P2 cell)
                            
                                V_neighbor = cell->neighbor(3)->measure(); 
                                fv_values.reinit (cell->neighbor(3));
                        
                                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                    q_point = fv_values.quadrature_point(i);
                                    A_R3(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                    A_R3(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                    A_R3(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                    A_R3(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                    A_R3(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                }
                                
                                A_R3(index+3,5) = A_R3(index,0); A_R3(index+3, 6) = A_R3(index,1); A_R3(index+3,7) = A_R3(index,2); 
                                A_R3(index+3,8) = A_R3(index,3); A_R3(index+3, 9) = A_R3(index,4);
                                
                                index++; 
                            }
                            
                            if (cell->neighbor_index(1) != -1 ) {
                            
                                // (P3 cell)
                            
                                V_neighbor = cell->neighbor(1)->measure(); 
                                fv_values.reinit (cell->neighbor(1));
                        
                                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                    q_point = fv_values.quadrature_point(i);
                                    A_R3(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                    A_R3(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                    A_R3(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                    A_R3(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                    A_R3(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                }
                                
                                A_R3(index+3,5) = A_R3(index,0); A_R3(index+3, 6) = A_R3(index,1); A_R3(index+3,7) = A_R3(index,2); 
                                A_R3(index+3,8) = A_R3(index,3); A_R3(index+3, 9) = A_R3(index,4);
                                
                                index++; 
                            }
                            
                            if (cell->neighbor_index(2) != -1 ) {
                            
                                // (P4 cell)
                            
                                V_neighbor = cell->neighbor(2)->measure(); 
                                fv_values.reinit (cell->neighbor(2));
                        
                                for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                                    q_point = fv_values.quadrature_point(i);
                                    A_R3(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                                    A_R3(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                                    A_R3(index,2) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(0) - WENO_poly_consts[c](2));
                                    A_R3(index,3) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1)*q_point(1) - WENO_poly_consts[c](3));
                                    A_R3(index,4) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0)*q_point(1) - WENO_poly_consts[c](4));
                                }
                                
                                A_R3(index+3,5) = A_R3(index,0); A_R3(index+3, 6) = A_R3(index,1); A_R3(index+3,7) = A_R3(index,2); 
                                A_R3(index+3,8) = A_R3(index,3); A_R3(index+3, 9) = A_R3(index,4);
                                
                            }
                            
                            // Form the constraint matrix 
                            
                            // Row 1, 2 and 3 - Zero normal velocity condition 
                            
                            C_R3(0,0) = nx1*(x_g1 - WENO_poly_consts[c](0)); C_R3(0,1) = nx1*(y_g1 - WENO_poly_consts[c](1)); 
                            C_R3(0,2) = nx1*(x_g1*x_g1 - WENO_poly_consts[c](2)); C_R3(0,3) = nx1*(y_g1*y_g1 - WENO_poly_consts[c](3));
                            C_R3(0,4) = nx1*(x_g1*y_g1 - WENO_poly_consts[c](4));
                            
                            C_R3(0,5) = ny1*(x_g1 - WENO_poly_consts[c](0)); C_R3(0,6) = ny1*(y_g1 - WENO_poly_consts[c](1)); 
                            C_R3(0,7) = ny1*(x_g1*x_g1 - WENO_poly_consts[c](2)); C_R3(0,8) = ny1*(y_g1*y_g1 - WENO_poly_consts[c](3));
                            C_R3(0,9) = ny1*(x_g1*y_g1 - WENO_poly_consts[c](4));
                            
                            C_R3(1,0) = nx2*(x_g2 - WENO_poly_consts[c](0)); C_R3(1,1) = nx2*(y_g2 - WENO_poly_consts[c](1)); 
                            C_R3(1,2) = nx2*(x_g2*x_g2 - WENO_poly_consts[c](2)); C_R3(1,3) = nx2*(y_g2*y_g2 - WENO_poly_consts[c](3));
                            C_R3(1,4) = nx2*(x_g2*y_g2 - WENO_poly_consts[c](4));
                            
                            C_R3(1,5) = ny2*(x_g2 - WENO_poly_consts[c](0)); C_R3(1,6) = ny2*(y_g2 - WENO_poly_consts[c](1)); 
                            C_R3(1,7) = ny2*(x_g2*x_g2 - WENO_poly_consts[c](2)); C_R3(1,8) = ny2*(y_g2*y_g2 - WENO_poly_consts[c](3));
                            C_R3(1,9) = ny2*(x_g2*y_g2 - WENO_poly_consts[c](4));
                            
                            C_R3(2,0) = nx1*(x_c - WENO_poly_consts[c](0)); C_R3(2,1) = nx2*(y_c - WENO_poly_consts[c](1)); 
                            C_R3(2,2) = nx1*(x_c*x_c - WENO_poly_consts[c](2)); C_R3(2,3) = nx2*(y_c*y_c - WENO_poly_consts[c](3));
                            C_R3(2,4) = nx1*(x_c*y_c - WENO_poly_consts[c](4));
                            
                            C_R3(2,5) = ny1*(x_c - WENO_poly_consts[c](0)); C_R3(2,6) = ny1*(y_c - WENO_poly_consts[c](1)); 
                            C_R3(2,7) = ny1*(x_c*x_c - WENO_poly_consts[c](2)); C_R3(2,8) = ny1*(y_c*y_c - WENO_poly_consts[c](3));
                            C_R3(2,9) = ny1*(x_c*y_c - WENO_poly_consts[c](4));
                            
                            // Row 4 and 5 - Zero derivative for tangential velocity 
                            
                            C_R3(3,0) = -nx1*ny1; C_R3(3,1) = -ny1*ny1; C_R3(3,2) = -2.0*nx1*ny1*x_g1; C_R3(3,3) = -2.0*ny1*ny1*y_g1; 
                            C_R3(3,4) = -nx1*ny1*y_g1 - ny1*ny1*x_g1; 
                            
                            C_R3(3,5) =  nx1*nx1; C_R3(3,6) =  ny1*nx1; C_R3(3,7) =  2.0*nx1*nx1*x_g1; C_R3(3,8) =  2.0*ny1*nx1*y_g1; 
                            C_R3(3,9) =  nx1*nx1*y_g1 + ny1*nx1*x_g1; 
                            
                            C_R3(4,0) = -nx2*ny2; C_R3(4,1) = -ny2*ny2; C_R3(4,2) = -2.0*nx2*ny2*x_g2; C_R3(4,3) = -2.0*ny2*ny2*y_g2; 
                            C_R3(4,4) = -nx2*ny2*y_g2 - ny2*ny2*x_g2; 
                            
                            C_R3(4,5) =  nx2*nx2; C_R3(4,6) =  ny2*nx2; C_R3(4,7) =  2.0*nx2*nx2*x_g2; C_R3(4,8) =  2.0*ny2*nx2*y_g2; 
                            C_R3(4,9) =  nx2*nx2*y_g2 + ny2*nx2*x_g2;
                            
                            CLS_R3_slip[wall_boundary_global_index_map[c]].initialize(A_R3, C_R3);
                            

                        } // End of non-corner boundary cell loop 
                    } // End of face_id 2 loop 
                } // End of cell face loop  
            } // End of boundary cell loop 
        } // End of boundary cell loop 
    } // End of cell loop 
    
    std::cout << "Done!" << std::endl;
    std::cout << "============================" << std::endl;
}

*/ 

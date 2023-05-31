#include "../include/Weno32.h"

// Precompute reconstruction matrices 

void Weno3_2D::precompute_matrices() {

    std::cout << "Computing the stencil matrices" << std::endl;

    unsigned int N_gp = 2;  // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);

    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    Point<2> q_point;
    Point<2> C; 
    Point<2> P; 
    
    double V_neighbor; 

    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    FullMatrix<double> A_R3; // Least Squares Matrix for r=3 stencil 
    FullMatrix<double> C_R3; // Constraint Matrix for r=3 stencil 
    
    FullMatrix<double> A_R21(2,2); // Matrix for r=2 stencil 1
    FullMatrix<double> A_R22(2,2); // Matrix for r=2 stencil 2
    FullMatrix<double> A_R23(2,2); // Matrix for r=2 stencil 3
    FullMatrix<double> A_R24(2,2); // Matrix for r=2 stencil 4
    
    double x0, y0;
    
    unsigned int ROWS, index; 
    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        x0 = WENO_poly_consts[c](0); 
        y0 = WENO_poly_consts[c](1);
        
        if ( !(cell->at_boundary()) ) {
            
            // =====================================================================
            // r = 3 stencil (Centered Stencil)
            // =====================================================================
            
            A_R3.reinit(4, 5); A_R3 = 0.0; 
            C_R3.reinit(4, 5); C_R3 = 0.0;
            
            // Fill A_R3 (Least squares r = 3 matrix)

            // Row 1 (S2 Cell)
            
            V_neighbor = cell->neighbor(0)->neighbor(3)->measure(); 
            fv_values.reinit (cell->neighbor(0)->neighbor(3));
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = fv_values.quadrature_point(i);
                A_R3(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                A_R3(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                A_R3(0,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                A_R3(0,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0)- - WENO_poly_consts[c](3));
                A_R3(0,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
            } 

            // Row 2 (S4 cell)

            V_neighbor = cell->neighbor(3)->neighbor(1)->measure(); 
            fv_values.reinit (cell->neighbor(3)->neighbor(1));
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = fv_values.quadrature_point(i);
                A_R3(1,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                A_R3(1,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                A_R3(1,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                A_R3(1,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0)- - WENO_poly_consts[c](3));
                A_R3(1,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
            }

            // Row 3 (S6 cell) 
            
            V_neighbor = cell->neighbor(1)->neighbor(2)->measure(); 
            fv_values.reinit (cell->neighbor(1)->neighbor(2));
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = fv_values.quadrature_point(i);
                A_R3(2,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                A_R3(2,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                A_R3(2,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                A_R3(2,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0)- - WENO_poly_consts[c](3));
                A_R3(2,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
            }

            // Row 4 (S8 cell) 
            
            V_neighbor = cell->neighbor(2)->neighbor(0)->measure(); 
            fv_values.reinit (cell->neighbor(2)->neighbor(0));
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = fv_values.quadrature_point(i);
                A_R3(3,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                A_R3(3,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                A_R3(3,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                A_R3(3,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0)- - WENO_poly_consts[c](3));
                A_R3(3,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
            }

            // Fill C_R3 (Constraint r = 3 matrix)

            // Row 1 (P1 cell)

            V_neighbor = cell->neighbor(0)->measure(); 
            fv_values.reinit (cell->neighbor(0));
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = fv_values.quadrature_point(i);
                C_R3(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                C_R3(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                C_R3(0,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                C_R3(0,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0)- - WENO_poly_consts[c](3));
                C_R3(0,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
            }

            // Row 2 (P2 cell)

            V_neighbor = cell->neighbor(3)->measure(); 
            fv_values.reinit (cell->neighbor(3));
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = fv_values.quadrature_point(i);
                C_R3(1,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                C_R3(1,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                C_R3(1,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                C_R3(1,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0)- - WENO_poly_consts[c](3));
                C_R3(1,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
            }

            // Row 3 (P3 cell) 
            
            V_neighbor = cell->neighbor(1)->measure(); 
            fv_values.reinit (cell->neighbor(1));
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = fv_values.quadrature_point(i);
                C_R3(2,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                C_R3(2,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                C_R3(2,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                C_R3(2,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0)- - WENO_poly_consts[c](3));
                C_R3(2,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
            }

            // Row 4 (P4 cell) 
            
            V_neighbor = cell->neighbor(2)->measure(); 
            fv_values.reinit (cell->neighbor(2));
            
            for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                q_point = fv_values.quadrature_point(i);
                C_R3(3,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                C_R3(3,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                C_R3(3,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                C_R3(3,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0)- - WENO_poly_consts[c](3));
                C_R3(3,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
            }

            CLS_R3[c].initialize(A_R3, C_R3);
            
            // =====================================================================
            // r = 2 stencil 1 (P1 and P2)
            // =====================================================================
            
            A_R21(0,0) = C_R3(0,0); A_R21(0,1) = C_R3(0,1); // Row 1 (P1 cell) 
            A_R21(1,0) = C_R3(1,0); A_R21(1,1) = C_R3(1,1); // Row 2 (P2 cell)

            LU_R21[c].initialize(A_R21);

            // =====================================================================
            // r = 2 stencil 2 (P2 and P3)
            // =====================================================================

            A_R22(0,0) = C_R3(1,0); A_R22(0,1) = C_R3(1,1); // Row 1 (P2 cell)
            A_R22(1,0) = C_R3(2,0); A_R22(1,1) = C_R3(2,1); // Row 2 (P3 cell)

            LU_R22[c].initialize(A_R22);

            // =====================================================================
            // r = 2 stencil 3 (P3 and P4)
            // =====================================================================

            A_R23(0,0) = C_R3(2,0); A_R23(0,1) = C_R3(2,1); // Row 1 (P3 cell)
            A_R23(1,0) = C_R3(3,0); A_R23(1,1) = C_R3(3,1); // Row 2 (P4 cell)

            LU_R23[c].initialize(A_R23);

            // =====================================================================
            // r = 2 stencil 4 (P4 and P1)
            // =====================================================================

            A_R24(0,0) = C_R3(3,0); A_R24(0,1) = C_R3(3,1); // Row 1 (P4 cell)
            A_R24(1,0) = C_R3(0,0); A_R24(1,1) = C_R3(0,1); // Row 2 (P1 cell)

            LU_R24[c].initialize(A_R24);

        } // End of interior cell loop 
        
        else {
            
            int p1, p2, p3, p4; 

            p1 = cell->neighbor_index(0); p2 = cell->neighbor_index(3);
            p3 = cell->neighbor_index(1); p4 = cell->neighbor_index(2); 

            // Mark cells at the corner  

            if ((p1==-1 && p2==-1) || (p2==-1 && p3 ==-1) || (p3==-1 && p4 ==-1) || (p4==-1 && p1 ==-1)) {
               is_corner_cell[c] = true;  
            }

            else {
                is_corner_cell[c] = false;
            }

            QGauss<2-1> face_quadrature_formula(2);
            FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);

            Tensor<1,2> face_normal_vector1; // Face normal vector
	        Tensor<1,2> face_normal_vector2; // Face normal vector

            double nx1, ny1;   // Face normal vectors
	        double nx2, ny2; 

            double x_g1, y_g1, x_g2, y_g2;  
            
            if (!(is_corner_cell[c])) {
                
                // =====================================================================
                // r = 3 centered stencil  
                // =====================================================================
                
                index = 0;  

                // Constraint matrix  
                C_R3.reinit(5, 5);

                C_R3 = 0.0; 
        
                // Fill C 

                if (cell->neighbor_index(0) != -1) {
                    
                    // P1 cell 

                    V_neighbor = cell->neighbor(0)->measure(); 
                    fv_values.reinit (cell->neighbor(0));
                
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        C_R3(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        C_R3(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                        C_R3(index,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                        C_R3(index,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0) - WENO_poly_consts[c](3));
                        C_R3(index,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
                    }
                    
                    index++; 
                }
                
                if (cell->neighbor_index(3) != -1) {

                    // Row 2 (P2 cell)

                    V_neighbor = cell->neighbor(3)->measure(); 
                    fv_values.reinit (cell->neighbor(3));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        C_R3(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        C_R3(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                        C_R3(index,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                        C_R3(index,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0) - WENO_poly_consts[c](3));
                        C_R3(index,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
                    }
                    
                    index++; 
                }
                
                if (cell->neighbor_index(1) != -1) {

                    // Row 3 (P3 cell) 
                    
                    V_neighbor = cell->neighbor(1)->measure(); 
                    fv_values.reinit (cell->neighbor(1));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        C_R3(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        C_R3(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                        C_R3(index,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                        C_R3(index,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0) - WENO_poly_consts[c](3));
                        C_R3(index,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
                    }
                    
                    index++; 
                
                }
                
                if (cell->neighbor_index(2) != -1) {

                    // Row 4 (P4 cell)

                    V_neighbor = cell->neighbor(2)->measure(); 
                    fv_values.reinit (cell->neighbor(2));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        C_R3(index,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        C_R3(index,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                        C_R3(index,2) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(0)-x0) - WENO_poly_consts[c](2));
                        C_R3(index,3) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(1)-y0)*(q_point(1)-y0) - WENO_poly_consts[c](3));
                        C_R3(index,4) += (1./V_neighbor)*fv_values.JxW (i)*((q_point(0)-x0)*(q_point(1)-y0) - WENO_poly_consts[c](4));
                    }
                    
                    index++; 
                }
                
                // initialize fv_face_values based on the boundary face 
                
                if (p1==-1) {
                    fv_face_values.reinit(cell, 0);
                }
                
                else if (p2 == -1) {
                    fv_face_values.reinit(cell, 3);
                }
                
                else if (p3 == -1) {
                    fv_face_values.reinit(cell, 1);
                }
                
                else {
                    fv_face_values.reinit(cell, 2);
                }
                
                face_normal_vector1 = fv_face_values.normal_vector(0); 
                nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                    
                face_normal_vector2 = fv_face_values.normal_vector(1);
                nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
                    
                x_g1 = fv_face_values.quadrature_point(0)(0); 
                y_g1 = fv_face_values.quadrature_point(0)(1);
                x_g2 = fv_face_values.quadrature_point(1)(0);
                y_g2 = fv_face_values.quadrature_point(1)(1);
                
                // Boundary Condition on Gauss Point 1 

                C_R3(index,0) = nx1; 
                C_R3(index,1) = ny1; 
                C_R3(index,2) = 2.0*(x_g1-x0)*nx1; 
                C_R3(index,3) = 2.0*(y_g1-y0)*ny1;
                C_R3(index,4) = (nx1*(y_g1-y0) + ny1*(x_g1-x0));

                // Boundary Condition on Gauss Point 2 

                C_R3(index+1,0) = nx2; 
                C_R3(index+1,1) = ny2; 
                C_R3(index+1,2) = 2.0*(x_g2-x0)*nx2; 
                C_R3(index+1,3) = 2.0*(y_g2-y0)*ny2;
                C_R3(index+1,4) = (nx2*(y_g2-y0) + ny2*(x_g2-x0)); 
                
                CLS_R3[c].initialize(C_R3); 

                // =====================================================================
                // r = 2 stencil 1 
                // =====================================================================
                
                if (p1 == -1) { // WEST 
                    // P2 cell 

                    A_R21 = 0.0; 
                    
                    V_neighbor = cell->neighbor(3)->measure(); 
                    fv_values.reinit (cell->neighbor(3));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R21(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R21(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    A_R21(1,0) = nx1;       A_R21(1,1) = ny1;       // Transmissive boundary
                }
                
                else if (p2 == -1) { // NORTH 
                    // P1 cell 

                    A_R21 = 0.0; 
                    
                    V_neighbor = cell->neighbor(0)->measure(); 
                    fv_values.reinit (cell->neighbor(0));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R21(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R21(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    A_R21(1,0) = nx1;       A_R21(1,1) = ny1;       // Transmissive boundary
                }
                
                else { // EAST and SOUTH 
                    // P1 cell 

                    A_R21 = 0.0; 
                    
                    V_neighbor = cell->neighbor(0)->measure(); 
                    fv_values.reinit (cell->neighbor(0));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R21(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R21(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    // P2 cell 
                    
                    V_neighbor = cell->neighbor(3)->measure(); 
                    fv_values.reinit (cell->neighbor(3));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R21(1,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R21(1,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                }
                
                LU_R21[c].initialize(A_R21); 
                
                // =====================================================================
                // r = 2 stencil 2
                // =====================================================================
                
                if (p3 == -1) { // EAST 
                    // P2 cell 

                    A_R22 = 0.0; 
                    
                    V_neighbor = cell->neighbor(3)->measure(); 
                    fv_values.reinit (cell->neighbor(3));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R22(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R22(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    A_R22(1,0) = nx1;       A_R22(1,1) = ny1;       // Transmissive boundary
                }
                
                else if (p2 == -1) { // NORTH 
                    // P1 cell 

                    A_R22 = 0.0; 
                    
                    V_neighbor = cell->neighbor(1)->measure(); 
                    fv_values.reinit (cell->neighbor(1));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R22(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R22(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    A_R22(1,0) = nx1;       A_R22(1,1) = ny1;       // Transmissive boundary
                }
                
                else { // WEST and SOUTH 
                    // P1 cell 

                    A_R22 = 0.0; 
                    
                    V_neighbor = cell->neighbor(3)->measure(); 
                    fv_values.reinit (cell->neighbor(3));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R22(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R22(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    // P2 cell 
                    
                    V_neighbor = cell->neighbor(1)->measure(); 
                    fv_values.reinit (cell->neighbor(1));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R22(1,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R22(1,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                }
                
                LU_R22[c].initialize(A_R22);
                
                // =====================================================================
                // r = 2 stencil 3
                // =====================================================================
                
                if (p3 == -1) { // EAST 
                    // P4 cell 

                    A_R23 = 0.0; 
                    
                    V_neighbor = cell->neighbor(2)->measure(); 
                    fv_values.reinit (cell->neighbor(2));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R23(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R23(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    A_R23(1,0) = nx1;       A_R23(1,1) = ny1;       // Transmissive boundary
                }
                
                else if (p4 == -1) { // SOUTH 
                    // P3 cell 

                    A_R23 = 0.0; 
                    
                    V_neighbor = cell->neighbor(1)->measure(); 
                    fv_values.reinit (cell->neighbor(1));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R23(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R23(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    A_R23(1,0) = nx1;       A_R23(1,1) = ny1;       // Transmissive boundary
                }
                
                else { // WEST and NORTH 
                    
                    // P3 cell 

                    A_R23 = 0.0; 
                    
                    V_neighbor = cell->neighbor(1)->measure(); 
                    fv_values.reinit (cell->neighbor(1));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R23(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R23(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    // P4 cell
                    
                    V_neighbor = cell->neighbor(2)->measure(); 
                    fv_values.reinit (cell->neighbor(2));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R23(1,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R23(1,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                }
                
                LU_R23[c].initialize(A_R23); 
                
                // =====================================================================
                // r = 2 stencil 4
                // =====================================================================
                
                if (p1 == -1) { // WEST
                    // P4 cell 

                    A_R24 = 0.0; 
                    
                    V_neighbor = cell->neighbor(2)->measure(); 
                    fv_values.reinit (cell->neighbor(2));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R24(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R24(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    A_R24(1,0) = nx1;       A_R24(1,1) = ny1;       // Transmissive boundary
                }
                
                else if (p4 == -1) { // SOUTH 
                    // P1 cell 

                    A_R24 = 0.0; 
                    
                    V_neighbor = cell->neighbor(0)->measure(); 
                    fv_values.reinit (cell->neighbor(0));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R24(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R24(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    A_R24(1,0) = nx1;       A_R24(1,1) = ny1;       // Transmissive boundary
                }
                
                else { // EAST and NORTH 
                    
                    // P4 cell 

                    A_R24 = 0.0; 
                    
                    V_neighbor = cell->neighbor(2)->measure(); 
                    fv_values.reinit (cell->neighbor(2));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R24(0,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R24(0,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                    // P1 cell
                    
                    V_neighbor = cell->neighbor(0)->measure(); 
                    fv_values.reinit (cell->neighbor(0));
                    
                    for (unsigned int i = 0; i < N_gp*N_gp; i++) {
                        q_point = fv_values.quadrature_point(i);
                        A_R24(1,0) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(0) - WENO_poly_consts[c](0));
                        A_R24(1,1) += (1./V_neighbor)*fv_values.JxW (i)*(q_point(1) - WENO_poly_consts[c](1));
                    }
                    
                }
                
                LU_R24[c].initialize(A_R24); 
                
            } // End of non-corner cell loop 
                
        } // End of boundary cell loop
        
        
    } // End of cell loop 
    
    std::cout << "Done!" << std::endl;
    std::cout << "===========================" << std::endl;

} // End of function  

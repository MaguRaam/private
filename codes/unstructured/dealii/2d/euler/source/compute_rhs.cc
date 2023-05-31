#include "../include/Euler.h"


// Compute the rhs vectors

void Euler_2D::compute_rhs() {

    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;
	QGauss<2-1> face_quadrature_formula(1);
	FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);
	Tensor<1,2> face_normal_vector; // Face normal vector

	// Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
	Point<2> face_center;

    double V_c, S_f; // Volume of the cell and surface area of the face
	
	double nx, ny;   // Face normal vectors
	Vector<double> UL(4); Vector<double> UR(4); // Solving the Riemann Problem
	Vector<double> WL(4); Vector<double> WR(4); // Solving the Riemann Problem
	bool boundary; 
    Vector<double> F(4);

	for (unsigned int c = 0; cell != endc; ++cell, ++c) {
    
		V_c = cell->measure();

        rhs1(c) = 0.0;
        rhs2(c) = 0.0;
        rhs3(c) = 0.0;
        rhs4(c) = 0.0;
        
        
        for (unsigned int f = 0; f < faces_per_cell; ++f) {
		
			face_center = cell->face(f)->center(); 
			S_f = cell->face(f)->measure();
            
			// Get some geometry info
    
			fv_face_values.reinit(cell, f);
                    
			face_normal_vector = fv_face_values.normal_vector(0); 
			nx = face_normal_vector[0]; ny = face_normal_vector[1];
                
			// Left face
			UL(0) = RHO[c];
			UL(1) = RHO_U[c];
			UL(2) = RHO_V[c];
			UL(3) = E[c];
                    
			WL = conserved_to_primitive(UL);
			
			if (cell->face(f)->at_boundary()) {
				
				if (cell->face(f)->boundary_id() == 0) {
					
					// Inlet
					
					WR(0) = 1.4; 
					WR(1) = M;
					WR(2) = 0.0; 
					WR(3) = 1.0;
				}
				
				if (cell->face(f)->boundary_id() == 1) {
					
					// Outlet 
					
					WR(0) = WL(0);  
					WR(1) = WL(1);  
					WR(2) = WL(2); 
					WR(3) = WL(3);
				}
				
				else {
					
					// Reflecting 
				
					WR(0) = WL(0);  
					WR(1) = WL(1) - 2.0*WL(1)*nx*nx - 2.0*WL(2)*nx*ny; 
					WR(2) = WL(2) - 2.0*WL(1)*nx*ny - 2.0*WL(2)*ny*ny;
					WR(3) = WL(3);
				}
			}
                
			else {
                    
				// Get the right state values
				
				UR(0) = RHO[cell->neighbor_index(f)];
				UR(1) = RHO_U[cell->neighbor_index(f)];
				UR(2) = RHO_V[cell->neighbor_index(f)];
				UR(3) = E[cell->neighbor_index(f)]; 
			
				WR = conserved_to_primitive(UR);
				
				boundary = false; 
			}
                
			F = local_Lax_Friedrichs_riemann_solver(WL, WR, nx, ny, face_center, boundary);

            // Add it to the rhs vectors

            rhs1(c) += (-1.0/V_c)*(F(0)*S_f);
            rhs2(c) += (-1.0/V_c)*(F(1)*S_f);
            rhs3(c) += (-1.0/V_c)*(F(2)*S_f);
            rhs4(c) += (-1.0/V_c)*(F(3)*S_f);
		} // End of face loop 
		
		
		
    } // End of cell loop 
}

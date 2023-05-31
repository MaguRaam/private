#include "../include/Weno32.h"

double evaluate_ader_polynomial(Vector<double> coeffs, Vector<double> coeffs_time, Vector<double> consts,  Point<2> P, double t) {

    // Inputs are the coefficients of the polynomial, the center of the cell and the point of evaluation

    double x = P(0); double y = P(1); 
    
    double x0 = consts(0); double y0 = consts(1); 
    
    // Use Honers Algorith 
    
    double retval = coeffs(0) + 
                    coeffs(1)*(x-x0) + 
                    coeffs(2)*(y-y0) + 
                    coeffs(3)*((x-x0)*(x-x0) - consts(2)) + 
                    coeffs(4)*((y-y0)*(y-y0) - consts(3)) +
                    coeffs(5)*((x-x0)*(y-y0) - consts(4)) + 
                    coeffs_time(0)*t +
                    coeffs_time(1)*t*t +
                    coeffs_time(2)*(x-consts(0))*t + 
                    coeffs_time(3)*(y-consts(1))*t ; 
    
    return retval; 

}

FullMatrix<double> Assemble_Space_Time_Galerkin_LHS_Matrix(Vector<double> Int) {
    
    FullMatrix<double> A(4,4);

    double r2_3 = 2./3.; double r1_3 = 1./3.; 

    A(0,0) = 0.5*Int[0];     A(0,1) = r2_3*Int[0];   A(0,2) = 0.0;            A(0,3) = 0.0;
    A(1,0) = r1_3*Int[0];    A(1,1) = 0.5*Int[0];    A(1,2) = 0.0;            A(1,3) = 0.0;
    A(2,0) = 0.0;            A(2,1) = 0.0;           A(2,2) = 0.5*Int[5];     A(2,3) = 0.5*Int[6];
    A(3,0) = 0.0;            A(3,1) = 0.0;           A(3,2) = 0.5*Int[6];     A(3,3) = 0.5*Int[9];
    
    
    return A; 

}

Vector<double> Assemble_Space_Time_Galerkin_RHS_Vector(Vector<double> F_space, Vector<double> F_time, 
                                                       Vector<double> G_space, Vector<double> G_time, 
                                                       Vector<double> Int) {
    double r1_3 = 1./3.;   double r2_3 = 2./3.;  

    Vector<double> b(4);

    b(0) = -(0.5*F_space[1]*Int[0] + F_space[3]*Int[1] + 0.5*F_space[5]*Int[2] + r1_3*F_time[2]*Int[0] + 
             0.5*G_space[2]*Int[0] + G_space[4]*Int[2] + 0.5*G_space[5]*Int[1] + r1_3*G_time[3]*Int[0] ); 

    b(1) = -(r1_3*F_space[1]*Int[0] + r2_3*F_space[3]*Int[1] + r1_3*F_space[5]*Int[2] + 0.25*F_time[2]*Int[0] + 
             r1_3*G_space[2]*Int[0] + r2_3*G_space[4]*Int[2] + r1_3*G_space[5]*Int[1] + 0.25*G_time[3]*Int[0] ); 

    b(2) = -(0.5*F_space[1]*Int[3] + F_space[3]*Int[7] + 0.5*F_space[5]*Int[8] + r1_3*F_time[2]*Int[3] + 
             0.5*G_space[2]*Int[3] + G_space[4]*Int[8] + 0.5*G_space[5]*Int[7] + r1_3*G_time[3]*Int[3] );  
  
    b(3) = -(0.5*F_space[1]*Int[4] + F_space[3]*Int[10] + 0.5*F_space[5]*Int[11] + r1_3*F_time[2]*Int[4] + 
             0.5*G_space[2]*Int[4] + G_space[4]*Int[11] + 0.5*G_space[5]*Int[10] + r1_3*G_time[3]*Int[4] );  
    
    return b; 

}


void Weno3_2D::solve_ader() {
    
    // Using ADER Method

    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;
    unsigned int vertices_per_cell = GeometryInfo<2>::vertices_per_cell;

    // Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

    // Nodal flux values 
    Vector<double> F1_nodal_space(9);  Vector<double> G1_nodal_space(9);
    Vector<double> F2_nodal_space(9);  Vector<double> G2_nodal_space(9);
    Vector<double> F3_nodal_space(9);  Vector<double> G3_nodal_space(9);
    Vector<double> F4_nodal_space(9);  Vector<double> G4_nodal_space(9);

    Vector<double> F1_nodal_time(6);  Vector<double> G1_nodal_time(6);
    Vector<double> F2_nodal_time(6);  Vector<double> G2_nodal_time(6);
    Vector<double> F3_nodal_time(6);  Vector<double> G3_nodal_time(6);
    Vector<double> F4_nodal_time(6);  Vector<double> G4_nodal_time(6);

    // Modal flux coefficients 
    Vector<double> F1_coeffs_space(6); Vector<double> F2_coeffs_space(6); 
    Vector<double> F3_coeffs_space(6); Vector<double> F4_coeffs_space(6); 
    Vector<double> G1_coeffs_space(6); Vector<double> G2_coeffs_space(6); 
    Vector<double> G3_coeffs_space(6); Vector<double> G4_coeffs_space(6);

    Vector<double> F1_coeffs_time(4); Vector<double> F2_coeffs_time(4); 
    Vector<double> F3_coeffs_time(4); Vector<double> F4_coeffs_time(4); 
    Vector<double> G1_coeffs_time(4); Vector<double> G2_coeffs_time(4); 
    Vector<double> G3_coeffs_time(4); Vector<double> G4_coeffs_time(4);

    // Time variation 
    std::vector < Vector<double> > coeffs_RHO_time(triangulation.n_active_cells()); 
    std::vector < Vector<double> > coeffs_RHO_U_time(triangulation.n_active_cells()); 
    std::vector < Vector<double> > coeffs_RHO_V_time(triangulation.n_active_cells());
    std::vector < Vector<double> > coeffs_E_time(triangulation.n_active_cells());
    std::vector < Vector <double> > Integration_Constants(triangulation.n_active_cells());
   
    for (unsigned int c = 0; c < triangulation.n_active_cells(); c++) {
        coeffs_RHO_time[c].reinit(4);
        coeffs_RHO_U_time[c].reinit(4); 
        coeffs_RHO_V_time[c].reinit(4);
        coeffs_E_time[c].reinit(4);
        Integration_Constants[c].reinit(12); 
    }

    Point<2> face_center; 
    Point<2> cell_center; 
    Point<2> vertex; 
    

    Vector<double> U(4); Vector<double> W(4);                                 // Getting the nodal values
    std::vector< Point<2> > node_list_space(9); 
      
    
    // =====================  Integration constants for space-time projections =====================
    
    unsigned int N_gp = 3;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    Point<2> q_point;

    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        

		double x, y; 

		fv_values.reinit(cell);
		
		double x0 = WENO_poly_consts[c](0); 
		double y0 = WENO_poly_consts[c](1);

		Integration_Constants[c](0) = cell->measure(); Integration_Constants[c](8) = 0.0;
		Integration_Constants[c](1) = 0.0;             Integration_Constants[c](6) = 0.0;
		Integration_Constants[c](2) = 0.0;             Integration_Constants[c](7) = 0.0;
		Integration_Constants[c](3) = 0.0;             Integration_Constants[c](9) = 0.0;
		Integration_Constants[c](4) = 0.0;             Integration_Constants[c](10) = 0.0;
		Integration_Constants[c](5) = 0.0;             Integration_Constants[c](11) = 0.0;
		

		
		for (unsigned int i = 0; i < N_gp*N_gp; i++) {

			q_point = fv_values.quadrature_point(i);

			x = q_point(0); y = q_point(1);  

		
			Integration_Constants[c](1) += (fv_values.JxW(i)*(x - x0));
			
			Integration_Constants[c](2) += (fv_values.JxW(i)*(y - y0));
			
			Integration_Constants[c](3) += (fv_values.JxW(i)*(x - x0));
			
			Integration_Constants[c](4) += (fv_values.JxW(i)*(y - y0));
			
			Integration_Constants[c](5) += (fv_values.JxW(i)*(x - x0)*(x - x0)); 
			
			Integration_Constants[c](6) += (fv_values.JxW(i)*(x - x0)*(y - y0));
			
			Integration_Constants[c](7) += (fv_values.JxW(i)*(x - x0)*(x - x0));
			
			Integration_Constants[c](8) += (fv_values.JxW(i)*(y - y0)*(x - x0));
			
			Integration_Constants[c](9) += (fv_values.JxW(i)*(y - y0)*(y - y0));
			
			Integration_Constants[c](10) += (fv_values.JxW(i)*(x - x0)*(y - y0));      
			
			Integration_Constants[c](11) += (fv_values.JxW(i)*(y - y0)*(y - y0));
			
		}
            
        
    }
    
    double x0, y0; 
    // =====================  Create the LHS Matrices for nodal to modal transcription =====================

    FullMatrix<double> A_space(9, 6);
    FullMatrix<double> A_time(6, 4);
    FullMatrix<double> A_space_time(4,4); 

    std::vector< Householder<double> > H_space(triangulation.n_active_cells());
    std::vector< Householder<double> > H_time(triangulation.n_active_cells());
    std::vector< LUdcmp > LU_space_time(triangulation.n_active_cells());
    
    // Assemble the purely spatial part of the matrix 
    
    cell = dof_handler.begin_active();

    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        cell_center = cell->center(); 
		
		x0 = WENO_poly_consts[c](0); y0 = WENO_poly_consts[c](1); 
		// Row 1 (Corresponds to cell centroid)

		A_space(0,0) = 1.0; 
		A_space(0,1) = 0.0; 
		A_space(0,2) = 0.0;
		A_space(0,3) = - WENO_poly_consts[c](2); 
		A_space(0,4) = - WENO_poly_consts[c](3);
		A_space(0,5) = - WENO_poly_consts[c](4);  

		// Row 2-5 (Corresponds to face centers)     

		for (unsigned int f = 0; f < faces_per_cell; f++) {
		
			face_center =  cell->face(f)->center(); 
	
			// Fill the matrix 
			
			A_space(f+1,0) = 1.0; 
			A_space(f+1,1) = face_center(0) - WENO_poly_consts[c](0); 
			A_space(f+1,2) = face_center(1) - WENO_poly_consts[c](1);
			A_space(f+1,3) = (face_center(0)-x0)*(face_center(0)-x0) - WENO_poly_consts[c](2); 
			A_space(f+1,4) = (face_center(1)-y0)*(face_center(1)-y0) - WENO_poly_consts[c](3);
			A_space(f+1,5) = (face_center(0)-x0)*(face_center(1)-y0) - WENO_poly_consts[c](4);  
		}

		// Row 6-9 (Corresponds to vertices)     

		for (unsigned int v = 0; v < vertices_per_cell; v++) {
		
			vertex =  cell->vertex(v); 
			
			// Fill the matrix 
			
			A_space(v+5,0) = 1.0; 
			A_space(v+5,1) = vertex(0) - WENO_poly_consts[c](0); 
			A_space(v+5,2) = vertex(1) - WENO_poly_consts[c](1);
			A_space(v+5,3) = (vertex(0)-x0)*(vertex(0)-x0) - WENO_poly_consts[c](2); 
			A_space(v+5,4) = (vertex(1)-y0)*(vertex(1)-y0) - WENO_poly_consts[c](3);
			A_space(v+5,5) = (vertex(0)-x0)*(vertex(1)-y0) - WENO_poly_consts[c](4); 
		}

		// Initialize the householder object 
		H_space[c].initialize(A_space);
        
    
    } // Finished creating LHS matrices for spatial nodes 
    
    // Assemble the purely temporal part of the matrix 
    
    cell = dof_handler.begin_active();

    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        
		x0 = WENO_poly_consts[c](0); y0 = WENO_poly_consts[c](1); 
		
		// Asseble the matrix (Row 1 - cell center)
		A_time(0,0) = 0.5; 
		A_time(0,1) = 0.5*0.5; 
		A_time(0,2) = 0.5*(x0 - WENO_poly_consts[c](0));
		A_time(0,3) = 0.5*(y0 - WENO_poly_consts[c](1));

		// Row 2-5 - face-centers

		for (unsigned int f = 0; f < faces_per_cell; f++) {

			face_center =  cell->face(f)->center(); 

			A_time(f+1,0) = 0.5; 
			A_time(f+1,1) = 0.5*0.5; 
			A_time(f+1,2) = 0.5*(face_center(0) - WENO_poly_consts[c](0));
			A_time(f+1,3) = 0.5*(face_center(1) - WENO_poly_consts[c](1));
		}
		
		// Row 6 - cell center, tau = 1.0
		A_time(5,0) = 1.0; 
		A_time(5,1) = 1.0; 
		A_time(5,2) = (x0 - WENO_poly_consts[c](0));
		A_time(5,3) = (y0 - WENO_poly_consts[c](1));
		
		H_time[c].initialize(A_time); 
        
    
    } // Finished creating LHS matrices for temporal nodes
    
    // Assemble the space time matrix
    
	cell = dof_handler.begin_active();

    for (unsigned int c = 0; cell != endc; ++cell, ++c) {

		A_space_time = Assemble_Space_Time_Galerkin_LHS_Matrix(Integration_Constants[c]);
		LU_space_time[c].initialize(A_space_time); 
    
    } // Finished creating space time matrices 
    
    

    // ===================== Variables required for updating in time =====================

    Vector<double> RHO_rhs_time(4);  Vector<double> RHO_U_rhs_time(4); Vector<double> RHO_V_rhs_time(4); Vector<double> E_rhs_time(4); 
    double time = 0.0;
    unsigned int count = 0;  
    unsigned int ITER = 0; unsigned MAX_ITER = 2; 
    
    double t_gauss[] = {-1.0/(2.0*std::sqrt(3)) + 0.5, 1.0/(2.0*std::sqrt(3)) + 0.5};  // Quadrature points in [ 0.0, 1.0] 

    // ===================== Variables required for corrector step =====================    

	double V_c, S_f; // Volume of the cell and surface area of the face
	
	double nx1, ny1;   // Face normal vectors
	double nx2, ny2; 

    // For obtaining face normal vectors
	QGauss<2-1> face_quadrature_formula(2);
	FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);

    Point<2> face_quadrature_point_1;
    Point<2> face_quadrature_point_2;
	
	Tensor<1,2> face_normal_vector1; // Face normal vector
	Tensor<1,2> face_normal_vector2; // Face normal vector
	
    Vector<double> UL1(4); Vector<double> UR1(4); // Solving the Riemann Problem
    Vector<double> UL2(4); Vector<double> UR2(4); // Solving the Riemann Problem
    Vector<double> UL3(4); Vector<double> UR3(4); // Solving the Riemann Problem
	Vector<double> UL4(4); Vector<double> UR4(4); // Solving the Riemann Problem
    
    Vector<double> WL1(4); Vector<double> WR1(4); // Solving the Riemann Problem
    Vector<double> WL2(4); Vector<double> WR2(4); // Solving the Riemann Problem
    Vector<double> WL3(4); Vector<double> WR3(4); // Solving the Riemann Problem
    Vector<double> WL4(4); Vector<double> WR4(4); // Solving the Riemann Problem
	
    Vector<double> F1(4); Vector<double> F2(4); // Flux at the face
    Vector<double> F3(4); Vector<double> F4(4); // Flux at the face
    Vector<double> F(4);
    
    bool boundary; unsigned int global_face_index;   
    
    unsigned int n_faces = triangulation.n_active_faces(); 

    std::vector< Vector<double> > Flux1(n_faces); 
    std::vector< Vector<double> > Flux2(n_faces);
    std::vector< Vector<double> > Flux3(n_faces);
    std::vector< Vector<double> > Flux4(n_faces);
    std::vector< bool > did_not_compute_flux_for_the_face(n_faces); 
    
    for (unsigned int f = 0; f < n_faces; f++) {
        Flux1[f].reinit(4);
        Flux2[f].reinit(4);
        Flux3[f].reinit(4);
        Flux4[f].reinit(4);
        did_not_compute_flux_for_the_face[f] = true; 
    }

    while (time < finalTime) {
		
		compute_time_step_based_on_cfl_number(time);
		time += dt;
		
		if (count%10 == 0) {
			restart(time, count); 
			post_process(count); 
            output_results(count);
        }
        
        std::cout << "time = " << time << ", Final time = " << finalTime<< std::endl;
        count++;

        // Get the coefficients 
        reconstruct(); 

        cell = dof_handler.begin_active();

        // =================================================  Predictor Step  =================================================
        
        for (unsigned int c = 0; cell != endc; ++cell, ++c) { 
            

			// Get the list of nodes 
		
			cell_center = cell->center();
			node_list_space[0] = cell_center;  

			// Row 2-5 (Corresponds to face centers)     

			for (unsigned int f = 0; f < faces_per_cell; f++) {
				node_list_space[f+1] = cell->face(f)->center();
			}

			// Row 6-9 (Corresponds to vertices)     

			for (unsigned int v = 0; v < vertices_per_cell; v++) { 
				node_list_space[v+5] = cell->vertex(v); 
			} 

			for (unsigned int i = 0; i < 9; i++) {
				
				// Evaluate Nodal values of Conserved Variables 
				U(0) = evaluate_weno_polynomial(coeffs_RHO[c],   WENO_poly_consts[c], node_list_space[i]);
				U(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], node_list_space[i]);
				U(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], node_list_space[i]);
				U(3) = evaluate_weno_polynomial(coeffs_E[c],     WENO_poly_consts[c], node_list_space[i]);

				// Convert to primitive variables 
				W = conserved_to_primitive(U);
		
				// Fill the vector
				F1_nodal_space(i) = dt*W(0)*W(1);                    G1_nodal_space(i) = dt*W(0)*W(2);                
				F2_nodal_space(i) = dt*(W(0)*W(1)*W(1) + W(3));      G2_nodal_space(i) = dt*W(0)*W(1)*W(2);     
				F3_nodal_space(i) = dt*W(0)*W(1)*W(2);               G3_nodal_space(i) = dt*(W(0)*W(2)*W(2) + W(3));        
				F4_nodal_space(i) = dt*W(1)*(U(3) + W(3));           G4_nodal_space(i) = dt*W(2)*(U(3) + W(3)); 
			}
				

			// Solve the least squares system to get flux coefficients (in space)

			H_space[c].least_squares(F1_coeffs_space, F1_nodal_space);  
			H_space[c].least_squares(F2_coeffs_space, F2_nodal_space); 
			H_space[c].least_squares(F3_coeffs_space, F3_nodal_space); 
			H_space[c].least_squares(F4_coeffs_space, F4_nodal_space);
			
			H_space[c].least_squares(G1_coeffs_space, G1_nodal_space);  
			H_space[c].least_squares(G2_coeffs_space, G2_nodal_space); 
			H_space[c].least_squares(G3_coeffs_space, G3_nodal_space); 
			H_space[c].least_squares(G4_coeffs_space, G4_nodal_space);

			// Begin Iterations to get time-dependent modes 

			ITER = 0;  

			coeffs_RHO_time[c] = 0.0; coeffs_RHO_U_time[c] = 0.0; coeffs_RHO_V_time[c] = 0.0; coeffs_E_time[c] = 0.0; 

			// Update nodal space values to match the flux coefficients 

			for (unsigned int i = 0; i < 6; i++) {
				
				face_center = node_list_space[i]; 
				
				if (i == 5) {
					face_center = node_list_space[0]; 
				}
			
				F1_nodal_space(i) = evaluate_weno_polynomial(F1_coeffs_space, WENO_poly_consts[c], face_center);  
				F2_nodal_space(i) = evaluate_weno_polynomial(F2_coeffs_space, WENO_poly_consts[c], face_center);
				F3_nodal_space(i) = evaluate_weno_polynomial(F3_coeffs_space, WENO_poly_consts[c], face_center);
				F4_nodal_space(i) = evaluate_weno_polynomial(F4_coeffs_space, WENO_poly_consts[c], face_center); 

				G1_nodal_space(i) = evaluate_weno_polynomial(G1_coeffs_space, WENO_poly_consts[c], face_center);  
				G2_nodal_space(i) = evaluate_weno_polynomial(G2_coeffs_space, WENO_poly_consts[c], face_center);
				G3_nodal_space(i) = evaluate_weno_polynomial(G3_coeffs_space, WENO_poly_consts[c], face_center);
				G4_nodal_space(i) = evaluate_weno_polynomial(G4_coeffs_space, WENO_poly_consts[c], face_center); 
						
			}  

			while (ITER < MAX_ITER) {
                                    
				// Get the time dependent nodal values of flux at tau = 0.5*dt

				for (unsigned int i = 0; i < 5; i++) {
					
					// Evaluate Nodal values of Conserved Variables 
					U(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], node_list_space[i], 0.5);
					U(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], node_list_space[i], 0.5);
					U(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], node_list_space[i], 0.5);
					U(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], node_list_space[i], 0.5);

					// Convert to primitive variables 
					W = conserved_to_primitive(U);

					// Evaluate corresponding nodal fluxes
					F1_nodal_time(i) = dt*W(0)*W(1);                    G1_nodal_time(i) = dt*W(0)*W(2);                
					F2_nodal_time(i) = dt*(W(0)*W(1)*W(1) + W(3));      G2_nodal_time(i) = dt*W(0)*W(1)*W(2);     
					F3_nodal_time(i) = dt*W(0)*W(1)*W(2);               G3_nodal_time(i) = dt*(W(0)*W(2)*W(2) + W(3));        
					F4_nodal_time(i) = dt*W(1)*(U(3) + W(3));           G4_nodal_time(i) = dt*W(2)*(U(3) + W(3)); 
				}

				// Get the time dependent nodal values of flux at tau = dt

				// Evaluate Nodal values of Conserved Variables 
				U(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], cell_center, 1.0);
				U(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], cell_center, 1.0);
				U(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], cell_center, 1.0);
				U(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], cell_center, 1.0);

				// Convert to primitive variables 
				W = conserved_to_primitive(U);

				// Evaluate corresponding nodal fluxes
				F1_nodal_time(5) = dt*W(0)*W(1);                    G1_nodal_time(5) = dt*W(0)*W(2);                
				F2_nodal_time(5) = dt*(W(0)*W(1)*W(1) + W(3));      G2_nodal_time(5) = dt*W(0)*W(1)*W(2);     
				F3_nodal_time(5) = dt*W(0)*W(1)*W(2);               G3_nodal_time(5) = dt*(W(0)*W(2)*W(2) + W(3));        
				F4_nodal_time(5) = dt*W(1)*(U(3) + W(3));           G4_nodal_time(5) = dt*W(2)*(U(3) + W(3)); 

				// Solve the least squares system to get flux coefficients (in time)

				// Subtract the pure spatial part from the time nodal values to solve the least squares system 

				for (unsigned int i = 0; i < 6; i++) {
					
					F1_nodal_time(i) = F1_nodal_time(i) - F1_nodal_space(i); 
					F2_nodal_time(i) = F2_nodal_time(i) - F2_nodal_space(i);
					F3_nodal_time(i) = F3_nodal_time(i) - F3_nodal_space(i);
					F4_nodal_time(i) = F4_nodal_time(i) - F4_nodal_space(i);

					G1_nodal_time(i) = G1_nodal_time(i) - G1_nodal_space(i);  
					G2_nodal_time(i) = G2_nodal_time(i) - G2_nodal_space(i);
					G3_nodal_time(i) = G3_nodal_time(i) - G3_nodal_space(i);
					G4_nodal_time(i) = G4_nodal_time(i) - G4_nodal_space(i);
				}  

				H_time[c].least_squares(F1_coeffs_time, F1_nodal_time);  
				H_time[c].least_squares(F2_coeffs_time, F2_nodal_time); 
				H_time[c].least_squares(F3_coeffs_time, F3_nodal_time); 
				H_time[c].least_squares(F4_coeffs_time, F4_nodal_time);
				
				H_time[c].least_squares(G1_coeffs_time, G1_nodal_time);  
				H_time[c].least_squares(G2_coeffs_time, G2_nodal_time); 
				H_time[c].least_squares(G3_coeffs_time, G3_nodal_time); 
				H_time[c].least_squares(G4_coeffs_time, G4_nodal_time);  
				
				// Take space-time Galerkin projections to imporve the solutions  


				// Get RHS Vectors 

				RHO_rhs_time   = Assemble_Space_Time_Galerkin_RHS_Vector(F1_coeffs_space, F1_coeffs_time, G1_coeffs_space, G1_coeffs_time, 
																		Integration_Constants[c]); 
				
				RHO_U_rhs_time = Assemble_Space_Time_Galerkin_RHS_Vector(F2_coeffs_space, F2_coeffs_time, G2_coeffs_space, G2_coeffs_time, 
																		Integration_Constants[c]);
				
				RHO_V_rhs_time = Assemble_Space_Time_Galerkin_RHS_Vector(F3_coeffs_space, F3_coeffs_time, G3_coeffs_space, G3_coeffs_time, 
																		Integration_Constants[c]);
			
				E_rhs_time     = Assemble_Space_Time_Galerkin_RHS_Vector(F4_coeffs_space, F4_coeffs_time, G4_coeffs_space, G4_coeffs_time, 
																		Integration_Constants[c]);   
				
				// Update time-dependent coefficients 

				LU_space_time[c].solve(RHO_rhs_time, coeffs_RHO_time[c]);        
				LU_space_time[c].solve(RHO_U_rhs_time, coeffs_RHO_U_time[c]);
				LU_space_time[c].solve(RHO_V_rhs_time, coeffs_RHO_V_time[c]);
				LU_space_time[c].solve(E_rhs_time, coeffs_E_time[c]);

				ITER++; 
						
			} // End of iterations
            
        } // End of cell loop    

        
        // =================================================  Corrector Step  =================================================
        
        cell = dof_handler.begin_active();

        for (unsigned int f = 0; f < n_faces; f++) {
            did_not_compute_flux_for_the_face[f] = true; 
        }
        

    	for (unsigned int c = 0; cell != endc; ++cell, ++c) {
    
		    V_c = cell->measure();

		    rhs1(c) = 0.0;
		    rhs2(c) = 0.0;
		    rhs3(c) = 0.0;
		    rhs4(c) = 0.0;
        
        
            for (unsigned int f = 0; f < faces_per_cell; ++f) {
            
                global_face_index = cell->face_index(f);
                S_f = cell->face(f)->measure();
            
                if(did_not_compute_flux_for_the_face[global_face_index]) {
                
                    // Get some geometry info
            
                    
                    fv_face_values.reinit(cell, f);
                    
                    face_normal_vector1 = fv_face_values.normal_vector(0); // Zero corresponds to the presence of only one quadrature point
                    nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                    
                    face_normal_vector2 = fv_face_values.normal_vector(1);
                    nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
                    
                    face_quadrature_point_1 = fv_face_values.quadrature_point(0); 
                    face_quadrature_point_2 = fv_face_values.quadrature_point(1); 
                
                    // Left face
                    UL1(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0]);
                    UL1(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0]);
                    UL1(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0]);
                    UL1(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0]);
                    
                    UL2(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0]);
                    UL2(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0]);
                    UL2(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0]);
                    UL2(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0]);

                    UL3(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1]);
                    UL3(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1]);
                    UL3(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1]);
                    UL3(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1]);
                    
                    UL4(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1]);
                    UL4(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1]);
                    UL4(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1]);
                    UL4(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1]);
                    
                    WL1 = conserved_to_primitive(UL1); WL2 = conserved_to_primitive(UL2);
                    WL3 = conserved_to_primitive(UL3); WL4 = conserved_to_primitive(UL4);
                
                    if (cell->face(f)->at_boundary()) {
                        
                        // Apply transmissive boundary conditions
                        
                        WR1 = WL1; WR2 = WL2; WR3 = WL3; WR4 = WL4;
                        
                        boundary = true; 
                    }

                
                    else {
                        
                        // Get the right state values
                        
                        UR1(0) = evaluate_ader_polynomial(coeffs_RHO[cell->neighbor_index(f)], coeffs_RHO_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1, t_gauss[0]);
                        
                        UR1(1) = evaluate_ader_polynomial(coeffs_RHO_U[cell->neighbor_index(f)], coeffs_RHO_U_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1, t_gauss[0]);
                        
                        UR1(2) = evaluate_ader_polynomial(coeffs_RHO_V[cell->neighbor_index(f)], coeffs_RHO_V_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1, t_gauss[0]);
                        

                        UR1(3) = evaluate_ader_polynomial(coeffs_E[cell->neighbor_index(f)], coeffs_E_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1, t_gauss[0]);
                        
                        
                        UR2(0) = evaluate_ader_polynomial(coeffs_RHO[cell->neighbor_index(f)], coeffs_RHO_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2, t_gauss[0]);
                        
                        UR2(1) = evaluate_ader_polynomial(coeffs_RHO_U[cell->neighbor_index(f)], coeffs_RHO_U_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2, t_gauss[0]);
                        
                        UR2(2) = evaluate_ader_polynomial(coeffs_RHO_V[cell->neighbor_index(f)], coeffs_RHO_V_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2, t_gauss[0]);
                        
                        UR2(3) = evaluate_ader_polynomial(coeffs_E[cell->neighbor_index(f)], coeffs_E_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2, t_gauss[0]); 

                        UR3(0) = evaluate_ader_polynomial(coeffs_RHO[cell->neighbor_index(f)], coeffs_RHO_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1, t_gauss[1]);
                        
                        UR3(1) = evaluate_ader_polynomial(coeffs_RHO_U[cell->neighbor_index(f)], coeffs_RHO_U_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1, t_gauss[1]);
                        
                        UR3(2) = evaluate_ader_polynomial(coeffs_RHO_V[cell->neighbor_index(f)], coeffs_RHO_V_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1, t_gauss[1]);
                        
                        UR3(3) = evaluate_ader_polynomial(coeffs_E[cell->neighbor_index(f)], coeffs_E_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1, t_gauss[1]);
                        
                        UR4(0) = evaluate_ader_polynomial(coeffs_RHO[cell->neighbor_index(f)], coeffs_RHO_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2, t_gauss[1]);
                        
                        UR4(1) = evaluate_ader_polynomial(coeffs_RHO_U[cell->neighbor_index(f)], coeffs_RHO_U_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2, t_gauss[1]);
                        
                        UR4(2) = evaluate_ader_polynomial(coeffs_RHO_V[cell->neighbor_index(f)], coeffs_RHO_V_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2, t_gauss[1]);
                        
                        UR4(3) = evaluate_ader_polynomial(coeffs_E[cell->neighbor_index(f)], coeffs_E_time[cell->neighbor_index(f)],
                                                          WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2, t_gauss[1]); 
                        
                        WR1 = conserved_to_primitive(UR1); WR2 = conserved_to_primitive(UR2);
                        WR3 = conserved_to_primitive(UR3); WR4 = conserved_to_primitive(UR4);
                        
                        boundary = false; 
                    }
        
                
                    Flux1[global_face_index] = local_Lax_Friedrichs_riemann_solver(WL1, WR1, nx1, ny1, face_quadrature_point_1, boundary);
                    Flux2[global_face_index] = local_Lax_Friedrichs_riemann_solver(WL2, WR2, nx2, ny2, face_quadrature_point_2, boundary);
                    Flux3[global_face_index] = local_Lax_Friedrichs_riemann_solver(WL3, WR3, nx1, ny1, face_quadrature_point_1, boundary);
                    Flux4[global_face_index] = local_Lax_Friedrichs_riemann_solver(WL4, WR4, nx2, ny2, face_quadrature_point_2, boundary);
                
                
                    did_not_compute_flux_for_the_face[global_face_index] = false; 
                }   
            
                else {
                    
                    Flux1[global_face_index] *= -1.0; 
                    Flux2[global_face_index] *= -1.0;
                    Flux3[global_face_index] *= -1.0;
                    Flux4[global_face_index] *= -1.0;
                }
            
            
                F(0) = 0.25*(Flux1[global_face_index](0) + Flux2[global_face_index](0) + Flux3[global_face_index](0) + Flux4[global_face_index](0)); 
                F(1) = 0.25*(Flux1[global_face_index](1) + Flux2[global_face_index](1) + Flux3[global_face_index](1) + Flux4[global_face_index](1)); 
                F(2) = 0.25*(Flux1[global_face_index](2) + Flux2[global_face_index](2) + Flux3[global_face_index](2) + Flux4[global_face_index](2)); 
                F(3) = 0.25*(Flux1[global_face_index](3) + Flux2[global_face_index](3) + Flux3[global_face_index](3) + Flux4[global_face_index](3)); 

                // Add it to the rhs vectors

                rhs1(c) += (-1.0/V_c)*(F(0)*S_f);
                rhs2(c) += (-1.0/V_c)*(F(1)*S_f);
                rhs3(c) += (-1.0/V_c)*(F(2)*S_f);
                rhs4(c) += (-1.0/V_c)*(F(3)*S_f);
            }
        }  // End of cell loop  

        // Update the variables  

		for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
			RHO(c) = RHO(c) + dt*rhs1(c);
			RHO_U(c) = RHO_U(c) + dt*rhs2(c);
			RHO_V(c) = RHO_V(c) + dt*rhs3(c);
			E(c) = E(c) + dt*rhs4(c);
		} 
		
    } // End of time loop 
    
    restart(time, count); 
    post_process(count);
    output_results(count); 
} // End of function 
 

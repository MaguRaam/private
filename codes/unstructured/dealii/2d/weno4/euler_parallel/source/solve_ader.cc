#include "../include/Weno432.h"

// Update the solution using ADER method

double evaluate_ader_polynomial(Vector<double>, Vector<double>, Vector<double>,  Point<2>, double, double);
FullMatrix<double> Assemble_Space_Time_Galerkin_LHS_Matrix(Vector<double>);
void Assemble_Space_Time_Galerkin_RHS_Vector(Vector<double>, Vector<double>, Vector<double>, Vector<double>, Vector<double>, Vector<double>&); 

void Weno4_2D::solve_ader() {

	Use_ader = true;

    auto start = std::chrono::system_clock::now();
    
    double r1_3 = 1./3.; 
    double r2_3 = 2./3.; 
   
    double h, h_sq; 

	double dx_h, dy_h, x0, y0;
	
	unsigned int g_i; 
    
    // Using ADER Method

    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;
    unsigned int vertices_per_cell = GeometryInfo<2>::vertices_per_cell;

    // Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell;

    // Nodal flux values 
    Vector<double> F1_nodal_space(13);  Vector<double> G1_nodal_space(13);
    Vector<double> F2_nodal_space(13);  Vector<double> G2_nodal_space(13);
    Vector<double> F3_nodal_space(13);  Vector<double> G3_nodal_space(13);
    Vector<double> F4_nodal_space(13);  Vector<double> G4_nodal_space(13);

    Vector<double> F1_nodal_time(15);  Vector<double> G1_nodal_time(15);
    Vector<double> F2_nodal_time(15);  Vector<double> G2_nodal_time(15);
    Vector<double> F3_nodal_time(15);  Vector<double> G3_nodal_time(15);
    Vector<double> F4_nodal_time(15);  Vector<double> G4_nodal_time(15);

    // Modal flux coefficients 
    Vector<double> F1_coeffs_space(10); Vector<double> F2_coeffs_space(10); 
    Vector<double> F3_coeffs_space(10); Vector<double> F4_coeffs_space(10); 
    Vector<double> G1_coeffs_space(10); Vector<double> G2_coeffs_space(10); 
    Vector<double> G3_coeffs_space(10); Vector<double> G4_coeffs_space(10);

    Vector<double> F1_coeffs_time(10); Vector<double> F2_coeffs_time(10); 
    Vector<double> F3_coeffs_time(10); Vector<double> F4_coeffs_time(10); 
    Vector<double> G1_coeffs_time(10); Vector<double> G2_coeffs_time(10); 
    Vector<double> G3_coeffs_time(10); Vector<double> G4_coeffs_time(10);

    // Time variation 
    std::vector < Vector<double> > coeffs_RHO_time(n_relevant_cells); 
    std::vector < Vector<double> > coeffs_RHO_U_time(n_relevant_cells); 
    std::vector < Vector<double> > coeffs_RHO_V_time(n_relevant_cells);
    std::vector < Vector<double> > coeffs_E_time(n_relevant_cells);
    std::vector < Vector <double> > Integration_Constants(n_relevant_cells);
   
    for (unsigned int c = 0; c < n_relevant_cells; c++) {
        coeffs_RHO_time[c].reinit(10);
        coeffs_RHO_U_time[c].reinit(10); 
        coeffs_RHO_V_time[c].reinit(10);
        coeffs_E_time[c].reinit(10);
        Integration_Constants[c].reinit(51); 
    }

    Point<2> face_center; 
    Point<2> cell_center; 
    Point<2> vertex; 
    Point<2> P; 

    Vector<double> U(4); Vector<double> W(4);                                 // Getting the nodal values
    std::vector< Point<2> > node_list_space(13); 
      
    
    // =====================  Integration constants for space-time projections =====================

    
    unsigned int N_gp = 5;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    Point<2> q_point;

	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		cell = local_index_to_iterator[c];
		
		h_sq = Cell[c].measure(); 
		h = std::sqrt(h_sq); 

		double x, y; 

		fv_values.reinit(cell);
		
		x0 = WENO_poly_consts[c](0); 
		y0 = WENO_poly_consts[c](1);

		Integration_Constants[c](0) = h_sq;            Integration_Constants[c](8) = 0.0;
		Integration_Constants[c](1) = 0.0;             Integration_Constants[c](6) = 0.0;
		Integration_Constants[c](2) = 0.0;             Integration_Constants[c](7) = 0.0;
		Integration_Constants[c](3) = 0.0;             Integration_Constants[c](9) = 0.0;
		Integration_Constants[c](4) = 0.0;             Integration_Constants[c](10) = 0.0;
		Integration_Constants[c](5) = 0.0;             Integration_Constants[c](11) = 0.0;
		Integration_Constants[c](24) = 0.0; 
		
		
		for (unsigned int i = 0; i < N_gp*N_gp; i++) {

			q_point = fv_values.quadrature_point(i);

			x = q_point(0); y = q_point(1); 
			
			dx_h = (x - x0)/h; 
			dy_h = (y - y0)/h;

			Integration_Constants[c](1) += (fv_values.JxW(i)*dx_h);
			
			Integration_Constants[c](2) += (fv_values.JxW(i)*dy_h);
			
			Integration_Constants[c](3) += (fv_values.JxW(i)*dx_h);
			
			Integration_Constants[c](4) += (fv_values.JxW(i)*dy_h);
			
			Integration_Constants[c](5) += (fv_values.JxW(i)*dx_h*dx_h); 
			
			Integration_Constants[c](6) += (fv_values.JxW(i)*dx_h*dy_h);
			
			Integration_Constants[c](7) += (fv_values.JxW(i)*dx_h*dx_h);
			
			Integration_Constants[c](8) += (fv_values.JxW(i)*dy_h*dx_h);
			
			Integration_Constants[c](9) += (fv_values.JxW(i)*dy_h*dy_h);
			
			Integration_Constants[c](10) += (fv_values.JxW(i)*dx_h*dy_h);      
			
			Integration_Constants[c](11) += (fv_values.JxW(i)*dy_h*dy_h);
			
			Integration_Constants[c](12) += fv_values.JxW (i)*(dx_h*dx_h - WENO_poly_consts[c](2));
			
			Integration_Constants[c](13) += fv_values.JxW (i)*(dy_h*dy_h - WENO_poly_consts[c](3));
			
			Integration_Constants[c](14) += fv_values.JxW (i)*(dx_h*dy_h - WENO_poly_consts[c](4));
			
			Integration_Constants[c](15) += fv_values.JxW (i)*dx_h*(dx_h*dx_h - WENO_poly_consts[c](2));
			
			Integration_Constants[c](16) += fv_values.JxW (i)*dx_h*(dy_h*dy_h - WENO_poly_consts[c](3));
			
			Integration_Constants[c](17) += fv_values.JxW (i)*dx_h*(dx_h*dy_h - WENO_poly_consts[c](4));
			
			Integration_Constants[c](18) += fv_values.JxW (i)*dy_h*(dx_h*dx_h - WENO_poly_consts[c](2));
			
			Integration_Constants[c](19) += fv_values.JxW (i)*dy_h*(dy_h*dy_h - WENO_poly_consts[c](3));
			
			Integration_Constants[c](20) += fv_values.JxW (i)*dy_h*(dx_h*dy_h - WENO_poly_consts[c](4));
			
			Integration_Constants[c](21) += fv_values.JxW (i)*(dx_h*dx_h - WENO_poly_consts[c](2))*(dx_h*dx_h - WENO_poly_consts[c](2));
			
			Integration_Constants[c](22) += fv_values.JxW (i)*(dx_h*dx_h - WENO_poly_consts[c](2))*(dy_h*dy_h - WENO_poly_consts[c](3));
			
			Integration_Constants[c](23) += fv_values.JxW (i)*(dx_h*dx_h - WENO_poly_consts[c](2))*(dx_h*dy_h - WENO_poly_consts[c](4));
			
			Integration_Constants[c](24) += fv_values.JxW (i)*(dy_h*dy_h - WENO_poly_consts[c](3))*(dy_h*dy_h - WENO_poly_consts[c](3));
			
			Integration_Constants[c](25) += fv_values.JxW (i)*(dy_h*dy_h - WENO_poly_consts[c](3))*(dx_h*dy_h - WENO_poly_consts[c](4));
			
			Integration_Constants[c](26) += fv_values.JxW (i)*(dx_h*dy_h - WENO_poly_consts[c](4))*(dx_h*dy_h - WENO_poly_consts[c](4));
			
			Integration_Constants[c](27) += fv_values.JxW (i)*(dx_h*dx_h);
			
			Integration_Constants[c](28) += fv_values.JxW (i)*(dx_h*dy_h);
			
			Integration_Constants[c](29) += fv_values.JxW (i)*(dy_h*dy_h);
			
			Integration_Constants[c](30) += fv_values.JxW (i)*(dx_h*dx_h)*dx_h;
			
			Integration_Constants[c](31) += fv_values.JxW (i)*(dx_h*dy_h)*dx_h;
			
			Integration_Constants[c](32) += fv_values.JxW (i)*(dy_h*dy_h)*dx_h;
			
			Integration_Constants[c](33) += fv_values.JxW (i)*(dx_h*dx_h)*dy_h;
			
			Integration_Constants[c](34) += fv_values.JxW (i)*(dx_h*dy_h)*dy_h;
			
			Integration_Constants[c](35) += fv_values.JxW (i)*(dy_h*dy_h)*dy_h;
			
			Integration_Constants[c](36) += fv_values.JxW (i)*dx_h*(dx_h*dx_h-WENO_poly_consts[c](2));
			
			Integration_Constants[c](37) += fv_values.JxW (i)*dy_h*(dx_h*dx_h-WENO_poly_consts[c](2));
			
			Integration_Constants[c](38) += fv_values.JxW (i)*(dx_h*dx_h)*(dx_h*dx_h-WENO_poly_consts[c](2));
			
			Integration_Constants[c](39) += fv_values.JxW (i)*(dx_h*dy_h)*(dx_h*dx_h-WENO_poly_consts[c](2));
			
			Integration_Constants[c](40) += fv_values.JxW (i)*(dy_h*dy_h)*(dx_h*dx_h-WENO_poly_consts[c](2));
			
			Integration_Constants[c](41) += fv_values.JxW (i)*dx_h*(dy_h*dy_h-WENO_poly_consts[c](3));
			
			Integration_Constants[c](42) += fv_values.JxW (i)*dy_h*(dy_h*dy_h-WENO_poly_consts[c](3));
			
			Integration_Constants[c](43) += fv_values.JxW (i)*(dx_h*dx_h)*(dy_h*dy_h-WENO_poly_consts[c](3));
			
			Integration_Constants[c](44) += fv_values.JxW (i)*(dx_h*dy_h)*(dy_h*dy_h-WENO_poly_consts[c](3));
			
			Integration_Constants[c](45) += fv_values.JxW (i)*(dy_h*dy_h)*(dy_h*dy_h-WENO_poly_consts[c](3));
			
			Integration_Constants[c](46) += fv_values.JxW (i)*dx_h*(dx_h*dy_h-WENO_poly_consts[c](4));
			
			Integration_Constants[c](47) += fv_values.JxW (i)*dy_h*(dx_h*dy_h-WENO_poly_consts[c](4));
			
			Integration_Constants[c](48) += fv_values.JxW (i)*(dx_h*dx_h)*(dx_h*dy_h-WENO_poly_consts[c](4));
			
			Integration_Constants[c](49) += fv_values.JxW (i)*(dx_h*dy_h)*(dx_h*dy_h-WENO_poly_consts[c](4));
			
			Integration_Constants[c](50) += fv_values.JxW (i)*(dy_h*dy_h)*(dx_h*dy_h-WENO_poly_consts[c](4));
			
		}
		
		for (unsigned int k = 0; k < 51; ++k) {
			Integration_Constants[c](k) = Integration_Constants[c](k)/h_sq;
		}
        
    }
    // =====================  Create the LHS Matrices for nodal to modal transcription =====================

    FullMatrix<double> A(13, 10);
    FullMatrix<double> A_time(15, 10);
    FullMatrix<double> A_space_time(10, 10);

    std::vector< Householder<double> > H_space(n_relevant_cells);
    std::vector< Householder<double> > H_time(n_relevant_cells);
    std::vector< LUdcmp > LU_space_time(n_relevant_cells);
    
    cell = dof_handler.begin_active();
    
    // ===================== Variables required for updating in time =====================

    Vector<double> RHO_rhs_time(10);  Vector<double> RHO_U_rhs_time(10); Vector<double> RHO_V_rhs_time(10); Vector<double> E_rhs_time(10); 
    unsigned int count = 0;  
    unsigned int ITER = 0; unsigned MAX_ITER = 4; 
    
    double t_gauss[] = {-1.0/(2.0*std::sqrt(3)) + 0.5, 1.0/(2.0*std::sqrt(3)) + 0.5};  // Quadrature points in [ 0.0, 1.0] 

    // ===================== Variables required for corrector step =====================    

	double V_c, S_f; // Volume of the cell and surface area of the face
	
	double nx1, ny1;   // Face normal vectors
	double nx2, ny2; 

    Point<2> face_quadrature_point_1; Point<2> periodic_face_quadrature_point_1;
    Point<2> face_quadrature_point_2; Point<2> periodic_face_quadrature_point_2;
	
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
    
    bool boundary; unsigned int local_face_index;

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
    
    // =====================  Create the LHS Matrices for nodal to modal transcription in space =====================
    
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		cell = local_index_to_iterator[c];
            
		h = Cell[c].h(); 
		
		// Row 1 (Corresponds to cell centroid)
		
		x0 = WENO_poly_consts[c](0); 
		y0 = WENO_poly_consts[c](1);
		
		cell_center(0) = x0; 
		cell_center(1) = y0; 
		
		A(0,0) = 1.0; 
		A(0,1) = 0.0; 
		A(0,2) = 0.0;
		A(0,3) = - WENO_poly_consts[c](2); 
		A(0,4) = - WENO_poly_consts[c](3);
		A(0,5) = - WENO_poly_consts[c](4); 
		A(0,6) = - WENO_poly_consts[c](5); 
		A(0,7) = - WENO_poly_consts[c](6);
		A(0,8) = - WENO_poly_consts[c](7);
		A(0,9) = - WENO_poly_consts[c](8);

		// Row 2-5 (Corresponds to face centers)     

		for (unsigned int f = 0; f < faces_per_cell; f++) {
		
			face_center =  cell->face(f)->center(); 
	
			// Fill the matrix 
			
			dx_h = (face_center(0) - x0)/h; dy_h = (face_center(1) - y0)/h;
			
			A(f+1,0) = 1.0; 
			A(f+1,1) = dx_h; 
			A(f+1,2) = dy_h;
			A(f+1,3) = dx_h*dx_h - WENO_poly_consts[c](2); 
			A(f+1,4) = dy_h*dy_h - WENO_poly_consts[c](3);
			A(f+1,5) = dx_h*dy_h - WENO_poly_consts[c](4);
			A(f+1,6) = dx_h*dx_h*dx_h - WENO_poly_consts[c](5); 
			A(f+1,7) = dy_h*dy_h*dy_h - WENO_poly_consts[c](6);
			A(f+1,8) = dx_h*dx_h*dy_h - WENO_poly_consts[c](7);
			A(f+1,9) = dx_h*dy_h*dy_h - WENO_poly_consts[c](8);

		}

		// Row 6-9 (Corresponds to vertices)     

		for (unsigned int v = 0; v < vertices_per_cell; v++) {
		
			vertex =  cell->vertex(v); 
			
			dx_h = (vertex(0) - x0)/h; dy_h = (vertex(1) - y0)/h;
			
			// Fill the matrix 
			
			A(v+5,0) = 1.0; 
			A(v+5,1) = dx_h; 
			A(v+5,2) = dy_h;
			A(v+5,3) = dx_h*dx_h - WENO_poly_consts[c](2); 
			A(v+5,4) = dy_h*dy_h - WENO_poly_consts[c](3);
			A(v+5,5) = dx_h*dy_h - WENO_poly_consts[c](4);
			A(v+5,6) = dx_h*dx_h*dx_h - WENO_poly_consts[c](5); 
			A(v+5,7) = dy_h*dy_h*dy_h - WENO_poly_consts[c](6);
			A(v+5,8) = dx_h*dx_h*dy_h - WENO_poly_consts[c](7);
			A(v+5,9) = dx_h*dy_h*dy_h - WENO_poly_consts[c](8);
		}
		
		// Row 10-13 (Corresponds to interior of the cell)     

		for (unsigned int v = 0; v < vertices_per_cell; v++) {
		
			vertex =  cell->vertex(v); 
			P(0) = 0.5*(vertex(0) + cell_center(0)); P(1) = 0.5*(vertex(1) + cell_center(1));
			
			dx_h = (P(0) - x0)/h; dy_h = (P(1) - y0)/h;
			
			// Fill the matrix 
			
			A(v+9,0) = 1.0; 
			A(v+9,1) = dx_h; 
			A(v+9,2) = dy_h;
			A(v+9,3) = dx_h*dx_h - WENO_poly_consts[c](2); 
			A(v+9,4) = dy_h*dy_h - WENO_poly_consts[c](3);
			A(v+9,5) = dx_h*dy_h - WENO_poly_consts[c](4);
			A(v+9,6) = dx_h*dx_h*dx_h - WENO_poly_consts[c](5); 
			A(v+9,7) = dy_h*dy_h*dy_h - WENO_poly_consts[c](6);
			A(v+9,8) = dx_h*dx_h*dy_h - WENO_poly_consts[c](7);
			A(v+9,9) = dx_h*dy_h*dy_h - WENO_poly_consts[c](8);
		}
		
		// Initialize the householder object 
		
		H_space[c].initialize(A);      
        
    } // Finished creating LHS matrices for spatial nodes 
    
    // =====================  Create the LHS Matrices for nodal to modal transcription in time =====================
    
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		cell = local_index_to_iterator[c];
        
		h = Cell[c].h(); 
		
		x0 = WENO_poly_consts[c](0); 
		y0 = WENO_poly_consts[c](1); 
		
		// tau = 1./3.
		
		// Asseble the matrix (Row 1 - cell center)
		A_time(0,0) = r1_3; 
		A_time(0,1) = (r1_3)*(r1_3); 
		A_time(0,2) = (r1_3)*(r1_3)*(r1_3);
		A_time(0,3) = 0.0;
		A_time(0,4) = 0.0;
		A_time(0,5) = 0.0;
		A_time(0,6) = 0.0;
		A_time(0,7) = (r1_3)*(- WENO_poly_consts[c](2));
		A_time(0,8) = (r1_3)*(- WENO_poly_consts[c](3));
		A_time(0,9) = (r1_3)*(- WENO_poly_consts[c](4));

		// Row 2-5 - face-centers

		for (unsigned int f = 0; f < faces_per_cell; f++) {

			face_center =  cell->face(f)->center(); 
			
			dx_h = (face_center(0) - x0)/h; dy_h = (face_center(1) - y0)/h;

			A_time(f+1,0) = r1_3; 
			A_time(f+1,1) = (r1_3)*(r1_3); 
			A_time(f+1,2) = (r1_3)*(r1_3)*(r1_3);
			A_time(f+1,3) = (r1_3)*dx_h;
			A_time(f+1,4) = (r1_3)*dy_h;
			A_time(f+1,5) = (r1_3)*(r1_3)*dx_h;
			A_time(f+1,6) = (r1_3)*(r1_3)*dy_h;
			A_time(f+1,7) = (r1_3)*(dx_h*dx_h - WENO_poly_consts[c](2));
			A_time(f+1,8) = (r1_3)*(dy_h*dy_h - WENO_poly_consts[c](3));
			A_time(f+1,9) = (r1_3)*(dx_h*dy_h - WENO_poly_consts[c](4));
		}
		
		for (unsigned int v = 0; v < vertices_per_cell; v++) {

			vertex = cell->vertex(v);
			
			dx_h = (vertex(0) - x0)/h; dy_h = (vertex(1) - y0)/h;
			
			A_time(v+5,0) = r1_3; 
			A_time(v+5,1) = (r1_3)*(r1_3); 
			A_time(v+5,2) = (r1_3)*(r1_3)*(r1_3);
			A_time(v+5,3) = (r1_3)*dx_h;
			A_time(v+5,4) = (r1_3)*dy_h;
			A_time(v+5,5) = (r1_3)*(r1_3)*dx_h;
			A_time(v+5,6) = (r1_3)*(r1_3)*dy_h;
			A_time(v+5,7) = (r1_3)*(dx_h*dx_h - WENO_poly_consts[c](2));
			A_time(v+5,8) = (r1_3)*(dy_h*dy_h - WENO_poly_consts[c](3));
			A_time(v+5,9) = (r1_3)*(dx_h*dy_h - WENO_poly_consts[c](4));
		}
		
		// tau = 2./3.
		
		// Asseble the matrix (Row 1 - cell center)
		A_time(9,0) = r2_3; 
		A_time(9,1) = (r2_3)*(r2_3); 
		A_time(9,2) = (r2_3)*(r2_3)*(r2_3);
		A_time(9,3) = 0.0;
		A_time(9,4) = 0.0;
		A_time(9,5) = 0.0;
		A_time(9,6) = 0.0;
		A_time(9,7) = (r2_3)*(- WENO_poly_consts[c](2));
		A_time(9,8) = (r2_3)*(- WENO_poly_consts[c](3));
		A_time(9,9) = (r2_3)*(- WENO_poly_consts[c](4));

		// Row 2-5 - face-centers

		for (unsigned int f = 0; f < faces_per_cell; f++) {

			face_center =  cell->face(f)->center(); 
			
			dx_h = (face_center(0) - x0)/h; dy_h = (face_center(1) - y0)/h;

			A_time(f+10,0) = r2_3; 
			A_time(f+10,1) = (r2_3)*(r2_3); 
			A_time(f+10,2) = (r2_3)*(r2_3)*(r2_3);
			A_time(f+10,3) = (r2_3)*dx_h;
			A_time(f+10,4) = (r2_3)*dy_h;
			A_time(f+10,5) = (r2_3)*(r2_3)*dx_h;
			A_time(f+10,6) = (r2_3)*(r2_3)*dy_h;
			A_time(f+10,7) = (r2_3)*(dx_h*dx_h - WENO_poly_consts[c](2));
			A_time(f+10,8) = (r2_3)*(dy_h*dy_h - WENO_poly_consts[c](3));
			A_time(f+10,9) = (r2_3)*(dx_h*dy_h - WENO_poly_consts[c](4));
		}
		
		// Row 15 - cell center, tau = 1.0
		A_time(14,0) = 1.0; 
		A_time(14,1) = 1.0; 
		A_time(14,2) = 1.0;
		A_time(14,3) = 0.0;
		A_time(14,4) = 0.0;
		A_time(14,5) = 0.0;
		A_time(14,6) = 0.0;
		A_time(14,7) = ( - WENO_poly_consts[c](2));
		A_time(14,8) = ( - WENO_poly_consts[c](3));
		A_time(14,9) = ( - WENO_poly_consts[c](4));
		
		H_time[c].initialize(A_time);      
	
    } // Finished creating LHS matrices for temporal nodes 
    
    // =====================  Create space-time Galerkin matrix =====================
    
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {
            
		A_space_time = Assemble_Space_Time_Galerkin_LHS_Matrix(Integration_Constants[c]);
		
		LU_space_time[c].initialize(A_space_time);      
        
    } // Finished creating LHS matrices for temporal nodes 

	Utilities::System::MemoryStats stat;
	Utilities::System::get_memory_stats(stat);
	pcout<<"solve Ader memory: "<<std::endl;
	pcout<<stat.VmRSS/std::pow(2,20)<<std::endl;
	pcout<<"total Memory consumption in Solve ader in GB: "<<24.0*stat.VmRSS/std::pow(2,20)<<std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "log.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Memory consumption in Solve ader per node in GB = " << stat.VmRSS/std::pow(2,20) << std::endl;
        fout_convergence.close();
	}
	 
	compute_time_step_based_on_cfl_number();

	unsigned int output_count;
	output_count = 0.1 / dt ;
	
	while (time < finalTime) {
		
		auto start_ader = std::chrono::system_clock::now();

        compute_time_step_based_on_cfl_number();

        if (count%output_count == 0 ) {

	        auto start = std::chrono::system_clock::now();

            output_results();

            auto end = std::chrono::system_clock::now();
            std::chrono::duration<double> elapsed_seconds = end-start;

            if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
 	            std::ofstream fout_convergence ;
 	            fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
                fout_convergence.precision(7) ;

                const std::string filename = "timer.dat";
                fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

                fout_convergence << "Time taken to write plot file = " << elapsed_seconds.count() << std::endl;
                fout_convergence.close();
            }

        }

		if (count%100 == 0 ) {
            restart();
        }

		if ((count+5)%100 == 0 ) {
            restart_r();
        }

		time += dt;
    	if ( (Utilities::MPI::this_mpi_process(mpi_communicator) == 0) && ( ( count%50 == 0 ) || std::fabs(time - finalTime) < 1e-8) ){
		    std::ofstream fout_convergence ;
    	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
    	    fout_convergence.precision(7) ;
	
    	    const std::string filename = "log.dat";
		    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    		fout_convergence <<"time = " << time <<",\t dt: "<<dt<< ",\t Final time = " << finalTime<< std::endl;
    	    fout_convergence.close();
		}
		pcout<<"time = " << time <<",\t dt: "<<dt<< ",\t Final time = " << finalTime<< std::endl;

        count++;

        // Get the coefficients 

		auto start_reconstruct = std::chrono::system_clock::now();

        reconstruct(); 

		auto end_reconstruct = std::chrono::system_clock::now();

		// =================================================  Predictor Step  =================================================
    	
		for (unsigned int c = 0; c < n_relevant_cells; ++c) {

			cell = local_index_to_iterator[c];
			
			// Get the list of nodes 		
			
			h = Cell[c].h();
							
			x0 = WENO_poly_consts[c](0); 
			y0 = WENO_poly_consts[c](1); 
			
			cell_center(0) = x0; 
			cell_center(1) = y0;

			node_list_space[0] = cell_center;  

			// Row 2-5 (Corresponds to face centers)     
			
			for (unsigned int f = 0; f < faces_per_cell; f++) {
				node_list_space[f+1] = cell->face(f)->center();
			}

			// Row 6-9 (Corresponds to vertices)     

			for (unsigned int v = 0; v < vertices_per_cell; v++) { 
				node_list_space[v+5] = cell->vertex(v); 
			}
			
			// Row 10-13 (Corresponds to interior points)
			
			for (unsigned int v = 0; v < vertices_per_cell; v++) {
				vertex =  cell->vertex(v);
				P(0) = 0.5*(vertex(0) + cell_center(0)); P(1) = 0.5*(vertex(1) + cell_center(1));
				node_list_space[v+9] = P; 
			}

			for (unsigned int i = 0; i < 13; i++) {
				
				// Evaluate Nodal values of Conserved Variables 
				U(0) = evaluate_weno_polynomial(coeffs_RHO[c],   WENO_poly_consts[c], node_list_space[i], h);
				U(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], node_list_space[i], h);
				U(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], node_list_space[i], h);
				U(3) = evaluate_weno_polynomial(coeffs_E[c],     WENO_poly_consts[c], node_list_space[i], h);

				// Convert to primitive variables 
				W = conserved_to_primitive(U);
		
				// Fill the vector
				F1_nodal_space(i) = (dt/h)*W(0)*W(1);                    G1_nodal_space(i) = (dt/h)*W(0)*W(2);                
				F2_nodal_space(i) = (dt/h)*(W(0)*W(1)*W(1) + W(3));      G2_nodal_space(i) = (dt/h)*W(0)*W(1)*W(2);     
				F3_nodal_space(i) = (dt/h)*W(0)*W(1)*W(2);               G3_nodal_space(i) = (dt/h)*(W(0)*W(2)*W(2) + W(3));        
				F4_nodal_space(i) = (dt/h)*W(1)*(U(3) + W(3));           G4_nodal_space(i) = (dt/h)*W(2)*(U(3) + W(3)); 
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

			for (unsigned int i = 0; i < 13; i++) {
					
					F1_nodal_space(i) = evaluate_weno_polynomial(F1_coeffs_space, WENO_poly_consts[c], node_list_space[i], h);  
					F2_nodal_space(i) = evaluate_weno_polynomial(F2_coeffs_space, WENO_poly_consts[c], node_list_space[i], h);
					F3_nodal_space(i) = evaluate_weno_polynomial(F3_coeffs_space, WENO_poly_consts[c], node_list_space[i], h);
					F4_nodal_space(i) = evaluate_weno_polynomial(F4_coeffs_space, WENO_poly_consts[c], node_list_space[i], h); 

					G1_nodal_space(i) = evaluate_weno_polynomial(G1_coeffs_space, WENO_poly_consts[c], node_list_space[i], h);  
					G2_nodal_space(i) = evaluate_weno_polynomial(G2_coeffs_space, WENO_poly_consts[c], node_list_space[i], h);
					G3_nodal_space(i) = evaluate_weno_polynomial(G3_coeffs_space, WENO_poly_consts[c], node_list_space[i], h);
					G4_nodal_space(i) = evaluate_weno_polynomial(G4_coeffs_space, WENO_poly_consts[c], node_list_space[i], h); 
							
			}  

			while (ITER < MAX_ITER) {
								
				// Get the time dependent nodal values of flux at tau = (1./3.)*dt

				for (unsigned int i = 0; i < 9; i++) {
					
					// Evaluate Nodal values of Conserved Variables 
					U(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], node_list_space[i], r1_3, h);
					U(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], node_list_space[i], r1_3, h);
					U(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], node_list_space[i], r1_3, h);
					U(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], node_list_space[i], r1_3, h);

					// Convert to primitive variables 
					W = conserved_to_primitive(U);

					// Evaluate corresponding nodal fluxes
					F1_nodal_time(i) = (dt/h)*W(0)*W(1);                    G1_nodal_time(i) = (dt/h)*W(0)*W(2);                
					F2_nodal_time(i) = (dt/h)*(W(0)*W(1)*W(1) + W(3));      G2_nodal_time(i) = (dt/h)*W(0)*W(1)*W(2);     
					F3_nodal_time(i) = (dt/h)*W(0)*W(1)*W(2);               G3_nodal_time(i) = (dt/h)*(W(0)*W(2)*W(2) + W(3));        
					F4_nodal_time(i) = (dt/h)*W(1)*(U(3) + W(3));           G4_nodal_time(i) = (dt/h)*W(2)*(U(3) + W(3)); 
				}
				
				// Get the time dependent nodal values of flux at tau = (2./3.)*dt

				for (unsigned int i = 0; i < 5; i++) {
					
					// Evaluate Nodal values of Conserved Variables 
					U(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], node_list_space[i], r2_3, h);
					U(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], node_list_space[i], r2_3, h);
					U(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], node_list_space[i], r2_3, h);
					U(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], node_list_space[i], r2_3, h);

					// Convert to primitive variables 
					W = conserved_to_primitive(U);

					// Evaluate corresponding nodal fluxes
					F1_nodal_time(i+9) = (dt/h)*W(0)*W(1);                    G1_nodal_time(i+9) = (dt/h)*W(0)*W(2);                
					F2_nodal_time(i+9) = (dt/h)*(W(0)*W(1)*W(1) + W(3));      G2_nodal_time(i+9) = (dt/h)*W(0)*W(1)*W(2);     
					F3_nodal_time(i+9) = (dt/h)*W(0)*W(1)*W(2);               G3_nodal_time(i+9) = (dt/h)*(W(0)*W(2)*W(2) + W(3));        
					F4_nodal_time(i+9) = (dt/h)*W(1)*(U(3) + W(3));           G4_nodal_time(i+9) = (dt/h)*W(2)*(U(3) + W(3)); 
				}

				// Get the time dependent nodal values of flux at tau = 1.0

				// Evaluate Nodal values of Conserved Variables 
				U(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], cell_center, 1.0, h);
				U(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], cell_center, 1.0, h);
				U(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], cell_center, 1.0, h);
				U(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], cell_center, 1.0, h);

				// Convert to primitive variables 
				W = conserved_to_primitive(U);

				// Evaluate corresponding nodal fluxes
				F1_nodal_time(14) = (dt/h)*W(0)*W(1);                    G1_nodal_time(14) = (dt/h)*W(0)*W(2);                
				F2_nodal_time(14) = (dt/h)*(W(0)*W(1)*W(1) + W(3));      G2_nodal_time(14) = (dt/h)*W(0)*W(1)*W(2);     
				F3_nodal_time(14) = (dt/h)*W(0)*W(1)*W(2);               G3_nodal_time(14) = (dt/h)*(W(0)*W(2)*W(2) + W(3));        
				F4_nodal_time(14) = (dt/h)*W(1)*(U(3) + W(3));           G4_nodal_time(14) = (dt/h)*W(2)*(U(3) + W(3)); 

				// Subtract the pure spatial part from the time nodal values to solve the least squares system 

				for (unsigned int i = 0; i < 9; i++) {
					
					F1_nodal_time(i) = F1_nodal_time(i) - F1_nodal_space(i); 
					F2_nodal_time(i) = F2_nodal_time(i) - F2_nodal_space(i);
					F3_nodal_time(i) = F3_nodal_time(i) - F3_nodal_space(i);
					F4_nodal_time(i) = F4_nodal_time(i) - F4_nodal_space(i);

					G1_nodal_time(i) = G1_nodal_time(i) - G1_nodal_space(i);  
					G2_nodal_time(i) = G2_nodal_time(i) - G2_nodal_space(i);
					G3_nodal_time(i) = G3_nodal_time(i) - G3_nodal_space(i);
					G4_nodal_time(i) = G4_nodal_time(i) - G4_nodal_space(i);
				}  
				
				
				for (unsigned int i = 0; i < 5; i++) {
					
					F1_nodal_time(i+9) = F1_nodal_time(i+9) - F1_nodal_space(i); 
					F2_nodal_time(i+9) = F2_nodal_time(i+9) - F2_nodal_space(i);
					F3_nodal_time(i+9) = F3_nodal_time(i+9) - F3_nodal_space(i);
					F4_nodal_time(i+9) = F4_nodal_time(i+9) - F4_nodal_space(i);

					G1_nodal_time(i+9) = G1_nodal_time(i+9) - G1_nodal_space(i);  
					G2_nodal_time(i+9) = G2_nodal_time(i+9) - G2_nodal_space(i);
					G3_nodal_time(i+9) = G3_nodal_time(i+9) - G3_nodal_space(i);
					G4_nodal_time(i+9) = G4_nodal_time(i+9) - G4_nodal_space(i);
				}
				
				F1_nodal_time(14) = F1_nodal_time(14) - F1_nodal_space(0); 
				F2_nodal_time(14) = F2_nodal_time(14) - F2_nodal_space(0);
				F3_nodal_time(14) = F3_nodal_time(14) - F3_nodal_space(0);
				F4_nodal_time(14) = F4_nodal_time(14) - F4_nodal_space(0);

				G1_nodal_time(14) = G1_nodal_time(14) - G1_nodal_space(0);  
				G2_nodal_time(14) = G2_nodal_time(14) - G2_nodal_space(0);
				G3_nodal_time(14) = G3_nodal_time(14) - G3_nodal_space(0);
				G4_nodal_time(14) = G4_nodal_time(14) - G4_nodal_space(0);

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

				Assemble_Space_Time_Galerkin_RHS_Vector(F1_coeffs_space, F1_coeffs_time, G1_coeffs_space, G1_coeffs_time, 
																		Integration_Constants[c], RHO_rhs_time); 
				
				Assemble_Space_Time_Galerkin_RHS_Vector(F2_coeffs_space, F2_coeffs_time, G2_coeffs_space, G2_coeffs_time, 
																		Integration_Constants[c], RHO_U_rhs_time);
				
				Assemble_Space_Time_Galerkin_RHS_Vector(F3_coeffs_space, F3_coeffs_time, G3_coeffs_space, G3_coeffs_time, 
																			Integration_Constants[c], RHO_V_rhs_time);
			
				Assemble_Space_Time_Galerkin_RHS_Vector(F4_coeffs_space, F4_coeffs_time, G4_coeffs_space, G4_coeffs_time, 
																		Integration_Constants[c], E_rhs_time); 
				
				// Update time-dependent coefficients 

				LU_space_time[c].solve(RHO_rhs_time, coeffs_RHO_time[c]);        
				LU_space_time[c].solve(RHO_U_rhs_time, coeffs_RHO_U_time[c]);
				LU_space_time[c].solve(RHO_V_rhs_time, coeffs_RHO_V_time[c]);
				LU_space_time[c].solve(E_rhs_time, coeffs_E_time[c]);

				ITER++; 
						
			} // End of iterations

			for (unsigned int f = 0; f < faces_per_cell; ++f) {

   	            nx1 = Cell[c].nx1(f); ny1 = Cell[c].ny1(f);
   	            nx2 = Cell[c].nx2(f); ny2 = Cell[c].ny2(f);
                
   	            face_quadrature_point_1 = Cell[c].face_quadrature_point1(f); 
                face_quadrature_point_2 = Cell[c].face_quadrature_point2(f); 
					
				h = std::sqrt( Cell[c].measure() ); 
				
				// Left face
				UL1(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0], h);
				UL1(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0], h);
				UL1(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0], h);
				UL1(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0], h);
			
				UL2(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0], h);
				UL2(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0], h);
				UL2(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0], h);
				UL2(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0], h);

				UL3(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1], h);
				UL3(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1], h);
				UL3(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1], h);
				UL3(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1], h);
					
				UL4(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1], h);
				UL4(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1], h);
				UL4(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1], h);
				UL4(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1], h);
				
				WL1 = conserved_to_primitive(UL1); WL2 = conserved_to_primitive(UL2);
				WL3 = conserved_to_primitive(UL3); WL4 = conserved_to_primitive(UL4);

	            if (WL1(0) < 0.0 || WL1(3) < 0.0 || WL2(0) < 0.0 || WL2(3) < 0.0 || WL3(0) < 0.0 || WL3(3) < 0.0 || WL4(0) < 0.0 || WL4(3) < 0.0) {

	                for (unsigned int i = 0; i < 10; i++) {
	                  coeffs_RHO_time[c](i) = 0.0; 
	                  coeffs_RHO_U_time[c](i) = 0.0;
	                  coeffs_RHO_V_time[c](i) = 0.0;
	                  coeffs_E_time[c](i) = 0.0;
	                }
	                for (unsigned int i = 1; i < 10; i++) {
	                  coeffs_RHO[c](i) = 0.0; 
	                  coeffs_RHO_U[c](i) = 0.0;
	                  coeffs_RHO_V[c](i) = 0.0;
	                  coeffs_E[c](i) = 0.0;
	                }
				}

			}  // check bound
		} // End of cell loop    
		
		// =================================================  Corrector Step  =================================================

		for (unsigned int f = 0; f < n_faces; f++) {
			did_not_compute_flux_for_the_face[f] = true; 
		}	

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

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
					
    	            nx1 = Cell[c].nx1(f); ny1 = Cell[c].ny1(f);
    	            nx2 = Cell[c].nx2(f); ny2 = Cell[c].ny2(f);
                
    	            face_quadrature_point_1 = Cell[c].face_quadrature_point1(f); 
	                face_quadrature_point_2 = Cell[c].face_quadrature_point2(f); 
					
					h = Cell[c].h(); 
				
					// Left face
					UL1(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0], h);
					UL1(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0], h);
					UL1(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0], h);
					UL1(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_1, t_gauss[0], h);
					
					UL2(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0], h);
					UL2(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0], h);
					UL2(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0], h);
					UL2(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_2, t_gauss[0], h);

					UL3(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1], h);
					UL3(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1], h);
					UL3(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1], h);
					UL3(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_1, t_gauss[1], h);
					
					UL4(0) = evaluate_ader_polynomial(coeffs_RHO[c],   coeffs_RHO_time[c],   WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1], h);
					UL4(1) = evaluate_ader_polynomial(coeffs_RHO_U[c], coeffs_RHO_U_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1], h);
					UL4(2) = evaluate_ader_polynomial(coeffs_RHO_V[c], coeffs_RHO_V_time[c], WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1], h);
					UL4(3) = evaluate_ader_polynomial(coeffs_E[c],     coeffs_E_time[c],     WENO_poly_consts[c], face_quadrature_point_2, t_gauss[1], h);
					
					WL1 = conserved_to_primitive(UL1); WL2 = conserved_to_primitive(UL2);
					WL3 = conserved_to_primitive(UL3); WL4 = conserved_to_primitive(UL4);
					
					// Get the right state values
			
					if (cell->face(f)->at_boundary()) {
                        
						WR1 = WL1; WR2 = WL2; WR3 = WL3; WR4 = WL4; 
						boundary = true;
    	            }
				
					else {
						
	    	        	DoFHandler<2>::cell_iterator neighbor = cell->neighbor(f);
			            neighbor->get_dof_indices(local_neighbor_dof_indices);
						unsigned int neighbor_c = global_to_local_index_map[local_neighbor_dof_indices[0]];
						
						h = Cell[neighbor_c].h(); 
						
						UR1(0) = evaluate_ader_polynomial(coeffs_RHO[neighbor_c], coeffs_RHO_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_1, t_gauss[0], h);
						
						UR1(1) = evaluate_ader_polynomial(coeffs_RHO_U[neighbor_c], coeffs_RHO_U_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_1, t_gauss[0], h);
						
						UR1(2) = evaluate_ader_polynomial(coeffs_RHO_V[neighbor_c], coeffs_RHO_V_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_1, t_gauss[0], h);
						

						UR1(3) = evaluate_ader_polynomial(coeffs_E[neighbor_c], coeffs_E_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_1, t_gauss[0], h);
						
						
						UR2(0) = evaluate_ader_polynomial(coeffs_RHO[neighbor_c], coeffs_RHO_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_2, t_gauss[0], h);
						
						UR2(1) = evaluate_ader_polynomial(coeffs_RHO_U[neighbor_c], coeffs_RHO_U_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_2, t_gauss[0], h);
						
						UR2(2) = evaluate_ader_polynomial(coeffs_RHO_V[neighbor_c], coeffs_RHO_V_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_2, t_gauss[0], h);
						
						UR2(3) = evaluate_ader_polynomial(coeffs_E[neighbor_c], coeffs_E_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_2, t_gauss[0], h); 

						UR3(0) = evaluate_ader_polynomial(coeffs_RHO[neighbor_c], coeffs_RHO_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_1, t_gauss[1], h);
						
						UR3(1) = evaluate_ader_polynomial(coeffs_RHO_U[neighbor_c], coeffs_RHO_U_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_1, t_gauss[1], h);
						
						UR3(2) = evaluate_ader_polynomial(coeffs_RHO_V[neighbor_c], coeffs_RHO_V_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_1, t_gauss[1], h);
						
						UR3(3) = evaluate_ader_polynomial(coeffs_E[neighbor_c], coeffs_E_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_1, t_gauss[1], h);
						
						UR4(0) = evaluate_ader_polynomial(coeffs_RHO[neighbor_c], coeffs_RHO_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_2, t_gauss[1], h);
						
						UR4(1) = evaluate_ader_polynomial(coeffs_RHO_U[neighbor_c], coeffs_RHO_U_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_2, t_gauss[1], h);
						
						UR4(2) = evaluate_ader_polynomial(coeffs_RHO_V[neighbor_c], coeffs_RHO_V_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_2, t_gauss[1], h);
						
						UR4(3) = evaluate_ader_polynomial(coeffs_E[neighbor_c], coeffs_E_time[neighbor_c],
														WENO_poly_consts[neighbor_c],  face_quadrature_point_2, t_gauss[1], h); 
						
						
						boundary = false;
						
						WR1 = conserved_to_primitive(UR1); WR2 = conserved_to_primitive(UR2);
						WR3 = conserved_to_primitive(UR3); WR4 = conserved_to_primitive(UR4);

					}
				
					Flux1[local_face_index] = rotated_HLLC_riemann_solver(WL1, WR1, nx1, ny1, face_quadrature_point_1, boundary);
					Flux2[local_face_index] = rotated_HLLC_riemann_solver(WL2, WR2, nx2, ny2, face_quadrature_point_2, boundary);
					Flux3[local_face_index] = rotated_HLLC_riemann_solver(WL3, WR3, nx1, ny1, face_quadrature_point_1, boundary);
					Flux4[local_face_index] = rotated_HLLC_riemann_solver(WL4, WR4, nx2, ny2, face_quadrature_point_2, boundary);

					did_not_compute_flux_for_the_face[local_face_index] = false; 
				}   
			
				else {
					
					Flux1[local_face_index] *= -1.0; 
					Flux2[local_face_index] *= -1.0;
					Flux3[local_face_index] *= -1.0;
					Flux4[local_face_index] *= -1.0;
				}
			
			
				F(0) = 0.25*(Flux1[local_face_index](0) + Flux2[local_face_index](0) + Flux3[local_face_index](0) + Flux4[local_face_index](0)); 
				F(1) = 0.25*(Flux1[local_face_index](1) + Flux2[local_face_index](1) + Flux3[local_face_index](1) + Flux4[local_face_index](1)); 
				F(2) = 0.25*(Flux1[local_face_index](2) + Flux2[local_face_index](2) + Flux3[local_face_index](2) + Flux4[local_face_index](2)); 
				F(3) = 0.25*(Flux1[local_face_index](3) + Flux2[local_face_index](3) + Flux3[local_face_index](3) + Flux4[local_face_index](3)); 

				// Add it to the rhs vectors

				rhs1(c) += (-1.0/V_c)*(F(0)*S_f);
				rhs2(c) += (-1.0/V_c)*(F(1)*S_f);
				rhs3(c) += (-1.0/V_c)*(F(2)*S_f);
				rhs4(c) += (-1.0/V_c)*(F(3)*S_f);
			}
        }  // End of cell loop   

        // Update the variables  

		for (unsigned int c = 0; c < n_locally_cells; ++c) {

			g_i = local_to_global_index_map[c];

			local_RHO(g_i) = local_RHO(g_i) + dt*rhs1(c);
			local_RHO_U(g_i) = local_RHO_U(g_i) + dt*rhs2(c);
			local_RHO_V(g_i) = local_RHO_V(g_i) + dt*rhs3(c);
			local_E(g_i) = local_E(g_i) + dt*rhs4(c);	
		}

	    auto start_com = std::chrono::system_clock::now();

		local_RHO.compress(VectorOperation::insert);
		local_RHO_U.compress(VectorOperation::insert);
		local_RHO_V.compress(VectorOperation::insert);
		local_E.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E; 

	    auto end_com = std::chrono::system_clock::now();

	    std::chrono::duration<double> elapsed_seconds_recon = end_reconstruct - start_reconstruct;
	    std::chrono::duration<double> elapsed_seconds_com = end_com - start_com;
	    std::chrono::duration<double> elapsed_seconds_ader = end_com - start_ader;

    	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0 && count < 10){
		    std::ofstream fout_convergence ;
    	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
    	    fout_convergence.precision(7) ;
	
    	    const std::string filename = "timer.dat";
		    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    		fout_convergence << "time taken by reconstruction = " << elapsed_seconds_recon.count() << std::endl;	
    		fout_convergence << "time taken by data transfer of ader = " << elapsed_seconds_com.count() << std::endl;
    		fout_convergence << "time taken by 1 step of ader = " << elapsed_seconds_ader.count() << std::endl;
    	    fout_convergence.close();
		}		
		
    } // End of time loop 

    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "timer.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "time taken by ader = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}
    
	output_results();
	restart();
 
} // End of function 

double evaluate_ader_polynomial(Vector<double> coeffs, Vector<double> coeffs_time, Vector<double> consts,  Point<2> P, double t, double h) {

    // Inputs are the coefficients of the polynomial, the center of the cell and the point of evaluation
	
	double dx = (P(0)-consts(0))/h; double dy = (P(1)-consts(1))/h; 
    
    // Use Honers Algorith 
    
    double space_part = coeffs(0) -  coeffs(3)*consts(2) - coeffs(4)*consts(3) - coeffs(5)*consts(4) - 
			            coeffs(6)*consts(5) - coeffs(7)*consts(6) - coeffs(8)*consts(7) - coeffs(9)*consts(8) +
                        dx*(coeffs(1) + dy*coeffs(5) + dx*(coeffs(3) + dy*coeffs(8) + dx*coeffs(6)) ) + 
                        dy*(coeffs(2) + dy*(coeffs(4) + dy*coeffs(7) + dx*coeffs(9) ) ); 
    
    
    double time_part = t*(coeffs_time(0) + 
	                      coeffs_time(3)*dx + 
						  coeffs_time(4)*dy + 
						  coeffs_time(7)*(dx*dx-consts(2)) +
						  coeffs_time(8)*(dy*dy-consts(3)) +
						  coeffs_time(9)*(dx*dy-consts(4)) +
						  t*(coeffs_time(1) + t*coeffs_time(2) + coeffs_time(5)*dx + coeffs_time(6)*dy)); 
    
    return (space_part + time_part); 

}

FullMatrix<double> Assemble_Space_Time_Galerkin_LHS_Matrix(Vector<double> Int) {
    
    FullMatrix<double> A(10,10);
    
    double r2_3 = 2./3.; double r1_3 = 1./3.; 
    
    // Phi1 = t 
    A(0,0) = 0.5*Int[0];     
    A(0,1) = r2_3*Int[0];   
    A(0,2) = 0.75*Int[0];
    A(0,3) = 0.0; 
    A(0,4) = 0.0; 
    A(0,5) = 0.0; 
    A(0,6) = 0.0;
    A(0,7) = 0.0; 
    A(0,8) = 0.0;
    A(0,9) = 0.0;
    
    // Phi2 = t^2 
    A(1,0) = r1_3*Int[0];     
    A(1,1) = 0.5*Int[0];   
    A(1,2) = 0.6*Int[0];
    A(1,3) = 0.0; 
    A(1,4) = 0.0; 
    A(1,5) = 0.0; 
    A(1,6) = 0.0;
    A(1,7) = 0.0; 
    A(1,8) = 0.0;
    A(1,9) = 0.0;
    
    // Phi3 = t^3 
    A(2,0) = 0.25*Int[0];     
    A(2,1) = 0.4*Int[0];   
    A(2,2) = 0.5*Int[0];
    A(2,3) = 0.0; 
    A(2,4) = 0.0; 
    A(2,5) = 0.0; 
    A(2,6) = 0.0;
    A(2,7) = 0.0; 
    A(2,8) = 0.0;
    A(2,9) = 0.0;
    
    // Phi4 = (x-c1)*t
    A(3,0) = 0.0;     
    A(3,1) = 0.0;   
    A(3,2) = 0.0;
    A(3,3) = 0.5*Int[5]; 
    A(3,4) = 0.5*Int[6]; 
    A(3,5) = r2_3*Int[5]; 
    A(3,6) = r2_3*Int[6];
    A(3,7) = 0.5*Int[15]; 
    A(3,8) = 0.5*Int[16];
    A(3,9) = 0.5*Int[17];
    
    // Phi5 = (y-c2)*t
    A(4,0) = 0.0;     
    A(4,1) = 0.0;   
    A(4,2) = 0.0;
    A(4,3) = 0.5*Int[6]; 
    A(4,4) = 0.5*Int[9]; 
    A(4,5) = r2_3*Int[6]; 
    A(4,6) = r2_3*Int[9];
    A(4,7) = 0.5*Int[18]; 
    A(4,8) = 0.5*Int[19];
    A(4,9) = 0.5*Int[20];
    
    // Phi6 = (x-c1)*t^2
    A(5,0) = 0.0;     
    A(5,1) = 0.0;   
    A(5,2) = 0.0;
    A(5,3) = r1_3*Int[5]; 
    A(5,4) = r1_3*Int[6]; 
    A(5,5) = 0.5*Int[5]; 
    A(5,6) = 0.5*Int[6];
    A(5,7) = r1_3*Int[15]; 
    A(5,8) = r1_3*Int[16];
    A(5,9) = r1_3*Int[17];
    
    // Phi7 = (y-c2)*t^2
    A(6,0) = 0.0;     
    A(6,1) = 0.0;   
    A(6,2) = 0.0;
    A(6,3) = r1_3*Int[6]; 
    A(6,4) = r1_3*Int[9]; 
    A(6,5) = 0.5*Int[6]; 
    A(6,6) = 0.5*Int[9];
    A(6,7) = r1_3*Int[18]; 
    A(6,8) = r1_3*Int[19];
    A(6,9) = r1_3*Int[20];
    
    // Phi8 = (x^2-c3)*t 
    A(7,0) = 0.0;     
    A(7,1) = 0.0;   
    A(7,2) = 0.0;
    A(7,3) = 0.5*Int[15]; 
    A(7,4) = 0.5*Int[18]; 
    A(7,5) = r2_3*Int[15]; 
    A(7,6) = r2_3*Int[18];
    A(7,7) = 0.5*Int[21]; 
    A(7,8) = 0.5*Int[22];
    A(7,9) = 0.5*Int[23];
    
    // Phi9 = (y^2-c4)*t 
    A(8,0) = 0.0;     
    A(8,1) = 0.0;   
    A(8,2) = 0.0;
    A(8,3) = 0.5*Int[16]; 
    A(8,4) = 0.5*Int[19]; 
    A(8,5) = r2_3*Int[16]; 
    A(8,6) = r2_3*Int[19];
    A(8,7) = 0.5*Int[22]; 
    A(8,8) = 0.5*Int[24];
    A(8,9) = 0.5*Int[25];
    
    // Phi10 = (xy-c4)*t
    A(9,0) = 0.0;     
    A(9,1) = 0.0;   
    A(9,2) = 0.0;
    A(9,3) = 0.5*Int[17]; 
    A(9,4) = 0.5*Int[20]; 
    A(9,5) = r2_3*Int[17]; 
    A(9,6) = r2_3*Int[20];
    A(9,7) = 0.5*Int[23]; 
    A(9,8) = 0.5*Int[25];
    A(9,9) = 0.5*Int[26];

    return A; 

}

void Assemble_Space_Time_Galerkin_RHS_Vector(Vector<double> FS, Vector<double> FT, 
                                                       Vector<double> GS, Vector<double> GT, 
                                                       Vector<double> Int, Vector<double>& b) {
    double r1_3 = 1./3.;   double r2_3 = 2./3.;  
    
    double F_x = FS[1]; double F_xx = FS[3]; double F_xy = FS[5]; double F_xxx = FS[6]; double F_xxy = FS[8]; double F_xyy = FS[9]; 
    double F_xt = FT[3]; double F_xtt = FT[5]; double F_xxt = FT[7]; double F_xyt = FT[9]; 
    
    double G_y = GS[2]; double G_yy = GS[4]; double G_xy = GS[5]; double G_yyy = GS[7]; double G_xxy = GS[8]; double G_xyy = GS[9]; 
    double G_yt = GT[4]; double G_ytt = GT[6]; double G_yyt = GT[8]; double G_xyt = GT[9]; 
    
    // Phi1 = t 
    
    b(0) = -(0.5*F_x*Int[0] + F_xx*Int[1] + 0.5*F_xy*Int[2] + 1.5*F_xxx*Int[27] + F_xxy*Int[28] + 0.5*F_xyy*Int[29] +  
             r1_3*F_xt*Int[0] + 0.25*F_xtt*Int[0] + r2_3*F_xxt*Int[1] + r1_3*F_xyt*Int[2])
           
           -(0.5*G_y*Int[0] + G_yy*Int[2] + 0.5*G_xy*Int[1] + 1.5*G_yyy*Int[29] + 0.5*G_xxy*Int[27] + G_xyy*Int[28] +  
             r1_3*G_yt*Int[0] + 0.25*G_ytt*Int[0] + r2_3*G_yyt*Int[2] + r1_3*G_xyt*Int[1]); 
          
    // Phi2 = t^2
           
    b(1) = -(r1_3*F_x*Int[0] + r2_3*F_xx*Int[1] + r1_3*F_xy*Int[2] + F_xxx*Int[27] + r2_3*F_xxy*Int[28] + r1_3*F_xyy*Int[29] +  
             0.25*F_xt*Int[0] + 0.2*F_xtt*Int[0] + 0.5*F_xxt*Int[1] + 0.25*F_xyt*Int[2]) 
           
           -(r1_3*G_y*Int[0] + r2_3*G_yy*Int[2] + r1_3*G_xy*Int[1] + G_yyy*Int[29] + r1_3*G_xxy*Int[27] + r2_3*G_xyy*Int[28] +  
             0.25*G_yt*Int[0] + 0.2*G_ytt*Int[0] + 0.5*G_yyt*Int[2] + 0.25*G_xyt*Int[1]); 
    
    
    // Phi3 = t^3
           
    b(2) = -(0.25*F_x*Int[0] + 0.5*F_xx*Int[1] + 0.25*F_xy*Int[2] + 0.75*F_xxx*Int[27] + 0.5*F_xxy*Int[28]+0.25*F_xyy*Int[29] 
           + 0.2*F_xt*Int[0] + 0.5*r1_3*F_xtt*Int[0] + 0.4*F_xxt*Int[1] + 0.2*F_xyt*Int[2]) 
           
           -(0.25*G_y*Int[0] + 0.5*G_yy*Int[2] + 0.25*G_xy*Int[1] + 0.75*G_yyy*Int[29] + 0.25*G_xxy*Int[27]+0.5*G_xyy*Int[28]
           + 0.2*G_yt*Int[0] + 0.5*r1_3*G_ytt*Int[0] + 0.4*G_yyt*Int[2] + 0.2*G_xyt*Int[1]);
         
    // Phi4 = (x-c1)*t

    b(3) = -(0.5*F_x*Int[3] + F_xx*Int[7] + 0.5*F_xy*Int[8] + 1.5*F_xxx*Int[30] + F_xxy*Int[31] + 0.5*F_xyy*Int[32] +  
             r1_3*F_xt*Int[3] + 0.25*F_xtt*Int[3] + r2_3*F_xxt*Int[7] + r1_3*F_xyt*Int[8])
           
            -(0.5*G_y*Int[3] + G_yy*Int[8] + 0.5*G_xy*Int[7] + 1.5*G_yyy*Int[32] + 0.5*G_xxy*Int[30] + G_xyy*Int[31] +  
             r1_3*G_yt*Int[3] + 0.25*G_ytt*Int[3] + r2_3*G_yyt*Int[8] + r1_3*G_xyt*Int[7]); 
    
    
    // Phi5 = (y-c2)*t
           
    b(4) = -(0.5*F_x*Int[4] + F_xx*Int[10] + 0.5*F_xy*Int[11] + 1.5*F_xxx*Int[33] + F_xxy*Int[34] + 0.5*F_xyy*Int[35] +  
             r1_3*F_xt*Int[4] + 0.25*F_xtt*Int[4] + r2_3*F_xxt*Int[10] + r1_3*F_xyt*Int[11])
           
            -(0.5*G_y*Int[4] + G_yy*Int[11] + 0.5*G_xy*Int[10] + 1.5*G_yyy*Int[35] + 0.5*G_xxy*Int[33] + G_xyy*Int[34] +  
             r1_3*G_yt*Int[4] + 0.25*G_ytt*Int[4] + r2_3*G_yyt*Int[11] + r1_3*G_xyt*Int[10]);
    
    // Phi6 = (x-c1)*t^2 
           
    b(5) = -(r1_3*F_x*Int[3] + r2_3*F_xx*Int[7] + r1_3*F_xy*Int[8] + F_xxx*Int[30] + r2_3*F_xxy*Int[31] + r1_3*F_xyy*Int[32] +  
             0.25*F_xt*Int[3] + 0.2*F_xtt*Int[3] + 0.5*F_xxt*Int[7] + 0.25*F_xyt*Int[8]) 
           
           -(r1_3*G_y*Int[3] + r2_3*G_yy*Int[8] + r1_3*G_xy*Int[7] + G_yyy*Int[35] + r1_3*G_xxy*Int[33] + r2_3*G_xyy*Int[34] +  
             0.25*G_yt*Int[3] + 0.2*G_ytt*Int[3] + 0.5*G_yyt*Int[8] + 0.25*G_xyt*Int[7]);
    
    // Phi7 = (y-c2)*t^2
    
    b(6) = -(r1_3*F_x*Int[4] + r2_3*F_xx*Int[10] + r1_3*F_xy*Int[11] + F_xxx*Int[33] + r2_3*F_xxy*Int[34] + r1_3*F_xyy*Int[35] +  
             0.25*F_xt*Int[4] + 0.2*F_xtt*Int[4] + 0.5*F_xxt*Int[10] + 0.25*F_xyt*Int[11]) 
           
           -(r1_3*G_y*Int[4] + r2_3*G_yy*Int[11] + r1_3*G_xy*Int[10] + G_yyy*Int[35] + r1_3*G_xxy*Int[33] + r2_3*G_xyy*Int[34] +  
             0.25*G_yt*Int[4] + 0.2*G_ytt*Int[4] + 0.5*G_yyt*Int[11] + 0.25*G_xyt*Int[10]);
    
    // Phi8 = (x^2-c3)*t
    b(7) = -(0.5*F_x*Int[12] + F_xx*Int[36] + 0.5*F_xy*Int[37] + 1.5*F_xxx*Int[38] + F_xxy*Int[39] + 0.5*F_xyy*Int[40] +  
             r1_3*F_xt*Int[12] + 0.25*F_xtt*Int[12] + r2_3*F_xxt*Int[36] + r1_3*F_xyt*Int[37])
           
           -(0.5*G_y*Int[12] + G_yy*Int[37] + 0.5*G_xy*Int[36] + 1.5*G_yyy*Int[40] + 0.5*G_xxy*Int[38] + G_xyy*Int[39] +  
             r1_3*G_yt*Int[12] + 0.25*G_ytt*Int[12] + r2_3*G_yyt*Int[37] + r1_3*G_xyt*Int[36]); 
    
    
    // Phi9 = (y^2-c4)*t
           
    
    b(8) = -(0.5*F_x*Int[13] + F_xx*Int[41] + 0.5*F_xy*Int[42] + 1.5*F_xxx*Int[43] + F_xxy*Int[44] + 0.5*F_xyy*Int[45] +  
             r1_3*F_xt*Int[13] + 0.25*F_xtt*Int[13] + r2_3*F_xxt*Int[41] + r1_3*F_xyt*Int[42])
           
           -(0.5*G_y*Int[13] + G_yy*Int[42] + 0.5*G_xy*Int[41] + 1.5*G_yyy*Int[45] + 0.5*G_xxy*Int[43] + G_xyy*Int[44] +  
             r1_3*G_yt*Int[13] + 0.25*G_ytt*Int[13] + r2_3*G_yyt*Int[42] + r1_3*G_xyt*Int[41]);
    
    // Phi10 = (xy-c5)*t
    
    
    b(9) = -(0.5*F_x*Int[14] + F_xx*Int[46] + 0.5*F_xy*Int[47] + 1.5*F_xxx*Int[48] + F_xxy*Int[49] + 0.5*F_xyy*Int[50] +  
             r1_3*F_xt*Int[14] + 0.25*F_xtt*Int[13] + r2_3*F_xxt*Int[46] + r1_3*F_xyt*Int[47])
           
           -(0.5*G_y*Int[14] + G_yy*Int[47] + 0.5*G_xy*Int[46] + 1.5*G_yyy*Int[50] + 0.5*G_xxy*Int[48] + G_xyy*Int[49] +  
             r1_3*G_yt*Int[14] + 0.25*G_ytt*Int[13] + r2_3*G_yyt*Int[47] + r1_3*G_xyt*Int[46]);
}


#include "../include/Weno32.h" 

void Weno3_2D::compute_cell_properties() {
	
	std::cout << "Computing properties for each cell" << "\n";

    // Iterate over all the cells 
    
    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    
    // For obtaining face normal vectors
	QGauss<2-1> face_quadrature_formula(2);
	FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);

    Point<2> face_quadrature_point_1;     
    Point<2> face_quadrature_point_2;  
    
    double S_f; // Volume of the cell and surface area of the face
	
	double nx1, ny1;   // Face normal vectors
	double nx2, ny2; 
	
	Tensor<1,2> face_normal_vector1; // Face normal vector
	Tensor<1,2> face_normal_vector2; // Face normal vector
    
    double area; 
    Point<2> quadrature_points1[4];
    Point<2> quadrature_points2[4]; 
    double face_normal_x1[4]; 
    double face_normal_x2[4]; 
    double face_normal_y1[4]; 
    double face_normal_y2[4];
    double surface_area[4];
    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
		
		area = cell->measure(); 
		
		for (unsigned int f = 0; f < 4; f++) {
			
			fv_face_values.reinit(cell, f);
			
			face_normal_vector1 = fv_face_values.normal_vector(0); // Zero corresponds to the presence of only one quadrature point
			nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
			
			face_normal_vector2 = fv_face_values.normal_vector(1);
			nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
			
			face_quadrature_point_1 = fv_face_values.quadrature_point(0); 
			face_quadrature_point_2 = fv_face_values.quadrature_point(1); 
			
			S_f = cell->face(f)->measure(); 
			
			// Fill the values 
			
			quadrature_points1[f] = face_quadrature_point_1;  
			quadrature_points2[f] = face_quadrature_point_2;  
			face_normal_x1[f] = nx1; ; 
			face_normal_x2[f] = nx2; 
			face_normal_y1[f] = ny1; 
			face_normal_y2[f] = ny2, 
			surface_area[f] = S_f; 
		}
		
		Cell[c].reinit(area, quadrature_points1, quadrature_points2, face_normal_x1, face_normal_x2, face_normal_y1, face_normal_y2, surface_area);
		
    }
    
	std::cout << "Done!" << std::endl;
    std::cout << "===========================" << std::endl;
}

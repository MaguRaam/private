#include "../include/Euler.h"

void Euler_2D::make_grid() {

    std::cout << "============================" << std::endl;
    std::cout << "Making grid" << std::endl;
    GridIn<2> grid_in;
	grid_in.attach_triangulation (triangulation);
	std::ifstream input_file("mesh.msh");
	grid_in.read_msh(input_file);

	std::cout << "Number of active cells: "
						<< triangulation.n_active_cells()
                        << std::endl;
    
    // Assign boundary ids 
    
    // For obtaining face normal vectors
	QGauss<2-1> face_quadrature_formula(1);
	FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_normal_vectors);
    
    double nx, ny;
    Tensor<1,2> face_normal_vector; 
    
    Triangulation<2>::active_cell_iterator cell = triangulation.begin_active();
	Triangulation<2>::active_cell_iterator endc = triangulation.end();
    
    Point<2> face_center;
    
    for (; cell!=endc; ++cell) {
		for (unsigned int f=0; f < GeometryInfo<2>::faces_per_cell; ++f) {
			if (cell->face(f)->at_boundary()) {
				fv_face_values.reinit(cell, f);
				face_normal_vector = fv_face_values.normal_vector(0);
				nx = face_normal_vector[0]; ny = face_normal_vector[1];
				
				if (nx == -1.0) {
					cell->face(f)->set_boundary_id(0); 
				}
				
				if (nx == 1.0) {
					cell->face(f)->set_boundary_id(1); 
				}
				
				else {
					cell->face(f)->set_boundary_id(2);
				}
			}
		}
	}
    
    std::cout << "Done!" << std::endl;
    std::cout << "============================" << std::endl;
    
    
} 

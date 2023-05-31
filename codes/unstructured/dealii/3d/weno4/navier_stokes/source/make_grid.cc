#include "../include/Weno432.h"

void Weno4_3D::make_grid() {

	double x_min = 0.0, x_max = 1.0;
	double y_min = 0.0, y_max = 1.0;
	double z_min = 0.0, z_max = 1.0;

	unsigned int n_cell_x = 15, n_cell_y = 15, n_cell_z = 15;

    std::vector<unsigned int> repetions(3); // No. of n_cells in x and y directions 
    repetions[0] = n_cell_x;
    repetions[1] = n_cell_y; 
    repetions[2] = n_cell_z; 
    
    bool colorize = false;  // Set boundary ids for the four boundaries 

    Point<3> P1(x_min, y_min, z_min);
    Point<3> P2(x_max, y_max, z_max);
    GridGenerator::subdivided_hyper_rectangle(triangulation, repetions, P1, P2, colorize);

//	GridTools::transform(&grid_transform, triangulation);

    Triangulation<3>::active_cell_iterator cell = triangulation.begin_active();
	Triangulation<3>::active_cell_iterator endc = triangulation.end();
	Point<3> face_center;
	unsigned int reflective = 0;
    for (; cell!=endc; ++cell)
	    for (unsigned int f=0; f < GeometryInfo<3>::faces_per_cell; ++f)
	        if (cell->face(f)->at_boundary())
	            cell->face(f)->set_boundary_id(2);  // zero dirichlet
            
	pcout<<"no of reflective boundary: "<<reflective<<std::endl;
/*
	{
	    std::ofstream out ("grid_2.eps");
	    GridOut grid_out;
	    grid_out.write_eps (triangulation, out);
	}

	{
	    std::ofstream out ("grid_1.vtk");
	    GridOut grid_out;
	    grid_out.write_vtk (triangulation, out);
	}
*/
}

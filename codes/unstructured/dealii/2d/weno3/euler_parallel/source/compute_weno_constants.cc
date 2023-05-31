#include "../include/Weno32.h"


void Weno3_2D::compute_weno_polynomial_constants() {
    
    pcout << "Computing WENO polynmial constants" << std::endl;

    unsigned int N_gp = 2;               // No. of quadrature points

	PETScWrappers::MPI::Vector locally_WENO_poly_consts_x(locally_owned_dofs,  MPI_COMM_WORLD);
	PETScWrappers::MPI::Vector locally_WENO_poly_consts_y(locally_owned_dofs,  MPI_COMM_WORLD);
	PETScWrappers::MPI::Vector locally_WENO_poly_consts_xx(locally_owned_dofs,  MPI_COMM_WORLD);
	PETScWrappers::MPI::Vector locally_WENO_poly_consts_yy(locally_owned_dofs,  MPI_COMM_WORLD);
	PETScWrappers::MPI::Vector locally_WENO_poly_consts_xy(locally_owned_dofs,  MPI_COMM_WORLD);

    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    Point<2> q_point;

    double V0, x_0, y_0, j_w;

	unsigned int c, g_i;

    for (; cell != endc; ++cell) {

	if (cell->is_locally_owned()){

        cell->get_dof_indices(local_dof_indices);
		c = global_to_local_index_map[local_dof_indices[0] ]; 
		g_i = local_dof_indices[0];
        V0 = Cell[c].measure();    

		double const_x = 0.0, const_y = 0.0;
		double const_xx = 0.0, const_yy = 0.0, const_xy = 0.0;

        for (unsigned int i = 0; i < N_gp*N_gp; i++) {
            q_point = Cell[c].cell_quadrature_point(i);
			j_w = Cell[c].jxw(i);
            const_x += (1./V0)*j_w*(q_point(0));
            const_y += (1./V0)*j_w*(q_point(1));
        }
        
        x_0 = const_x; 
        y_0 = const_y;

        for (unsigned int i = 0; i < N_gp*N_gp; i++) {

            q_point = Cell[c].cell_quadrature_point(i);
			j_w = Cell[c].jxw(i);
            const_xx += (1./V0)*j_w*(q_point(0)-x_0)*(q_point(0)-x_0);
            const_yy += (1./V0)*j_w*(q_point(1)-y_0)*(q_point(1)-y_0);
            const_xy += (1./V0)*j_w*(q_point(0)-x_0)*(q_point(1)-y_0);
        }

		locally_WENO_poly_consts_x(g_i) = const_x;
		locally_WENO_poly_consts_y(g_i) = const_y;
		locally_WENO_poly_consts_xx(g_i) = const_xx;
		locally_WENO_poly_consts_yy(g_i) = const_yy;
		locally_WENO_poly_consts_xy(g_i) = const_xy;
        
    }

	}
	locally_WENO_poly_consts_x.compress(VectorOperation::insert);
	locally_WENO_poly_consts_y.compress(VectorOperation::insert);
	locally_WENO_poly_consts_xx.compress(VectorOperation::insert);
	locally_WENO_poly_consts_yy.compress(VectorOperation::insert);
	locally_WENO_poly_consts_xy.compress(VectorOperation::insert);

	WENO_poly_consts_x = locally_WENO_poly_consts_x;
	WENO_poly_consts_y = locally_WENO_poly_consts_y;
	WENO_poly_consts_xx = locally_WENO_poly_consts_xx;
	WENO_poly_consts_yy = locally_WENO_poly_consts_yy;
	WENO_poly_consts_xy = locally_WENO_poly_consts_xy;
    
    pcout << "Done!" << std::endl;

}

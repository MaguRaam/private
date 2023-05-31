#include "../include/Weno432.h"


void Weno4_3D::compute_weno_polynomial_constants() {

    pcout << "Computing WENO polynmial constants" << std::endl;

    unsigned int N_gp = 2;               // No. of quadrature points
	unsigned int n_quad_points = N_gp*N_gp*N_gp;
    Point<3> q_point;

    double V0, x0, y0, z0, j_w, dx, dy, dz;
	
	double h;

	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

        V0 = Cell[c].measure();

		h = Cell[c].h();

		for (unsigned int i = 0; i < 19; i++) {
			WENO_poly_consts[c](i) = 0.0;
		}

		for (unsigned int i = 0; i < n_quad_points; i++) {
            q_point = Cell[c].cell_quadrature_point(i);
			j_w = Cell[c].jxw(i);

			WENO_poly_consts[c](0) += (1./V0)*j_w*(q_point(0));
			WENO_poly_consts[c](1) += (1./V0)*j_w*(q_point(1));
			WENO_poly_consts[c](2) += (1./V0)*j_w*(q_point(2));
		}
		
		x0 = WENO_poly_consts[c](0); 
		y0 = WENO_poly_consts[c](1);
		z0 = WENO_poly_consts[c](2);

		for (unsigned int i = 0; i < n_quad_points; i++) {

            q_point = Cell[c].cell_quadrature_point(i);
			j_w = Cell[c].jxw(i);

			dx = (q_point(0) - x0)/h;
			dy = (q_point(1) - y0)/h;
			dz = (q_point(2) - z0)/h;

			WENO_poly_consts[c](3)  += (1./V0)*j_w*dx*dx;
			WENO_poly_consts[c](4)  += (1./V0)*j_w*dy*dy;
			WENO_poly_consts[c](5)  += (1./V0)*j_w*dz*dz;
			WENO_poly_consts[c](6)  += (1./V0)*j_w*dx*dy;
			WENO_poly_consts[c](7)  += (1./V0)*j_w*dx*dz;
			WENO_poly_consts[c](8)  += (1./V0)*j_w*dy*dz;
			WENO_poly_consts[c](9)  += (1./V0)*j_w*dx*dx*dx;
			WENO_poly_consts[c](10) += (1./V0)*j_w*dy*dy*dy;
			WENO_poly_consts[c](11) += (1./V0)*j_w*dz*dz*dz;
			WENO_poly_consts[c](12) += (1./V0)*j_w*dx*dx*dy;
			WENO_poly_consts[c](13) += (1./V0)*j_w*dx*dx*dz;
			WENO_poly_consts[c](14) += (1./V0)*j_w*dy*dy*dx;
			WENO_poly_consts[c](15) += (1./V0)*j_w*dy*dy*dz;
			WENO_poly_consts[c](16) += (1./V0)*j_w*dz*dz*dx;
			WENO_poly_consts[c](17) += (1./V0)*j_w*dz*dz*dy;
			WENO_poly_consts[c](18) += (1./V0)*j_w*dx*dy*dz;
		}
    }

    pcout << "Done!" << std::endl;

}


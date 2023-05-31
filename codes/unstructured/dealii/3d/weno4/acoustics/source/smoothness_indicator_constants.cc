#include "../include/Weno432.h" 


// Compute constants for smoothness indicators 

void Weno4_3D::compute_IS_constants() {
    
    pcout << "Computing smoothness indicator constants" << std::endl;
    
    unsigned int N_gp = 4;               // No. of quadrature points
    QGauss<3> quadrature_formula(N_gp);

	FEValues<3> fv_values (mapping, fv, quadrature_formula, update_quadrature_points | update_JxW_values);

//    FEValues<3> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);	

    DoFHandler<3>::active_cell_iterator cell;

    Point<3> q_point;
	
	double x0, y0, z0; 
	double x, y, z, dx, dy, dz; 
	double h, h3;
    
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		cell = local_index_to_iterator[c];

		IS_constants[c] = 0.0; 

		fv_values.reinit(cell);
		h3 = 0.0; 
		for (unsigned int i=0; i<fv_values.n_quadrature_points; ++i)
            h3 += fv_values.JxW (i);	

		IS_constants[c](0) = h3;
		
		x0 = WENO_poly_consts[c](0); 
		y0 = WENO_poly_consts[c](1);
		z0 = WENO_poly_consts[c](2);
		
		h = Cell[c].h();

		for (unsigned int i = 0; i < N_gp*N_gp*N_gp; i++) {

			q_point = fv_values.quadrature_point(i); 
			
			x = q_point(0); y = q_point(1); z = q_point(2);
			dx = (x-x0)/h;
			dy = (y-y0)/h;
			dz = (z-z0)/h;

			IS_constants[c](1)  += fv_values.JxW (i)*dx;                                 // x 
			IS_constants[c](2)  += fv_values.JxW (i)*dy;                                 // y
			IS_constants[c](3)  += fv_values.JxW (i)*dz;                                 // z
			IS_constants[c](4)  += fv_values.JxW (i)*(dx*dx);                       // x^2
			IS_constants[c](5)  += fv_values.JxW (i)*(dy*dy);                       // y^2
			IS_constants[c](6)  += fv_values.JxW (i)*(dz*dz);                       // z^2
			IS_constants[c](7)  += fv_values.JxW (i)*(dx*dy);                       // xy
			IS_constants[c](8)  += fv_values.JxW (i)*(dx*dz);                       // xz
			IS_constants[c](9)  += fv_values.JxW (i)*(dy*dz);                       // yz
			IS_constants[c](10)  += fv_values.JxW (i)*(dx*dx*dx);                // x^3
			IS_constants[c](11)  += fv_values.JxW (i)*(dy*dy*dy);                // y^3 
			IS_constants[c](12)  += fv_values.JxW (i)*(dz*dz*dz);                // z^3 
			IS_constants[c](13)  += fv_values.JxW (i)*(dx*dx*dy);                // x^2y
			IS_constants[c](14)  += fv_values.JxW (i)*(dx*dx*dz);                // x^2z
			IS_constants[c](15)  += fv_values.JxW (i)*(dy*dy*dx);                // xy^2
			IS_constants[c](16)  += fv_values.JxW (i)*(dy*dy*dz);                // zy^2
			IS_constants[c](17)  += fv_values.JxW (i)*(dz*dz*dx);                // xz^2
			IS_constants[c](18)  += fv_values.JxW (i)*(dz*dz*dy);                // yz^2
			IS_constants[c](19)  += fv_values.JxW (i)*(dx*dy*dz);                // xyz
			IS_constants[c](20) += fv_values.JxW (i)*(dx*dx*dx*dx);         // x^4 
			IS_constants[c](21) += fv_values.JxW (i)*(dy*dy*dy*dy);         // y^4
			IS_constants[c](22) += fv_values.JxW (i)*(dz*dz*dz*dz);         // z^4
			IS_constants[c](23) += fv_values.JxW (i)*(dx*dx*dx*dy);         // x^3y
			IS_constants[c](24) += fv_values.JxW (i)*(dx*dx*dx*dz);         // x^3y
			IS_constants[c](25) += fv_values.JxW (i)*(dy*dy*dy*dx);         // y^3x
			IS_constants[c](26) += fv_values.JxW (i)*(dy*dy*dy*dz);         // y^3z
			IS_constants[c](27) += fv_values.JxW (i)*(dz*dz*dz*dx);         // z^3x
			IS_constants[c](28) += fv_values.JxW (i)*(dz*dz*dz*dy);         // z^3y
			IS_constants[c](29) += fv_values.JxW (i)*(dx*dx*dy*dy);         // x^2y^2
			IS_constants[c](30) += fv_values.JxW (i)*(dx*dx*dz*dz);         // x^2z^2
			IS_constants[c](31) += fv_values.JxW (i)*(dy*dy*dz*dz);         // y^2z^2
		}
    }
    
    pcout << "Done!" << std::endl;  
}

#include "../include/Weno32.h" 


// Compute constants for smoothness indicators 

void Weno3_2D::compute_IS_constants() {
    
    std::cout << "Computing smoothness indicator constants" << std::endl;
    
    unsigned int N_gp = 5;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    Point<2> q_point;
	
	double x0, y0; 
	double x, y; 

    for (unsigned int c = 0; cell != endc; ++cell, ++c) {

		
		IS_constants[c] = 0.0; 
		
		IS_constants[c](0) = cell->measure();

		fv_values.reinit(cell);
		
		x0 = WENO_poly_consts[c](0); 
		y0 = WENO_poly_consts[c](1);

		for (unsigned int i = 0; i < N_gp*N_gp; i++) {

			q_point = fv_values.quadrature_point(i); 
			
			x = q_point(0); y = q_point(1); 

			IS_constants[c](1)  += fv_values.JxW (i)*(x-x0);                                // x 
			IS_constants[c](2)  += fv_values.JxW (i)*(y-y0);                                // y
			IS_constants[c](3)  += fv_values.JxW (i)*((x-x0)*(x-x0));                       // x^2
			IS_constants[c](4)  += fv_values.JxW (i)*((y-y0)*(y-y0));                       // y^2
			IS_constants[c](5)  += fv_values.JxW (i)*((x-x0)*(y-y0));                       // xy

		}
	
    }
    
    std::cout << "Done!" << std::endl;  
    std::cout << "===========================" << std::endl;
}

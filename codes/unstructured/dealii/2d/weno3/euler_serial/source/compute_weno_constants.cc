#include "../include/Weno32.h"


void Weno3_2D::compute_weno_polynomial_constants() {
    
    std::cout << "Computing WENO polynmial constants" << std::endl;

    unsigned int N_gp = 2;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    Point<2> q_point;

    double V0, x_0, y_0;

    for (unsigned int c = 0; cell != endc; ++cell, ++c) {

        V0 = cell->measure();

        fv_values.reinit(cell);
        
        WENO_poly_consts[c](0) = 0.0;
        WENO_poly_consts[c](1) = 0.0; 
        WENO_poly_consts[c](2) = 0.0;
        WENO_poly_consts[c](3) = 0.0;
        WENO_poly_consts[c](4) = 0.0;
        

        for (unsigned int i = 0; i < N_gp*N_gp; i++) {
            q_point = fv_values.quadrature_point(i);

            WENO_poly_consts[c](0) += (1./V0)*fv_values.JxW (i)*(q_point(0));
            WENO_poly_consts[c](1) += (1./V0)*fv_values.JxW (i)*(q_point(1));
        }
        
        x_0 = WENO_poly_consts[c](0); 
        y_0 = WENO_poly_consts[c](1);

        for (unsigned int i = 0; i < N_gp*N_gp; i++) {

            q_point = fv_values.quadrature_point(i);

            WENO_poly_consts[c](2) += (1./V0)*fv_values.JxW (i)*(q_point(0)-x_0)*(q_point(0)-x_0);
            WENO_poly_consts[c](3) += (1./V0)*fv_values.JxW (i)*(q_point(1)-y_0)*(q_point(1)-y_0);
            WENO_poly_consts[c](4) += (1./V0)*fv_values.JxW (i)*(q_point(0)-x_0)*(q_point(1)-y_0);
        }
        
    }
    
    std::cout << "Done!" << std::endl;
    std::cout << "===========================" << std::endl;
}

#include "../include/Weno32.h"

// Find smoothness indicators 

// Second order polynomial 

double compute_second_order_smoothness_indicator(Vector<double> coeffs, Vector<double> IS, double h) {

    double u_x   = coeffs(0); double u_y   = coeffs(1); 
    
    double const_term = IS(0)*(h*u_x*u_x + h*u_y*u_y);
    
    return (const_term); 
}

// Third Order polynomial 

double compute_third_order_smoothness_indicator(Vector<double> coeffs, Vector<double> IS,  double h) {
	
	double h3 = h*h*h; 
    
    double u_x   = coeffs(0); double u_y   = coeffs(1); double u_xx  = coeffs(2); double u_yy  = coeffs(3);
    double u_xy  = coeffs(4); 
    
    double const_term = IS(0)*(4*h3*u_xx*u_xx + h3*u_xy*u_xy + 4*h3*u_yy*u_yy + h*u_x*u_x + h*u_y*u_y);

    double x_term = IS(1)*(4*h*u_x*u_xx + 2*h*u_xy*u_y);

    double y_term = IS(2)*(2*h*u_x*u_xy + 4*h*u_y*u_yy); 

    double x2_term = IS(3)*(4*h*u_xx*u_xx + h*u_xy*u_xy);

    double y2_term = IS(4)*(h*u_xy*u_xy + 4*h*u_yy*u_yy); 

    double xy_term = IS(5)*(4*h*u_xx*u_xy + 4*h*u_xy*u_yy);
    
    return (const_term + x_term + y_term + x2_term + y2_term + xy_term); 

}



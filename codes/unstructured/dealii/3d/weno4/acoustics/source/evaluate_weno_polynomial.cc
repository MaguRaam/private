#include "../include/Weno432.h"


// Evaluate the WENO polynomial at a point 


double evaluate_weno_polynomial(const Vector<double>& coeffs,const Vector<double>& consts,const Point<3>& P, double h) {
    
    double x = P(0); double y = P(1); double z = P(2);
	
	double dx = (x-consts(0))/h; double dy = (y-consts(1))/h; double dz = (z-consts(2))/h;
    
    // Use Honers Algorith 
    
    double retval = coeffs(0) -  coeffs(4)*consts(3) - coeffs(5)*consts(4) - coeffs(6)*consts(5) -
			        coeffs(7)*consts(6) - coeffs(8)*consts(7) - coeffs(9)*consts(8) -  coeffs(10)*consts(9) - 
					coeffs(11)*consts(10) - coeffs(12)*consts(11) - coeffs(13)*consts(12) - coeffs(14)*consts(13) - 
					coeffs(15)*consts(14) - coeffs(16)*consts(15) - coeffs(17)*consts(16) - coeffs(18)*consts(17) - coeffs(19)*consts(18) +
                    dx*(coeffs(1) + dy*coeffs(7) + dz*coeffs(8) + dx*(coeffs(4) + dy*coeffs(13) + dz*coeffs(14) + dx*coeffs(10)) ) + 
                    dy*(coeffs(2) + dz*coeffs(9) + dy*(coeffs(5) + dx*coeffs(15) + dz*coeffs(16) + dy*coeffs(11) ) ) + 
					dz*(coeffs(3) + dz*(coeffs(6) + dx*coeffs(17) + dy*coeffs(18) + dz*coeffs(12) ) ) + dx*dy*dz*coeffs(19); 

    return retval; 
}

Vector<double> evaluate_conservative_gradient(const Vector<double>& coeffs, const Vector<double>& consts, Point<3> P, double h) {
    
    double x = P(0); double y = P(1); double z = P(2);
	
	double dx = (x-consts(0))/h; double dy = (y-consts(1))/h; double dz = (z-consts(2))/h;
	double h_inv = 1.0/h;
    
	Vector<double> grad(3);

    // Use Honers Algorith 

	// Gradient in x direction
    grad(0) = h_inv * (coeffs(1) + 2.0*dx*coeffs(4) + dy*coeffs(7) + dz*coeffs(8) + 3.0*dx*dx*coeffs(10) 
					+ 2.0*dx*dy*coeffs(13) + 2.0*dx*dz*coeffs(14) + dy*dy*coeffs(15) + dz*dz*coeffs(17) + dy*dz*coeffs(19) );

    grad(1) = h_inv * (coeffs(2) + 2.0*dy*coeffs(5) + dx*coeffs(7) + dz*coeffs(9) + 3.0*dy*dy*coeffs(11) 
					+ dx*dx*coeffs(13) + 2.0*dx*dy*coeffs(15) + 2.0*dy*dz*coeffs(16) + dz*dz*coeffs(18) + dx*dz*coeffs(19) );

    grad(2) = h_inv * (coeffs(3) + 2.0*dz*coeffs(6) + dx*coeffs(8) + dy*coeffs(9) + 3.0*dz*dz*coeffs(12) 
					+ dx*dx*coeffs(14) + dy*dy*coeffs(16) + 2.0*dx*dz*coeffs(17) + 2.0*dy*dz*coeffs(18) + dx*dy*coeffs(19) );

    return grad; 
}

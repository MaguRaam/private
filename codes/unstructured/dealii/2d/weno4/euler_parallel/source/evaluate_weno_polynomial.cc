#include "../include/Weno432.h"


// Evaluate the WENO polynomial at a point 


double evaluate_weno_polynomial(Vector<double> coeffs, Vector<double> consts,  Point<2> P, double h) {
    
    double x = P(0); double y = P(1); 
	
	double dx = (x-consts(0))/h; double dy = (y-consts(1))/h; 
    
    // Use Honers Algorith 
    
    double retval = coeffs(0) -  coeffs(3)*consts(2) - coeffs(4)*consts(3) - coeffs(5)*consts(4) - 
			        coeffs(6)*consts(5) - coeffs(7)*consts(6) - coeffs(8)*consts(7) - coeffs(9)*consts(8) +
                    dx*(coeffs(1) + dy*coeffs(5) + dx*(coeffs(3) + dy*coeffs(8) + dx*coeffs(6)) ) + 
                    dy*(coeffs(2) + dy*(coeffs(4) + dy*coeffs(7) + dx*coeffs(9) ) ); 

    return retval; 
}

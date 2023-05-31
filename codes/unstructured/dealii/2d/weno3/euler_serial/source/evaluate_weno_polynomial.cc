#include "../include/Weno32.h"


// Evaluate the WENO polynomial at a point 


double evaluate_weno_polynomial(Vector<double> coeffs, Vector<double> consts,  Point<2> P) {
    
    double x = P(0); double x0 = consts(0); 
    double y = P(1); double y0 = consts(1); 

    return (coeffs(0) + 
            coeffs(1)*(x - x0) + 
            coeffs(2)*(y - y0) +
            coeffs(3)*((x-x0)*(x-x0) - consts(2)) + 
            coeffs(4)*((y-y0)*(y-y0) - consts(3)) + 
            coeffs(5)*((x-x0)*(y-y0) - consts(4)));
} 

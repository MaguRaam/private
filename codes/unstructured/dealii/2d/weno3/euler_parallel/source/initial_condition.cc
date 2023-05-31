#include "../include/Weno32.h"


Vector<double> initial_condition(Point <2> P) {
	
    double x = P(0); double y = P(1); 
    
    double a; 
    
	Vector<double> W(4); 
    
    if (y <= 1./2.) {
        W(0) = 2.0;
        W(1) = 0.0; 
        W(3) = 2.0*y+1.0; 
        a = std::sqrt(GAMMA*W(3)/W(0));
        W(2) = -0.025*a*cos(8.0*M_PI*x);
    }
    
    else {
        W(0) = 1.0;
        W(1) = 0.0; 
        W(3) = y+1.5; 
        a = std::sqrt(GAMMA*W(3)/W(0));
        W(2) = -0.025*a*cos(8.0*M_PI*x);
    }
    
    return W;
}

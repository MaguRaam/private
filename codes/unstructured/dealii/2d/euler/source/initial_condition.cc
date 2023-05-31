#include "../include/Euler.h"

// Set the initial condtions  

Vector<double> initial_condition(Point <2> P) {
    
	
    W(0) = 1.4; 
    W(1) = M; 
    W(2) = 0.0; 
    W(3) = 1.0; 
    
	
	
    return W; 
}

#include "../include/Weno432.h"

Vector<double> initial_condition(const Point <3>& p, const double h) {
	
	Vector<double> W(5);

	double x = p(0), y = p(1), z = p(2);

	W(0) = 1.0;
//	W(1) = cos(x*2.0*M_PI)*cos(y*2.0*M_PI)*cos(z*2.0*M_PI);
	W(1) = sin(x*2.0*M_PI)*sin(y*2.0*M_PI)*sin(z*2.0*M_PI);
	W(2) = sin(x*2.0*M_PI)*sin(y*2.0*M_PI)*sin(z*2.0*M_PI);
	W(3) = sin(x*2.0*M_PI)*sin(y*2.0*M_PI)*sin(z*2.0*M_PI);
	W(4) = 5.0;


/*	
   	W(0) = 1.4;
   	W(1) = M;
   	W(2) = 0.0;
   	W(3) = 0.0;
   	W(4) = 1.0;
*/
	return W;
}

#include "../include/Weno432.h"

Vector<double> initial_condition(const Point <3>& p, const double h) {
	
	Vector<double> U(2);

	double x = p(0), y = p(1), z = p(2);

	U(0) = 0.0;
	U(1) = 0.0;

	return U;
}

#include "../include/Weno432.h"

Vector<double> initial_condition(Point <2> P) {
    
	// finalTime = 0.2;

	double x = P(0); double y = P(1);
	Vector<double> W(4);

	double rho[] = { 1.0, 2.0, 1.0625, 0.5313};
	double u[]   = { 0.0, 0.0,    0.0,    0.0};
	double v[]   = {-0.3, 0.3, 0.8145, 0.4276};
	double p[]   = { 1.0, 1.0,    0.4,    0.4};


	// Zone 1

	if (x >= 0.5 && y >= 0.5) {
		W[0] = rho[0];
		W[1] =   u[0];
		W[2] =   v[0];
		W[3] =   p[0];
	}

	// Zone2

	if (x < 0.5 && y >= 0.5) {
		W[0] = rho[1];
		W[1] =   u[1];
		W[2] =   v[1];
		W[3] =   p[1];
	}

	// Zone 3

	if (x < 0.5 && y < 0.5) {
		W[0] = rho[2];
		W[1] =   u[2];
		W[2] =   v[2];
		W[3] =   p[2];
	}

	// Zone 4
    
	if (x >= 0.5 && y < 0.5) {
		W[0] = rho[3];
		W[1] =   u[3];
		W[2] =   v[3];
		W[3] =   p[3];
    }

    
    return W;
}

/*
 * quadrature.cc
 *      Author: sunder
 */

#include "quadrature.hh"

//----------------------------------------------------------------------------
// Get n point Gauss-Legendre quadrature points and weights in the
// interval [x1,x2]
// Source: Numerical Recipies, Press et. al., Edition 3
//----------------------------------------------------------------------------

void QGauss(double x1, double x2, Array1D<double>& x, Array1D<double>& w, int n) {

	// Get the quadrature points by solving the Legendre equations

	const double EPS = 1.0e-14;
	double z1, z, xm, xl, pp, p3, p2, p1;
	int m = (n+1)/2;
	xm = 0.5*(x2 + x1);
	xl = 0.5*(x2 - x1);
	for (int i=0;i<m;i++) {
		z = cos(M_PI*(i + 0.75)/(n + 0.5));
		do {
			p1 = 1.0;
			p2 = 0.0;
			for (int j = 0;j < n;j++) {
				p3 = p2;
				p2 = p1;
				p1 = ((2.0*j + 1.0)*z*p2 -j*p3)/(j+1);
			}
			pp = n*(z*p1 - p2)/(z*z-1.0);
			z1=z;
			z=z1-p1/pp;
		} while (std::abs(z-z1) > EPS);
		x[i] = xm-xl*z;
		x[n-1-i] = xm+xl*z;
		w[i] = 2.0*xl/((1.0-z*z)*pp*pp);
		w[n-1-i] = w[i];
	}
}

//----------------------------------------------------------------------------
// Except for QGaussTriangle, the source for all other quadrature formulae
// for triangles is: QUADRATURE_RULES_TRI
// https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tri/quadrature_rules_tri.html
//----------------------------------------------------------------------------

//----------------------------------------------------------------------------
// Centroid: order 1, degree of precision 1
//----------------------------------------------------------------------------

void QCentroid(Array1D<double>& xi, Array1D<double>& eta, Array1D<double>& w, int& N) {
	N = 1;

	xi.reinit(N);
	eta.reinit(N);
	w.reinit(N);

	 xi[0] = 0.33333333333333333333;
	eta[0] = 0.33333333333333333333;
	  w[0] = 1.00000000000000000000;
}



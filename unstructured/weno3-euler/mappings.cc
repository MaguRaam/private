/*
 * mappings.cc
 *      Author: sunder
 */

#include "mappings.hh"


//-----------------------------------------------------------------------
// Given a point in reference coordinates (xi,eta) on the triangle
// (x1,y1), (x2,y2), (x3,y3), get physical coordinates (x,y)
//-----------------------------------------------------------------------

void reference_to_physical(double x1, double y1, double x2, double y2, double x3, double y3,
		                 double xi, double eta, double& x, double& y) {

    x = x1 + (x2 - x1)*xi + (x3 - x1)*eta;
    y = y1 + (y2 - y1)*xi + (y3 - y1)*eta;

}

//-----------------------------------------------------------------------
// Given a point in physical coordinates (x,y) on the triangle
// (x1,y1), (x2,y2), (x3,y3), get reference coordinates (xi,eta)
//-----------------------------------------------------------------------


void physical_to_reference(double x1, double y1, double x2, double y2, double x3, double y3,
		                 double x, double y, double& xi, double& eta) {

    double r1_detJ = 1.0/((x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1));

    xi  = r1_detJ*(x3*y1 - x1*y3 + x*(y3 - y1) + y*(x1 - x3));
    eta = r1_detJ*(x1*y2 - x2*y1 + x*(y1 - y2) + y*(x2 - x1));
}

//-----------------------------------------------------------------------
// Given a point on the interval [0,1], find the mapped point on the
// line joining (1,0) and (0,1). This is useful for finding quadrature
// points on the inclined face of the reference element
//-----------------------------------------------------------------------

void straight_face_to_inclined_face(double u, double& xi, double& eta) {

	xi  = 1.0 - u;
    eta = u;
}

//-----------------------------------------------------------------------
// For a point in the reference plane (xi,eta), get the value of k'th
// basis function. These basis functions are orthogonal on a unit
// triangle
//-----------------------------------------------------------------------

double basis(double xi,double eta, int k) {
	double phi = 0.0;

	switch (k) {
	case 0:
		phi = 1.0;
		break;
	case 1:
		phi = 2.0*xi - 1.0 + eta;
		break;
	case 2:
		phi = -1.0 + 3.0*eta;
		break;
	case 3:
		phi = 1.0 - 2.0*eta + eta*eta  - 6.0*xi + 6.0*xi*eta + 6.0*xi*xi;
		break;
	case 4:
		phi = 5.0*eta*eta + 10.0*xi*eta - 6.0*eta - 2.0*xi + 1.0;
		break;
	case 5:
		phi = 1.0 - 8.0*eta + 10.0*eta*eta;
		break;
	default:
		phi = 0.0;
		std::cerr << "Basis not implemented" << std::endl;
		std::exit(EXIT_FAILURE);
		break;
	}

	return phi;
}


//-----------------------------------------------------------------------
// For a triangle in the reference plane, having coordinate (xi0,eta0),
// (xi1,eta1), (xi2,eta2) - find the average of k'th basis function
// over the triangle. This function is very useful for constructing the
// reconstruction matrices
//-----------------------------------------------------------------------

double basis_average_over_triangle(double x1, double y1, double x2, double y2, double x3, double y3, int k) {

	// Quadrature points and weights on reference triangle

	int q;
	const int N_gp = 6; // STRANG5, order 6, degree of precision 4

	double xiGP[]  = {0.816847572980459,0.091576213509771,0.091576213509771,0.108103018168070,0.445948490915965,0.445948490915965};
	double etaGP[] = {0.091576213509771,0.816847572980459,0.091576213509771,0.445948490915965,0.108103018168070,0.445948490915965};
	double wGP[]   = {0.109951743655322,0.109951743655322,0.109951743655322,0.223381589678011,0.223381589678011,0.223381589678011};

	double xi, eta;
	double integral = 0.0;


	for (q = 0; q < N_gp; ++q) {
		reference_to_physical(x1,y1,x2,y2,x3,y3, xiGP[q],etaGP[q],xi,eta);

		integral += wGP[q]*basis(xi,eta,k);
	}

	return integral;
}

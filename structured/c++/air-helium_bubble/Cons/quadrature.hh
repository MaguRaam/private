/*
 * quadrature.hh
 *      Author: sunder
 */

#ifndef QUADRATURE_HH_
#define QUADRATURE_HH_

#include "headers.hh"

// Gauss Quadrature weights and other transforms

//template<typename T>
//void isoparametric_transform(T, T, const std::vector<T>&, const std::vector<T>&, T&,  T&);

// Mid-point rule
template<typename T>
void QMidPoint(std::vector<T>&, std::vector<T>&, int);

// Trapezoidal rule
template<typename T>
void QTrapezoid(std::vector<T>&, std::vector<T>&, int);

// Simpson's 1/3 rule
template<typename T>
void QSimpson(std::vector<T>&, std::vector<T>&, int);

// Simpson's 3/8 rule
template<typename T>
void QSimpson3_8(std::vector<T>&, std::vector<T>&, int);

// Gauss - Legendre Quadrature
template<typename T>
void QGauss(std::vector<T>&, std::vector<T>&, int);

// 2D Isoparametric transformation

template<typename T>
void isoparametric_transform(T u, T v, const std::vector<T>&  xi, const std::vector<T>&  yi, T& x,  T& y) {

	// Shape functions
	T N[4];

    N[0] = (0.5-u)*(0.5-v);
    N[1] = (0.5+u)*(0.5-v);
    N[2] = (0.5-u)*(0.5+v);
    N[3] = (0.5+u)*(0.5+v);

    x = 0.0; y = 0.0;

    for (unsigned int i = 0; i < 4; i++) {
        x += xi[i]*N[i];
        y += yi[i]*N[i];
    }
}

// Mid-point rule

template<typename T>
void QMidPoint(std::vector<T>& x, std::vector<T>& w, int n) {

	assert(n==1);

    w[0] =  1.0; x[0] =  0.0;

}

// Trapezoidal rule

template<typename T>
void QTrapezoid(std::vector<T>& x, std::vector<T>& w, int n) {

	assert(n == 2);

	w[0] =  1./2.; w[1] = 1./2.;
	x[0] = -1./2.; x[1] = 1./2.;
}

// Simpson's 1/3 Rule

template<typename T>
void QSimpson(std::vector<T>& x, std::vector<T>& w, int n) {

	assert(n==3);

	w[0] =  1./6.; w[1] = 2./3.; w[2] = 1./6.;
	x[0] = -1./2.; x[1] =   0.0; x[2] = 1./2.;
}

// Simpson's 3/8 Rule

template<typename T>
void QSimpson3_8(std::vector<T>& x, std::vector<T>& w, int n) {

	assert(n==4);

	w[0] =  1./8.; w[1] =  3./8.; w[2] = 3./8.; w[3] = 1./8.;
	x[0] = -1./2.; x[1] = -1./6.; x[2] = 1./6.; x[3] = 1./2.;
}


// Gauss-Legendre Quadrature

template<typename T>
void QGauss(std::vector<T>& x, std::vector<T>& w, int n) {

	switch (n){
	case 1:
		w[0] = 1.0000000000000000; x[0] =  0.00000000000000000;
		break;
	case 2:
        w[0] = 0.5000000000000000; x[0] = -0.28867513459481287;
        w[1] = 0.5000000000000000; x[1] =  0.28867513459481287;
        break;
	case 3:
        w[0] = 0.4444444444444444; x[0] =  0.0000000000000000;
        w[1] = 0.2777777777777778; x[1] = -0.3872983346207417;
        w[2] = 0.2777777777777778; x[2] =  0.3872983346207417;
        break;
	case 4:
        w[0] = 0.3260725774312730; x[0] =  0.1699905217924282;
        w[1] = 0.3260725774312730; x[1] = -0.1699905217924282;
        w[2] = 0.1739274225687270; x[2] =  0.4305681557970263;
        w[3] = 0.1739274225687270; x[3] = -0.4305681557970263;
        break;
	case 5:
        w[0] = 0.28444444444444444; x[0] =  0.0000000000000000;
        w[1] = 0.23931433524968326; x[1] = -0.2692346550528416;
        w[2] = 0.23931433524968326; x[2] =  0.2692346550528416;
        w[3] = 0.11846344252809456; x[3] = -0.4530899229693320;
        w[4] = 0.11846344252809456; x[4] =  0.4530899229693320;
        break;
	case 6:
        w[0] = 0.1803807865240693; x[0] =  0.33060469323313224;
        w[1] = 0.1803807865240693; x[1] = -0.33060469323313224;
        w[2] = 0.2339569672863455; x[2] =  0.11930959304159845;
        w[3] = 0.2339569672863455; x[3] = -0.11930959304159845;
        w[4] = 0.0856622461895852; x[4] =  0.46623475710157600;
        w[5] = 0.0856622461895852; x[5] = -0.46623475710157600;
        break;
	default:
		// Get the quadrature points by solving the Legendre equations
		T x1 = -0.5; T x2 = 0.5;
		const T EPS = 1.0e-14;
		T z1, z, xm, xl, pp, p3, p2, p1;
		unsigned int m = (n+1)/2;
		xm = 0.5*(x2 + x1);
		xl = 0.5*(x2 - x1);
		for (unsigned int i=0;i<m;i++) {
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
}

#endif /* QUADRATURE_HH_ */

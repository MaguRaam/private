/*
 * pde.cc
 *      Author: sunder
 */

#include "hype.hh"

//----------------------------------------------------------------------------
// Convert a conserved variable to a primitive variable
// (Density,momentum,total energy) -> (Density,velocity,pressure)
//----------------------------------------------------------------------------

void PDECons2Prim(const Array1D<double>& Q, Array1D<double>& V) {

	double irho = 1.0/Q[0];

	V[0] = Q[0];
	V[1] = irho*Q[1];
	V[2] = irho*Q[2];
	V[3] = (GAMMA -1.0)*( Q[3] - 0.5*irho*(Q[1]*Q[1] + Q[2]*Q[2]) );
}

//----------------------------------------------------------------------------
// Convert a primitive variable to conserved variable
// (Density,velocity,pressure) -> (Density,momentum,total energy)
//----------------------------------------------------------------------------

void PDEPrim2Cons(const Array1D<double>& V, Array1D<double>& Q) {

	double e = (V[3])/(GAMMA - 1.0);
	double k = 0.5*V[0]*(V[1]*V[1] + V[2]*V[2]);

	Q[0] = V[0];
	Q[1] = V[0]*V[1];
	Q[2] = V[0]*V[2];
	Q[3] = k + e;
}

//----------------------------------------------------------------------------
// Conservative part of the flux. The input should be conserved variable.
// This function should also return the maximum absolute
// eigenvalue of the PDE for the given input state.
//----------------------------------------------------------------------------

double PDEFlux(const Array1D<double>& Q, double nx, double ny, Array1D<double>& F) {

	double rho = Q[0];
	double irho = 1.0/rho;
	double u = irho*Q[1];
	double v = irho*Q[2];
	double E = Q[3];
    double p = (GAMMA -1.0)*( E - 0.5*rho*(u*u + v*v) );
    double un = u*nx+v*ny;

    if (rho < rho_floor) {
        std::cerr << "Negative density, density = " <<  rho << std::endl;
        std::exit(EXIT_FAILURE);
    }

    if (p < prs_floor) {
        std::cerr << "Negative pressure, pressure = " << p << std::endl;
        std::exit(EXIT_FAILURE);
    }

    F[0] = rho*un;
    F[1] = rho*u*un + nx*p;
    F[2] = rho*v*un + ny*p;
    F[3] = un*(E + p);

    // Return the maximum eigenvalue

    return (std::abs(un) + std::sqrt(irho*GAMMA*p));
}

//----------------------------------------------------------------------------
// LLF Riemann solver to find upwind flux at faces
//----------------------------------------------------------------------------

double LLFRiemannSolver(const Array1D<double>& QL, const Array1D<double>& QR, double nx, double ny, Array1D<double>& Flux) {

	int iVar;

	Array1D<double> FL(nVar), FR(nVar);

	double smaxl = PDEFlux(QL, nx, ny, FL);
	double smaxr = PDEFlux(QR, nx, ny, FR);

	double smax = std::max(smaxl, smaxr);

	for (iVar = 0; iVar < nVar; ++iVar) {
		Flux[iVar] = 0.5*( FL[iVar] + FR[iVar] - smax*(QR[iVar] - QL[iVar]) );
	}

	return smax;
}



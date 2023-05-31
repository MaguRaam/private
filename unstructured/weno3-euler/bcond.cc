/*
 * bcond.cc
 *      Author: sunder
 */

#include "hype.hh"

// Commonly used boundary conditions
void slip_wall(const Array1D<double>&, double, double, Array1D<double>&);
void supersonic_outflow(const Array1D<double>&, Array1D<double>&);
void subsonic_outflow(const Array1D<double>&, double, Array1D<double>&);
void free_stream(const Array1D<double>&, Array1D<double>&);

//-----------------------------------------------------------------------
// Get the right hand side of the boundary
//-----------------------------------------------------------------------

void get_right_state(const Array1D<double>& QL, const Array1D<double>& Vb,
		int bcond, double nx, double ny, Array1D<double>& QR) {

	switch (bcond) {
	case 1:
		supersonic_outflow(QL, QR);
		break;
	case 2:
		slip_wall(QL,nx,ny,QR);
		break;
	case 3:
		subsonic_outflow(QL,Vb[3],QR);
		break;
	case 4:
		free_stream(Vb,QR);
		break;
	default:
		supersonic_outflow(QL, QR);
		//std::cerr << "Boundary condition " << bcond << " not implemented. Exiting" << std::endl;
		//std::exit(EXIT_FAILURE);
	}
}

//-----------------------------------------------------------------------
// Slip wall boundary condition. Reflect the momentum to make normal
// component zero
//-----------------------------------------------------------------------

void slip_wall(const Array1D<double>& QL, double nx, double ny, Array1D<double>& QR) {

	double un = QL[1]*nx + QL[2]*ny;

	QR[0] = QL[0];
	QR[1] = QL[1] - un*nx;
	QR[2] = QL[2] - un*ny;
	QR[3] = QL[3];
}

//-----------------------------------------------------------------------
// Supersonic outflow or transmissive or continuative boundary conditions
// Copy values from inside the domain.
//-----------------------------------------------------------------------

void supersonic_outflow(const Array1D<double>& QL, Array1D<double>& QR) {

	QR[0] = QL[0];
	QR[1] = QL[1];
	QR[2] = QL[2];
	QR[3] = QL[3];
}

//-----------------------------------------------------------------------
// Subsonic outflow (Back pressure)
//-----------------------------------------------------------------------

void subsonic_outflow(const Array1D<double>& QL, double p_inf,  Array1D<double>& QR) {

	Array1D<double> VL(nVar);
	Array1D<double> VR(nVar);

	PDECons2Prim(QL, VL);

	VR[0] = VL[0];
	VR[1] = VL[1];
	VR[2] = VL[2];
	VR[3] = p_inf; // Only fix the pressure

	PDEPrim2Cons(VR, QR);
}

//-----------------------------------------------------------------------
// Freestream
//-----------------------------------------------------------------------

void free_stream(const Array1D<double>& Vb, Array1D<double>& QR) {
	PDEPrim2Cons(Vb, QR);
}

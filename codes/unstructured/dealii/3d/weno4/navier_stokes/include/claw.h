/*
 * claw.h
 *      Author: sunder
 */

#ifndef CLAW_H_
#define CLAW_H_

#include "Headers.h"

void Allocate_2D_R(double**& m, int d1, int d2);
void Factor(  double** a, int n, int*& npivot,  double& det, bool& sing);
void Solve(  double** a, int n,   double*& b, int*& npivot);
void LUInverse(  double**& a_inv,   double** a, int n);

class Conservation_Law {
	double specific_heat_cp;
	double specific_heat_cv;
	double prandtl_no;
	double specific_heat_ratio;

public:
	/* Constructors */

	Conservation_Law();
	Conservation_Law(double, double, double);

	/* Physical properties of the gas */

	double c_p() const;
	double c_v() const;
	double gamma() const;
	double R() const;
	double Pr() const;
	double viscosity(const Vector<double>&);
	double viscosity(const double);
	double heat_conductivity(const Vector<double>&);
	double heat_conductivity(const double);


	/* Auxiliary properties of the gas */

	double Pressure(const Vector<double>&);
	double Temperature(const Vector<double>&);
	double speed_of_sound(const Vector<double>&);

	/* Conversions and flux calculations */

//	void primitive_to_conserved(const Vector<double>&, Vector<double>&);
	void primitive_to_conserved(const Vector<double>& W, Vector<double>& U) const;
	void conserved_to_primitive(const Vector<double>&, Vector<double>&) const;
	void Usual_To_Conservative(double*& Var_Q, double* Var_U, double gamma_l, double P_Infty_l);

	void compute_pure_convective_flux(const Vector<double>&, const double, const double, const double, const double, Vector<double>&);
	void compute_pure_viscous_flux(const Vector<double>&, const FullMatrix<double>&, const double, const double, Vector<double>&);
	void compute_viscous_flux_zero_flux(const Vector<double>&, const FullMatrix<double>&, const double, const double, Vector<double>&);

	void FluxFunction(Vector<double>& FV, double* Q, double nx, double ny, double nz);

	void temperature_gradient(const Vector<double>&, const FullMatrix<double>&, double*);

	void compute_x_y_flux(const Vector<double>&, const FullMatrix<double>&, double, double, Vector<double>&, Vector<double>&);
	void compute_physical_flux(const Vector<double>&, const FullMatrix<double>&, double, double, double, double, Vector<double>&);


	/* Riemann solvers  */
	Vector<double> local_Lax_Friedrichs_riemann_solver(const Vector<double>& UL, const Vector<double>& UR, double nx, double ny, double nz, const Point<3>& P, bool boundary);

	Vector<double> HLLC_riemann_solver(const Vector<double>& UL, const Vector<double>& UR, double nx, double ny, double nz, const Point<3>& P, bool boundary);

	Vector<double> rotated_HLLC_riemann_solver(const Vector<double>& WL, const Vector<double>& WR, double nx, double ny, double nz, const Point<3>& P, bool boundary);

	void compute_normal_viscous_flux(const Vector<double>& U,
		const std::vector<Vector<double> >& gradU,
		double nx, double ny, double nz, 
		Vector<double>& Flux);

	Vector<double> viscous_reiman_solver(const Vector<double>& UL, const Vector<double>& UR, const std::vector<Vector<double> >& gradUL, const std::vector<Vector<double> >& gradUR, double nx, double ny, double nz);

	Vector<double> NS_riemann_solver(const Vector<double>& UL, const Vector<double>& UR, const std::vector<Vector<double> >& gradUL, const std::vector<Vector<double> >& gradUR, double nx, double ny, double nz, Point<3> P, bool boundary);


	/* Boundary conditions */


	Vector<double> compute_convective_flux_at_no_slip_boundary(const Vector<double>&, double, double,  double = 0.0, double = 0.0);
	Vector<double> compute_viscous_flux_at_no_slip_adiabatic_boundary(const Vector<double>&, const FullMatrix<double>&, double, double, double =0.0, double = 0.0);

	Vector<double> compute_convective_flux_at_inflow_outflow_boundary(const Vector<double>&, const Vector<double>&, double, double);
	Vector<double> compute_viscous_flux_at_inflow_outflow_boundary(const Vector<double>&, const FullMatrix<double>&, const Vector<double>&, double, double);

	void additional_data(
		const Vector<double>& U,
		const std::vector< Vector<double> >& grad_U,
		double& p, std::vector< Vector<double> >& stress);

};





#endif /* CLAW_H_ */

/*
 * riemann_solver.cc
 *      Author: sunder
 */  

#include "riemann_solver.hh"

/* Definitions for the Riemann_Solver class */

//----------------------------------------------------------------------------
// Default constructor: Initializes both the phases to ideal gas with
// specific heat ratio 1.4 
//----------------------------------------------------------------------------

Riemann_Solver::Riemann_Solver():
	p1(0.0),
	p2(0.0),
	g1(1.4),
	g2(1.4)
{}


//----------------------------------------------------------------------------
// Main construtor: Takes the stiffness constants of the phases and their 
// specific heat ratios 
//----------------------------------------------------------------------------

Riemann_Solver::Riemann_Solver(const double& p1_, const double& p2_, const double& g1_, const double& g2_):
	p1(p1_),
	p2(p2_),
	g1(g1_),
	g2(g2_)
{}

//----------------------------------------------------------------------------
// Convert conserved variables to primitive variables 
//----------------------------------------------------------------------------

void Riemann_Solver::conserved_to_primitive(const std::vector<double>& Q, std::vector<double>& V) {

	double g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-Q[4])*(g1 -1.0) + Q[4]*(g2-1.0));
	double P_inf = ((g -1.0)/g)*( g1*p1*Q[4]/(g1 - 1.0) + g2*p2*(1.0 - Q[4])/(g2 - 1.0) );

	V[0] = Q[0];
    V[1] = Q[1]/Q[0];
	V[2] = Q[2]/Q[0];
    V[3] = (g -1.0)*( Q[3] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2])/Q[0] )  - g*P_inf;
    V[4] = Q[4];
}

//----------------------------------------------------------------------------
// Convert primitive variables to conserved variables 
//----------------------------------------------------------------------------

void Riemann_Solver::primitive_to_conserved(const std::vector<double>& V, std::vector<double>& Q) {
	
	
	double g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-V[4])*(g1 -1.0) + V[4]*(g2-1.0));
	double P_inf = ((g -1.0)/g)*( g1*p1*V[4]/(g1 - 1.0) + g2*p2*(1.0 - V[4])/(g2 - 1.0) );

    double e = (V[3] + g*P_inf)/(g - 1.0);
    double k = 0.5*V[0]*(V[1]*V[1] + V[2]*V[2]);

    Q[0] = V[0];
    Q[1] = V[0]*V[1];
    Q[2] = V[0]*V[2];
    Q[3] = k + e;
	Q[4] = V[4];
}

//----------------------------------------------------------------------------
// Find the conservative flux components F_x and F_y
//----------------------------------------------------------------------------

double Riemann_Solver::conservative_flux(const std::vector<double>& Q, 
										 const double& nx, const double& ny,
										 const double& x, const double& y,
										 std::vector<double>& F) {
	double g = 1.0 + (g1-1.0)*(g2-1.0)/((1.0-Q[4])*(g1 -1.0) + Q[4]*(g2-1.0));
	double P_inf = ((g -1.0)/g)*( g1*p1*Q[4]/(g1 - 1.0) + g2*p2*(1.0 - Q[4])/(g2 - 1.0) );
	double p = (g -1.0)*( Q[3] - 0.5*(Q[1]*Q[1] + Q[2]*Q[2])/Q[0] )  - g*P_inf;

	if (p < 0.0) {
		std::cerr << "Neagative pressure, p = " << p << std::endl;
		std::cerr << "At x = " << x << ", y = " << y << std::endl;  
		std::exit(EXIT_FAILURE);
	}

	if (Q[0] < 0.0) {
		std::cerr << "Neagative density, rho = " <<  Q[0] << std::endl;
		std::cerr << "At x = " << x << ", y = " << y << std::endl; 
		std::exit(EXIT_FAILURE);
	}

	/*

	if (Q[4] < 0.0 || Q[4] > 1.0) {
		std::cerr << "Phi out of range, phi = " <<  Q[4] << std::endl;
		std::cerr << "At x = " << x << ", y = " << y << std::endl; 
		std::exit(EXIT_FAILURE);
	}
	*/

	F[0] = nx*Q[1] + ny*Q[2];
	F[1] = nx*(Q[1]*Q[1]/Q[0] + p) + ny*(Q[1]*Q[2]/Q[0]);
	F[2] = nx*(Q[1]*Q[2]/Q[0]) + ny*(Q[2]*Q[2]/Q[0] + p);
	F[3] = nx*(Q[1]*(Q[3] + p)/Q[0]) + ny*(Q[2]*(Q[3] + p)/Q[0]);
	F[4] = 0.0;
	
	// Also obtain the maximum eigen value 
	
	double s_max = std::abs(Q[1]*nx + Q[2]*ny)/Q[0] + std::sqrt(g*(p + P_inf)/Q[0]);
	
	
	return s_max;
}

//----------------------------------------------------------------------------
// Find the non-conservative flux components nF
//----------------------------------------------------------------------------

void Riemann_Solver::non_conservative_flux(const std::vector<double>& Q,
									   const std::vector<double>& grad_Q_x,
									   const std::vector<double>& grad_Q_y,
									      std::vector<double>& nF) {

	double phi_x = grad_Q_x[4]; double phi_y = grad_Q_y[4];
	
	// Now find the fluxes 
	
	nF[0] = 0.0;
	nF[1] = 0.0;
	nF[2] = 0.0; 
	nF[3] = 0.0; 
	nF[4] = (Q[1]*phi_x + Q[2]*phi_y)/Q[0]; 
}

//----------------------------------------------------------------------------
// Find the non-conservative matrix B
//----------------------------------------------------------------------------

void Riemann_Solver::matrix_B(const std::vector<double>& Q, const double& nx, const double& ny, double** B) {
	
	for (int i = 0; i < n_comp; ++i) 
		for (int j = 0; j < n_comp; ++j)
			B[i][j] = 0.0; 
		
	B[4][4] =  nx*Q[1] + ny*Q[2];       
}

//------------------------------------------------------------------------------------------------------------------
// No need to change anything beyond this point 
//-------------------------------------------------------------------------------------------------------------------

//----------------------------------------------------------------------------
// The Roe matrix B̃(Ua,Ub) between two generic states Qa and Qb is numerically 
// via Gaussian quadrature, by the function RoeMatrix
//----------------------------------------------------------------------------

void Riemann_Solver::roe_matrix(const std::vector<double>& Qa, const std::vector<double>& Qb,
							const double& nx, const double& ny,
							double** BRoe) {

		// Gauss quadrature points in the interval [0,1] 
	
	const int N_gp = 3;
	double s_gp[] = {0.5-std::sqrt(15.)/10., 0.5, 0.5+std::sqrt(15.)/10.};
	double w_gp[] = {5./18., 8./18., 5./18.};
	
	double** B; Utilities::allocate_mem_2d(B, n_comp, n_comp);

	std::vector<double> Q(n_comp); 
	
	// First make the BRoe matrix zero
	
	for (int i = 0; i < n_comp; ++i) {
		for (int j = 0; j < n_comp; ++j) {
			BRoe[i][j] = 0.0; 
		}
	}
	
	for (int q = 0; q < N_gp; ++q) {
		
		for (int c = 0; c < n_comp; ++c) 
			Q[c] = Qa[c] + s_gp[q]*(Qb[c] - Qa[c]);
		
		matrix_B(Q, nx, ny, B);
		
		for (int i = 0; i < n_comp; ++i) {
			for (int j = 0; j < n_comp; ++j) {
				BRoe[i][j] += w_gp[q]*B[i][j];	
			}
		}
	}
	
	Utilities::free_mem_2d(B, n_comp);
}

//----------------------------------------------------------------------------
// LLF Riemann solver
//----------------------------------------------------------------------------

void mat_vec_mult(double** A, const std::vector<double>& x, std::vector<double>& b, int size) {
	
	for (int i = 0; i < size; ++i) { 
		b[i] = 0.0; 
		
		for (int j = 0; j < size; ++j) {
			b[i] += A[i][j]*x[j];
		}
		
	}
}

double Riemann_Solver::llf_riemann_solver(const std::vector<double>& QL, 
								      const std::vector<double>& QR,
									  const double& nx, 
									  const double& ny,
									  const double& x,
									  const double& y,
								      std::vector<double>& Flux, std::vector<double>& D) {

	
	
	std::vector<double> FL(n_comp); std::vector<double> FR(n_comp); 
	
	double s_max_l = conservative_flux(QL, nx, ny, x, y, FL); 
	double s_max_r = conservative_flux(QR, nx, ny, x, y, FR);
	
	double s_max = std::max(s_max_l, s_max_r);
	
	// Conservative jump 
	
	std::vector<double> Q_jump(n_comp); 
	
	for (int i = 0; i < n_comp; ++i) {
		Q_jump[i] = QR[i] - QL[i]; 
		Flux[i] = 0.5*(FR[i] + FL[i] - s_max*Q_jump[i]);
	}
	
	// Path conservative jump 
	
	double** B; Utilities::allocate_mem_2d(B, n_comp, n_comp); 

	roe_matrix(QL, QR, nx, ny, B); 

	// Find the diffusive flux 
	
	mat_vec_mult(B, Q_jump, D, n_comp);
	
	Utilities::free_mem_2d(B, n_comp);
	
	return s_max; 
}

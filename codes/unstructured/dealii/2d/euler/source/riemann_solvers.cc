#include "../include/Riemann_Solvers.h"
#include "../include/Headers.h"
#include "../include/Exceptions.h"

// Riemann Solver Function definitions  

// Simple flux functions and conversions 

// ============================================================================
Vector<double> conserved_to_primitive(Vector<double> U) {

    Vector<double> W(4); 
    
    W[0] = U[0];
    W[1] = U[1]/U[0];
    W[2] = U[2]/U[0]; 
    W[3] = (GAMMA-1.0)*(U[3] - 0.5*((U[1]*U[1] + U[2]*U[2])/U[0]));
    
    return W; 
}

Vector<double> primitive_to_conserved(Vector<double> W) {
    Vector<double> U(4);

    double e = W[3]/((GAMMA - 1.0)*W[0]);   // Internal energy, e = p/((gamma-1)*rho)
    double k = 0.5*(W[1]*W[1] + W[2]*W[2]); // Kinetic energy, k = (1/2)*(u^2 + v^2)

    U[0] = W[0];
    U[1] = W[0]*W[1];
    U[2] = W[0]*W[2];
    U[3] = W[0]*(k + e);

	return U;
}

Vector<double> flux_primitive(Vector<double> W) {

    Vector<double> F(4);

    double e = W[3]/((GAMMA - 1.0)*W[0]);   // Internal energy, e = p/((gamma-1)*rho)
    double k = 0.5*(W[1]*W[1] + W[2]*W[2]); // Kinetic energy, k = (1/2)*(u^2 + v^2)
    
    double E =  W[0]*(k + e);

    F[0] = W[0]*W[1];                   // F[0] = rho*u
    F[1] = W[0]*W[1]*W[1] + W[3];       // F[1] = rho*u^2 + p
    F[2] = W[0]*W[1]*W[2];              // F[2] = rho*u*v
    F[3] = W[1]*(E + W[3]);             // F[3] = u*(E + p)

    return F;
}

Vector<double> flux_conserved(Vector<double> U) {

    Vector<double> F(4);
    Vector<double> W(4); 
    
    W = conserved_to_primitive(U); 

    double e = W[3]/((GAMMA - 1.0)*W[0]);   // Internal energy, e = p/((gamma-1)*rho)
    double k = 0.5*(W[1]*W[1] + W[2]*W[2]); // Kinetic energy, k = (1/2)*(u^2 + v^2)
    
    double E =  W[0]*(k + e);

    F[0] = W[0]*W[1];                   // F[0] = rho*u
    F[1] = W[0]*W[1]*W[1] + W[3];       // F[1] = rho*u^2 + p
    F[2] = W[0]*W[1]*W[2];              // F[2] = rho*u*v
    F[3] = W[1]*(E + W[3]);             // F[3] = u*(E + p)

    return F;
}

double speed_of_sound_primitive_variable(Vector<double> W) {
    
    return std::sqrt(GAMMA*W[3]/W[0]);
}

double speed_of_sound_conserved_variable(Vector<double> U) {

    double p = (GAMMA-1.0)*(U[3] - 0.5*((U[1]*U[1] + U[2]*U[2])/U[0]));
    
    return std::sqrt(GAMMA*p/U[0]);
}

Vector<double> multiply_with_rotation_matrix(Vector<double> U, double nx, double ny) {
    
    Vector<double> U_hat(4); 
    
    U_hat[0] = U[0]; 
    U_hat[1] = U[1]*nx + U[2]*ny; 
    U_hat[2] = -U[1]*ny + U[2]*nx;
    U_hat[3] = U[3]; 
    
    return U_hat; 
}

Vector<double> multiply_with_inv_rotation_matrix(Vector<double> U, double nx, double ny) {
    
    Vector<double> U_hat(4); 
    
    U_hat[0] = U[0]; 
    U_hat[1] = U[1]*nx - U[2]*ny; 
    U_hat[2] = U[1]*ny + U[2]*nx;
    U_hat[3] = U[3]; 
    
    return U_hat; 
}

// ============================================================================

// Riemann Solvers 

// ============================================================================

// 1. Exact Riemann Solver 

// Functions to find exact pstar and ustar

double fK(double p, Vector<double> WK) {

    double rhoK = WK[0];
    double pK = WK[3];
    double retval; 

    double A_K = (2.0)/((GAMMA + 1.0)*rhoK);
    double B_K = ((GAMMA - 1)/(GAMMA + 1))*pK;
    double a_K = speed_of_sound_primitive_variable(WK);

    if (p > pK) {                                  // Shock
        retval = (p-pK)*std::sqrt(A_K/(p + B_K));
    }
    else {                                         // Rarefraction
        retval = (2*a_K/(GAMMA-1))*(std::pow((p/pK),((GAMMA-1)/(2*GAMMA))) - 1.0);
    }
    
    return retval;
}

double fK_deriv(double p, Vector<double> WK) {

    double rhoK = WK[0];
    double pK = WK[3];
    double retval; 

    double A_K = (2.0)/((GAMMA + 1)*rhoK);
    double B_K = ((GAMMA - 1)/(GAMMA + 1))*pK;
    double a_K = speed_of_sound_primitive_variable(WK);

    if (p > pK) {                                  // Shock
        retval = std::sqrt((A_K/(B_K + p)))*(1 - ((p - pK)/(2*(B_K + p)))); 
    }

    else {                                         // Rarefraction
        retval = (1.0/(rhoK*a_K))*((std::pow((p/pK), ((-1.0-GAMMA)/(2*GAMMA))))); 
    }

    return retval;
}

double star_pu_function(double p, Vector<double> WL, Vector<double> WR) {

    double u_L = WL[1];
    double u_R = WR[1];

    double fL = fK(p, WL);
    double fR = fK(p, WR);

    double f = fL + fR + (u_R - u_L);

    return f;
}

double star_pu_function_deriv(double p, Vector<double> WL, Vector<double> WR) {


    double fL_deriv = fK_deriv(p, WL);
    double fR_deriv = fK_deriv(p, WR);

    double f_deriv = fL_deriv + fR_deriv;

    return f_deriv;
}

double solve_for_pstar_newton_method(double(*f)(double, Vector<double>, Vector<double>), double(*f_deriv)(double, Vector<double>, Vector<double>), 
                                     Vector<double> WL, Vector<double> WR, double p_guess) {
    
    double p_star = p_guess; 
    //double tol=1.48e-08; 
    unsigned int maxiter = 50; 
    unsigned int iter = 0; 
    
    while (iter < maxiter) {
    
        p_star = p_guess - f(p_guess, WL, WR)/f_deriv(p_guess, WL, WR);
        p_guess = p_star; 
        
        iter ++; 
    }
    
    return p_star; 
    
}

void get_pstar_and_ustar(Vector<double> WL, Vector<double> WR, double& pstar, double& ustar) {

    double u_L = WL[1]; 
    double p_L = WL[3]; 
    double a_L = speed_of_sound_primitive_variable(WL); 

    double u_R = WR[1];
    double p_R = WR[3];
    double a_R = speed_of_sound_primitive_variable(WR);
   
    double g1 = (GAMMA - 1.0)/(2.0*GAMMA);

    // Two–Rarefaction approximation

    double term1 = a_L + a_R - 0.5*(GAMMA-1)*(u_R - u_L);
    double term2 = a_L/(std::pow(p_L, g1)) + a_R/(std::pow(p_R, g1));
    double p0 = std::pow((term1/term2),(1./g1));
    
    pstar = solve_for_pstar_newton_method(star_pu_function, star_pu_function_deriv, WL, WR, p0);
    ustar = 0.5*(u_L + u_R) + 0.5*(fK(pstar, WR) - fK(pstar, WL));
    
}

Vector<double> sample(double p_star, double u_star, Vector<double> WL, Vector<double> WR) {

    Vector<double> W(4);
    
    double S = 0.0; 

    // GAMMA and related constants

    double g1 = (GAMMA - 1.0)/(2.0*GAMMA);
    double g2 = (GAMMA + 1.0)/(2.0*GAMMA);
    double g3 = 2.0*GAMMA/(GAMMA - 1.0);
    double g4 = 2.0/(GAMMA - 1.0);
    double g5 = 2.0/(GAMMA + 1.0);
    double g6=(GAMMA-1.0)/(GAMMA + 1.0);
    double g7=(GAMMA - 1.0)/2.0;
	// g8 = GAMMA - 1.0;

    // Simplify notation - declare few extra variables

    double rho_L = WL[0]; 
    double u_L = WL[1];
    double v_L = WL[2];
    double p_L = WL[3];

    double rho_R = WR[0];
    double u_R = WR[1];
    double v_R = WR[2]; 
    double p_R = WR[3]; 

    // Acoustic Speeds

    double a_L = std::sqrt((p_L*GAMMA)/(rho_L));
    double a_star_L = a_L*std::pow((p_star/p_L), g1);
    double a_R = std::sqrt((p_R*GAMMA)/(rho_R));
    double a_star_R = a_R*std::pow((p_star/p_R), g1);


	// Shock and rarefaction  speeds

    double S_L = u_L - a_L*(std::sqrt(g2*(p_star/p_L) + g1));         // Left shock
    double S_HL = u_L - a_L;                                          // Left rarefaction head
    double S_TL = u_star - a_star_L;                                  // Left rarefaction tail

    double S_R = u_R + a_R*(std::sqrt(g2*(p_star/p_R) + g1));         // Right shock
    double S_HR = u_R + a_R;                                          // Right rarefaction head
    double S_TR = u_star + a_star_R;                                  // Right rarefaction tail


    if (S < u_star) {     // Left side of contact
        
        if (p_star > p_L) {   // Left side shock
            
            if(S_L <= S) {
                W[0] = rho_L*(((p_star/p_L)+g6)/(g6*(p_star/p_L) + 1));
                W[1] = u_star; 
                W[3] = p_star;
            }

            else if(S < S_L) {
                W[0] = rho_L;
                W[1] = u_L;
                W[3] = p_L;
            }
            
        }

        else if (p_L >= p_star) {  // Left side rarefaction
            
            if (S <= S_HL) {
                W[0] = rho_L;
                W[1] = u_L;
                W[3] = p_L;
            }

            else if (S_HL < S && S <= S_TL) {
                W[0] = rho_L*std::pow((g5 + (g6/a_L)*(u_L - S)), g4);
                W[1] = g5*(a_L + g7*u_L + S);
                W[3] = p_L*(std::pow((g5 + (g6/a_L)*(u_L - S)),g3));
            }

            else if (S_TL < S && S <= u_star) {
                W[0] = rho_L*std::pow((p_star/p_L),(1.0/GAMMA));
                W[1] = u_star;
                W[3] = p_star;
            }
        }
    }

    else if (S >= u_star) {  // Right side of contact
        
        if (p_star > p_R) {   // Right side shock
            
            if (u_star <= S && S< S_R) {
                W[0] = rho_R*(((p_star/p_R)+g6)/(g6*(p_star/p_R) + 1));
                W[1] = u_star;
                W[3] = p_star;
            }

            else if(S >= S_R) { 
                W[0] = rho_R;
                W[1] = u_R;
                W[3] = p_R;
            }
        }
        
        else if (p_R >= p_star) {  // Right side rarefaction
            
            if (S >= S_HR) {
                W[0] = rho_R;
                W[1] = u_R;
                W[3] = p_R;
            }

            else if (S_TR <= S && S <= S_HR) {
                W[0] = rho_R*std::pow((g5 - (g6/a_R)*(u_R - S)),g4);
                W[1] = g5*(-a_R + g7*u_R + S);
                W[3] = p_R*std::pow((g5 - (g6/a_R)*(u_R - S)),g3);
            }

            else if (u_star <= S && S<= S_TR) {
                W[0] = rho_R*std::pow((p_star/p_R), (1.0/GAMMA));
                W[1] = u_star;
                W[3] = p_star;
            }
        }
    }
    
    if (0.0 < u_star) {
        W[2] = v_L; 
    }
    
    else {
        W[2] = v_R; 
    }

    return W; 
}


Vector<double> exact_riemann_solver(Vector<double> WL, Vector<double> WR, double nx, double ny, Point<2> P, bool boundary) {
    
    if (WL[0] < 0.0 || WR[0] < 0.0 || WL[3] < 0.0 || WR[3] < 0.0) {
        ThrowNegativePressureDensityError(WL, WR, P, boundary); 
    }

    Vector <double> W(4); 
    double pstar, ustar; 
    
    Vector<double> W_hatL(4); Vector<double> W_hatR(4); Vector<double> W_hat(4);  
    
    W_hatL = multiply_with_rotation_matrix(WL, nx, ny);
    W_hatR = multiply_with_rotation_matrix(WR, nx, ny);
    
    get_pstar_and_ustar(W_hatL, W_hatR, pstar, ustar); 

    W_hat = sample(pstar, ustar, W_hatL, W_hatR);
    
    Vector<double> F_hat(4), F(4); 
    
    F_hat = flux_primitive(W_hat); 
    F = multiply_with_inv_rotation_matrix(F_hat, nx, ny);
    
    return F; 
}

// ============================================================================

// 2. HLLC Riemann Solver 

Vector<double> HLLC_riemann_solver(Vector<double> WL, Vector<double> WR, double nx, double ny, Point<2> P, bool boundary) {

    if (WL[0] < 0.0 || WR[0] < 0.0 || WL[3] < 0.0 || WR[3] < 0.0) {
        ThrowNegativePressureDensityError(WL, WR, P, boundary); 
    }
    
    Vector<double> F(4), F_hat(4), W_hatL(4), W_hatR(4);

    W_hatL = multiply_with_rotation_matrix(WL, nx, ny);
    W_hatR = multiply_with_rotation_matrix(WR, nx, ny);
	
    //double g3 = (GAMMA + 1.0)/(2*GAMMA);

	// Extract states
	double rho_L = W_hatL[0]; double u_L = W_hatL[1]; double v_L = W_hatL[2]; double p_L = W_hatL[3];
	double rho_R = W_hatR[0]; double u_R = W_hatR[1]; double v_R = W_hatR[2]; double p_R = W_hatR[3];

	double a_L = std::sqrt(GAMMA*p_L/rho_L); double a_R = std::sqrt(GAMMA*p_R/rho_R);
    /*
	// Pressure Estimates

	double rho_bar  = 0.5*(rho_L + rho_R); double a_bar = 0.5*(a_L + a_R);
	double p_pvrs = 0.5*(p_L + p_R) - 0.5*(u_R - u_L)*rho_bar*a_bar;
	double p_star = std::max(0.0, p_pvrs);

	// Wave speed estimates

	double q_L, q_R;

	if (p_star <= p_L) {
		q_L = 1.0;
	}
	else {
		q_L = std::sqrt(1.0 + g3*((p_star/p_L) - 1.0));
	}

	if (p_star <= p_R) {
		q_R = 1.0;
	}
	else {
		q_R = std::sqrt(1.0 + g3*((p_star/p_R) - 1.0));
	}

	double S_L = u_L - a_L*q_L; double S_R = u_R + a_R*q_R;

	double S_star = (p_R - p_L + rho_L*u_L*(S_L - u_L) - rho_R*u_R*(S_R - u_R))
			/(rho_L*(S_L - u_L) - rho_R*(S_R- u_R));
    */ 
    double S_L = std::min((u_R - a_R), (u_L - a_L)) ; 
	double S_R = std::max((u_L + a_L), (u_R + a_R)) ;
	double S_star = ( p_R - p_L + rho_L*u_L*(S_L-u_L) - rho_R*u_R*(S_R - u_R) )/(rho_L*(S_L - u_L) - rho_R*(S_R - u_R)) ;

	// HLLC Flux

    if (0.0 <= S_L) {
        F_hat = flux_primitive(W_hatL); 
    }
    
    else if (S_L <= 0.0 && 0.0 <= S_star) {
    
        Vector<double> U_hatL(4), U_starL(4), F_hatL(4); 
        
        U_hatL = primitive_to_conserved(W_hatL);
        F_hatL = flux_primitive(W_hatL); 
        
        U_starL[0] = rho_L*((S_L - u_L)/(S_L - S_star));
        U_starL[1] = rho_L*((S_L - u_L)/(S_L - S_star))*S_star;
        U_starL[2] = rho_L*((S_L - u_L)/(S_L - S_star))*v_L; 
        U_starL[3] = rho_L*((S_L - u_L)/(S_L - S_star))*(U_hatL[3]/rho_L + (S_star - u_L)*(S_star + p_L/(rho_L*(S_L - u_L))));
        
        F_hat[0] = F_hatL[0] + S_L*(U_starL[0] - U_hatL[0]);    
        F_hat[1] = F_hatL[1] + S_L*(U_starL[1] - U_hatL[1]);
        F_hat[2] = F_hatL[2] + S_L*(U_starL[2] - U_hatL[2]);
        F_hat[3] = F_hatL[3] + S_L*(U_starL[3] - U_hatL[3]);
        
    }
    
    else if (S_star <= 0.0 && 0.0 <= S_R) {
        
        Vector<double> U_hatR(4), U_starR(4), F_hatR(4); 
        
        U_hatR = primitive_to_conserved(W_hatR);
        F_hatR = flux_primitive(W_hatR); 
        
        U_starR[0] = rho_R*((S_R - u_R)/(S_R - S_star));
        U_starR[1] = rho_R*((S_R - u_R)/(S_R - S_star))*S_star;
        U_starR[2] = rho_R*((S_R - u_R)/(S_R - S_star))*v_R; 
        U_starR[3] = rho_R*((S_R - u_R)/(S_R - S_star))*(U_hatR[3]/rho_R + (S_star - u_R)*(S_star + p_R/(rho_R*(S_R - u_R))));
        
        F_hat[0] = F_hatR[0] + S_R*(U_starR[0] - U_hatR[0]);    
        F_hat[1] = F_hatR[1] + S_R*(U_starR[1] - U_hatR[1]);
        F_hat[2] = F_hatR[2] + S_R*(U_starR[2] - U_hatR[2]);
        F_hat[3] = F_hatR[3] + S_R*(U_starR[3] - U_hatR[3]);
    
    }
    
    else if (0.0 >= S_R) {
        F_hat = flux_primitive(W_hatR); 
    }

    F = multiply_with_inv_rotation_matrix(F_hat, nx, ny);
	
    return F;
}

// ============================================================================

// 3. Roe Riemann Solver 

Vector<double> Roe_riemann_solver(Vector<double> WL, Vector<double> WR, double nx, double ny, Point<2> P, bool boundary) {
    
    if (WL[0] < 0.0 || WR[0] < 0.0 || WL[3] < 0.0 || WR[3] < 0.0) {
        ThrowNegativePressureDensityError(WL, WR, P, boundary); 
    }
    
    //WARNING: No entropy fix. May misbehave in presence of transonic/sonic rarefactions
    
    Vector<double> F(4), F_hat(4), W_hatL(4), W_hatR(4);

    W_hatL = multiply_with_rotation_matrix(WL, nx, ny);
    W_hatR = multiply_with_rotation_matrix(WR, nx, ny);

    // Extract states
	double rho_L = W_hatL[0]; double u_L = W_hatL[1]; double v_L = W_hatL[2]; double p_L = W_hatL[3];
	double rho_R = W_hatR[0]; double u_R = W_hatR[1]; double v_R = W_hatR[2]; double p_R = W_hatR[3];

    double H_L = (0.5*u_L*u_L + p_L/((GAMMA - 1)*rho_L)) + p_L/rho_L; 
    double H_R = (0.5*u_R*u_R + p_R/((GAMMA - 1)*rho_R)) + p_R/rho_R;

    // Roe averages
    double rho_bar = std::sqrt(rho_L*rho_R);
    double u_bar = (std::sqrt(rho_L)*u_L + std::sqrt(rho_R)*u_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));
    double v_bar = (std::sqrt(rho_L)*v_L + std::sqrt(rho_R)*v_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));     
    double H_bar = (std::sqrt(rho_L)*H_L + std::sqrt(rho_R)*H_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));
    double a_bar = std::sqrt((GAMMA - 1.0)*(H_bar - 0.5*(u_bar*u_bar + v_bar*v_bar)));

    // Compute the eigenvalues
    double lambda[] = {u_bar - a_bar, u_bar, u_bar, u_bar + a_bar}; 

    // Compute the right eigenvectors 
    Vector<double> K1(4); Vector<double> K2(4); Vector<double> K3(4); Vector<double> K4(4);

    K1[0] = 1.0; K1[1] = u_bar - a_bar; K1[2] = v_bar; K1[3] = H_bar - u_bar*a_bar; 
    K2[0] = 1.0; K2[1] = u_bar; K2[2] = v_bar; K2[3] = 0.5*(u_bar*u_bar + v_bar*v_bar); 
    K3[0] = 0.0; K3[1] = 0.0; K3[2] = 1.0; K3[3] = v_bar; 
    K4[0] = 1.0; K4[1] = u_bar + a_bar; K4[2] = v_bar; K4[3] = H_bar + u_bar*a_bar; 

    // Compute the wave strengths
    double alpha[4]; 

    alpha[0] = (0.5/a_bar*a_bar)*(p_R - p_L - rho_bar*a_bar*(u_R - u_L));
    alpha[1] = (rho_R - rho_L) - (1./(a_bar*a_bar))*(p_R - p_L);
    alpha[2] = rho_bar*(v_R - v_L);
    alpha[3] = (0.5/a_bar*a_bar)*(p_R - p_L + rho_bar*a_bar*(u_R - u_L));

    // Compute the intecell flux 
    Vector<double> F_hatL(4); Vector<double> F_hatR(4); 

    F_hatL = flux_primitive(W_hatL); F_hatR = flux_primitive(W_hatR);

    F_hat[0] = 0.5*(F_hatL[0] + F_hatR[0]) - 
               0.5*(alpha[0]*std::abs(lambda[0])*K1[0] + alpha[1]*std::abs(lambda[1])*K2[0] +
                    alpha[2]*std::abs(lambda[2])*K3[0] + alpha[3]*std::abs(lambda[3])*K4[0]); 
    
    F_hat[1] = 0.5*(F_hatL[1] + F_hatR[1]) - 
               0.5*(alpha[0]*std::abs(lambda[0])*K1[1] + alpha[1]*std::abs(lambda[1])*K2[1] +
                    alpha[2]*std::abs(lambda[2])*K3[1] + alpha[3]*std::abs(lambda[3])*K4[1]); 
    
    F_hat[2] = 0.5*(F_hatL[2] + F_hatR[2]) - 
               0.5*(alpha[0]*std::abs(lambda[0])*K1[2] + alpha[1]*std::abs(lambda[1])*K2[2] +
                    alpha[2]*std::abs(lambda[2])*K3[2] + alpha[3]*std::abs(lambda[3])*K4[2]); 

    F_hat[3] = 0.5*(F_hatL[3] + F_hatR[3]) - 
               0.5*(alpha[0]*std::abs(lambda[0])*K1[3] + alpha[1]*std::abs(lambda[1])*K2[3] +
                    alpha[2]*std::abs(lambda[2])*K3[3] + alpha[3]*std::abs(lambda[3])*K4[3]); 

    // Multiply with inverse rotation matrix 

    F = multiply_with_inv_rotation_matrix(F_hat, nx, ny);

    return F; 
}

// ============================================================================

// 4. Local Lax Friedrichs Riemann Solver 

Vector<double> local_Lax_Friedrichs_riemann_solver(Vector<double> WL, Vector<double> WR, double nx, double ny, Point<2> P, bool boundary) {
    
    if (WL[0] < 0.0 || WR[0] < 0.0 || WL[3] < 0.0 || WR[3] < 0.0) {
        ThrowNegativePressureDensityError(WL, WR, P, boundary); 
    }
    
    Vector<double> F(4), F_hat(4), W_hatL(4), W_hatR(4);

    W_hatL = multiply_with_rotation_matrix(WL, nx, ny);
    W_hatR = multiply_with_rotation_matrix(WR, nx, ny);

    // Extract states
	double rho_L = W_hatL[0]; double u_L = W_hatL[1]; double v_L = W_hatL[2]; double p_L = W_hatL[3];
	double rho_R = W_hatR[0]; double u_R = W_hatR[1]; double v_R = W_hatR[2]; double p_R = W_hatR[3];

    double E_L = rho_L*(0.5*(u_L*u_L + v_L*v_L) + p_L/((GAMMA - 1.0)*rho_L));   
    double E_R = rho_R*(0.5*(u_R*u_R + v_R*v_R) + p_R/((GAMMA - 1.0)*rho_R));

    double a_L = speed_of_sound_primitive_variable(W_hatL); 
    double a_R = speed_of_sound_primitive_variable(W_hatR); 

    // Compute dissipation 
    double max_vel_L = a_L + std::abs(u_L); double max_vel_R = a_R + std::abs(u_R);
    double alpha = std::max(max_vel_L, max_vel_R);  

    // Compute flux 
    Vector<double> F_hatL(4); Vector<double> F_hatR(4); 
    F_hatL = flux_primitive(W_hatL); F_hatR = flux_primitive(W_hatR); 

    F_hat[0] = 0.5*(F_hatL[0] + F_hatR[0])  - 0.5*alpha*(rho_R - rho_L);
    F_hat[1] = 0.5*(F_hatL[1] + F_hatR[1])  - 0.5*alpha*(rho_R*u_R - rho_L*u_L);
    F_hat[2] = 0.5*(F_hatL[2] + F_hatR[2])  - 0.5*alpha*(rho_R*v_R - rho_L*v_L);
    F_hat[3] = 0.5*(F_hatL[3] + F_hatR[3])  - 0.5*alpha*(E_R - E_L);

    // Multiply with inverse rotation matrix 
    F = multiply_with_inv_rotation_matrix(F_hat, nx, ny);

    return F; 
}

// ============================================================================

// 5. HLL Riemann Solver 

Vector<double> HLL_riemann_solver(Vector<double> WL, Vector<double> WR, double nx, double ny, Point<2> P, bool boundary) {

    if (WL[0] < 0.0 || WR[0] < 0.0 || WL[3] < 0.0 || WR[3] < 0.0) {
        ThrowNegativePressureDensityError(WL, WR, P, boundary); 
    }

    // WARNING: Resolution of contact surfaces can be very inaccurate
    
    Vector<double> F(4), F_hat(4), W_hatL(4), W_hatR(4);

    W_hatL = multiply_with_rotation_matrix(WL, nx, ny);
    W_hatR = multiply_with_rotation_matrix(WR, nx, ny);

    // Extract states
	double rho_L = W_hatL[0]; double u_L = W_hatL[1]; double v_L = W_hatL[2]; double p_L = W_hatL[3];
	double rho_R = W_hatR[0]; double u_R = W_hatR[1]; double v_R = W_hatR[2]; double p_R = W_hatR[3];

    double H_L = (0.5*u_L*u_L + p_L/((GAMMA - 1)*rho_L)) + p_L/rho_L; 
    double H_R = (0.5*u_R*u_R + p_R/((GAMMA - 1)*rho_R)) + p_R/rho_R;

    double a_L = speed_of_sound_primitive_variable(W_hatL); 
    double a_R = speed_of_sound_primitive_variable(W_hatR);

    // Roe averages
    double u_bar = (std::sqrt(rho_L)*u_L + std::sqrt(rho_R)*u_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));
    double v_bar = (std::sqrt(rho_L)*v_L + std::sqrt(rho_R)*v_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));     
    double H_bar = (std::sqrt(rho_L)*H_L + std::sqrt(rho_R)*H_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));
    double a_bar = std::sqrt((GAMMA - 1.0)*(H_bar - 0.5*(u_bar*u_bar + v_bar*v_bar)));

    // Wavespeeds 
    double S_L = std::min(a_L, u_bar - a_bar); double S_R = std::min(a_R, u_bar + a_bar); 

    if (S_L > 0.0) {
        F_hat = flux_primitive(W_hatL);     
    }

    else if (S_R < 0.0) {
        F_hat = flux_primitive(W_hatR); 
    }

    else {
        Vector<double> U_hatL(4); Vector<double> U_hatR(4);
        Vector<double> F_hatL(4); Vector<double> F_hatR(4);

        U_hatL = primitive_to_conserved(W_hatL); U_hatR = primitive_to_conserved(W_hatR);
        F_hatL = flux_primitive(W_hatL); F_hatR = flux_primitive(W_hatR); 

        F_hat[0] = (S_R/(S_R - S_L))*F_hatL[0] - (S_L/(S_R - S_L))*F_hatR[0] +((S_R*S_L)/(S_R - S_L))*(U_hatR[0] - U_hatL[0]);
        F_hat[1] = (S_R/(S_R - S_L))*F_hatL[1] - (S_L/(S_R - S_L))*F_hatR[1] +((S_R*S_L)/(S_R - S_L))*(U_hatR[1] - U_hatL[1]);
        F_hat[2] = (S_R/(S_R - S_L))*F_hatL[2] - (S_L/(S_R - S_L))*F_hatR[2] +((S_R*S_L)/(S_R - S_L))*(U_hatR[2] - U_hatL[2]);
        F_hat[3] = (S_R/(S_R - S_L))*F_hatL[3] - (S_L/(S_R - S_L))*F_hatR[3] +((S_R*S_L)/(S_R - S_L))*(U_hatR[3] - U_hatL[3]);
    }

    // Multiply with inverse rotation matrix 
    F = multiply_with_inv_rotation_matrix(F_hat, nx, ny);

    return F; 
}

// ============================================================================

// 6. HLLI Riemann Solver 

Vector<double> HLLI_riemann_solver(Vector<double> WL, Vector<double> WR, double nx, double ny, Point<2> P, bool boundary) {

    if (WL[0] < 0.0 || WR[0] < 0.0 || WL[3] < 0.0 || WR[3] < 0.0) {
        ThrowNegativePressureDensityError(WL, WR, P, boundary); 
    }
    
    Vector<double> F(4), F_hat(4), W_hatL(4), W_hatR(4);

    W_hatL = multiply_with_rotation_matrix(WL, nx, ny);
    W_hatR = multiply_with_rotation_matrix(WR, nx, ny);

    // Extract states
	double rho_L = W_hatL[0]; double u_L = W_hatL[1]; double v_L = W_hatL[2]; double p_L = W_hatL[3];
	double rho_R = W_hatR[0]; double u_R = W_hatR[1]; double v_R = W_hatR[2]; double p_R = W_hatR[3];

    double H_L = (0.5*u_L*u_L + p_L/((GAMMA - 1)*rho_L)) + p_L/rho_L; 
    double H_R = (0.5*u_R*u_R + p_R/((GAMMA - 1)*rho_R)) + p_R/rho_R;

    double a_L = speed_of_sound_primitive_variable(W_hatL); 
    double a_R = speed_of_sound_primitive_variable(W_hatR);

    // Roe averages
    double u_bar = (std::sqrt(rho_L)*u_L + std::sqrt(rho_R)*u_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));
    double v_bar = (std::sqrt(rho_L)*v_L + std::sqrt(rho_R)*v_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));     
    double H_bar = (std::sqrt(rho_L)*H_L + std::sqrt(rho_R)*H_R)/(std::sqrt(rho_L) + std::sqrt(rho_R));
    double a_bar = std::sqrt((GAMMA - 1.0)*(H_bar - 0.5*(u_bar*u_bar + v_bar*v_bar)));

    // Wavespeeds 
    double S_L = std::min(a_L, u_bar - a_bar); double S_R = std::min(a_R, u_bar + a_bar); 

    if (S_L > 0.0) {
        F_hat = flux_primitive(W_hatL);     
    }

    else if (S_R < 0.0) {
        F_hat = flux_primitive(W_hatR); 
    }

    else { 
        // Compute HLL Flux 
        
        Vector<double> U_hatL(4); Vector<double> U_hatR(4);
        Vector<double> F_hatL(4); Vector<double> F_hatR(4);

        U_hatL = primitive_to_conserved(W_hatL); U_hatR = primitive_to_conserved(W_hatR);
        F_hatL = flux_primitive(W_hatL); F_hatR = flux_primitive(W_hatR); 

        F_hat[0] = (S_R/(S_R - S_L))*F_hatL[0] - (S_L/(S_R - S_L))*F_hatR[0] +((S_R*S_L)/(S_R - S_L))*(U_hatR[0] - U_hatL[0]);
        F_hat[1] = (S_R/(S_R - S_L))*F_hatL[1] - (S_L/(S_R - S_L))*F_hatR[1] +((S_R*S_L)/(S_R - S_L))*(U_hatR[1] - U_hatL[1]);
        F_hat[2] = (S_R/(S_R - S_L))*F_hatL[2] - (S_L/(S_R - S_L))*F_hatR[2] +((S_R*S_L)/(S_R - S_L))*(U_hatR[2] - U_hatL[2]);
        F_hat[3] = (S_R/(S_R - S_L))*F_hatL[3] - (S_L/(S_R - S_L))*F_hatR[3] +((S_R*S_L)/(S_R - S_L))*(U_hatR[3] - U_hatL[3]);

        // Compute Left and Right Eigenvectors 

        Vector<double> U_star(4); 

        U_star[0] = ((S_R*U_hatR[0]) - (S_L*U_hatL[0]) - (F_hatR[0] - F_hatL[0]))/(S_R - S_L);
        U_star[1] = ((S_R*U_hatR[1]) - (S_L*U_hatL[1]) - (F_hatR[1] - F_hatL[1]))/(S_R - S_L);
        U_star[2] = ((S_R*U_hatR[2]) - (S_L*U_hatL[2]) - (F_hatR[2] - F_hatL[2]))/(S_R - S_L);
        U_star[3] = ((S_R*U_hatR[3]) - (S_L*U_hatL[3]) - (F_hatR[3] - F_hatL[3]))/(S_R - S_L);

        double a_star = speed_of_sound_conserved_variable(U_star); 

        Vector<double> L_Eig_Vec(4); Vector<double> R_Eig_Vec(4);

        R_Eig_Vec[0] = 1.0; 
        R_Eig_Vec[1] = U_star[1]/U_star[0]; 
        R_Eig_Vec[2] = U_star[2]/U_star[0];
        R_Eig_Vec[3] = 0.5*(R_Eig_Vec[1]*R_Eig_Vec[1] + R_Eig_Vec[2]*R_Eig_Vec[2]);

        L_Eig_Vec[0] = (a_star*a_star - (GAMMA-1.0)*R_Eig_Vec[3])/(a_star*a_star); 
        L_Eig_Vec[1] = ((GAMMA-1.0)*R_Eig_Vec[1])/(a_star*a_star);; 
        L_Eig_Vec[2] = ((GAMMA-1.0)*R_Eig_Vec[2])/(a_star*a_star);
        L_Eig_Vec[3] = (1.0-GAMMA)/(a_star*a_star);

        double l_norm = L_Eig_Vec.l2_norm(); double r_norm = R_Eig_Vec.l2_norm(); 

        for (unsigned int i = 0; i < 4; i++) {
            L_Eig_Vec[i] = L_Eig_Vec[i]/l_norm; 
            R_Eig_Vec[i] = R_Eig_Vec[i]/r_norm;
        }

        double delta = 1.0 - std::min(U_star[1]/U_star[0], 0.0)/S_L - std::max(U_star[1]/U_star[0], 0.0)/S_R; 

        U_hatR.sadd(1.0, -1.0, U_hatL);

        // Form dot product 

        double dot_prod = 0.0;  

        for (unsigned int i = 0; i < 4; i++) {
            dot_prod += L_Eig_Vec[i]*U_hatR[i]; 
        }

        // Correct the fluxes 
        F_hat[0] = F_hat[0] - ((S_R*S_L)/(S_R-S_L))*delta*dot_prod*R_Eig_Vec[0]; 
        F_hat[1] = F_hat[1] - ((S_R*S_L)/(S_R-S_L))*delta*dot_prod*R_Eig_Vec[1];
        F_hat[2] = F_hat[2] - ((S_R*S_L)/(S_R-S_L))*delta*dot_prod*R_Eig_Vec[2];
        F_hat[3] = F_hat[3] - ((S_R*S_L)/(S_R-S_L))*delta*dot_prod*R_Eig_Vec[3];

    }

    // Multiply with inverse rotation matrix 
    F = multiply_with_inv_rotation_matrix(F_hat, nx, ny);

    return F; 
}

// 7. Rotated HLLC Riemann Solver 

Vector<double> rotated_HLLC_riemann_solver(Vector<double> WL, Vector<double> WR, double nx, double ny, Point<2> P, bool boundary) {
    
    if (WL[0] < 0.0 || WR[0] < 0.0 || WL[3] < 0.0 || WR[3] < 0.0) {
        ThrowNegativePressureDensityError(WL, WR, P, boundary); 
    }
    
    double alpha1, alpha2, n1x, n1y, n2x, n2y; 
    
	double du = WR(1) - WL(1); 
    double dv = WR(2) - WL(2); 
	double dq = std::sqrt(du*du + dv*dv) ;
    
    if(dq < 1.0e-10) { 
        n1x = nx; n1y = ny; 
    }
    
	else { 
        n1x = du/dq; n1y = dv/dq; 
    }
    
    alpha1 = (n1x*nx + n1y*ny) ;
	
    if(alpha1 < 0) { 
        n1x = -n1x; n1y = -n1y; alpha1 = -alpha1; 
    }
    
	n2x = -n1y ; n2y = n1x ;
	alpha2 = (n2x*nx + n2y*ny) ; 
	
    if(alpha2 < 0) { 
        n2x = -n2x; n2y = -n2y; alpha2 = -alpha2; 
    }
    
    Vector<double> F1(4); Vector<double> F2(4); Vector<double> F(4);; 
    
    F1 = HLLC_riemann_solver(WL, WR, n1x, n1y, P,  boundary);
    F2 = HLLC_riemann_solver(WL, WR, n2x, n2y, P,  boundary);
    
    for (unsigned int i = 0; i < 4; i++) {
        F(i) = alpha1*F1(i) + alpha2*F2(i); 
    }

    return F; 
} 

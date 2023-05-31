/*
 * weno.cc
 *      Author: sunder
 */

#include "weno.hh"

double sgn(double a) {
	if (a>0.0) return 1.0;
	else if (a==0.0) return 0.0;
	else return -1.0; 
}

double minmod(double a, double b) {

	return 0.5*(sgn(a)+ sgn(b))*std::min(std::abs(a), std::abs(b)); 
	
	//return 0.5*(a+b);
	
}

// Second order WENO (linear portion of r=3 centered WENO limiting)

void get_weno2_coefficients(const std::vector<double>& U_x,
		   const std::vector<double>& U_y,
		   const std::vector<double>& U_xy,
		   std::vector<double>& coeffs) {
	/*
	static const double epsilon = 1.0e-12;
    static const double central_cell_wt = 100.0;
    int p = 4;
    double u_x, u_y;

    //----------------------- get u_x and u_xx -----------------------

    // Left stencil
    double uL_x = -2.0*U_x[1] + 0.5*U_x[0] + 1.5*U_x[2];
	double uL_xx = 0.5*U_x[0] -     U_x[1] + 0.5*U_x[2];

	// Central Stencil
	double uC_x  = 0.5*U_x[3] -          0.5*U_x[1];
	double uC_xx = 0.5*U_x[1] - U_x[2] + 0.5*U_x[3];

    // Right stencil
	double uR_x  = -1.5*U_x[2] + 2.0*U_x[3] - 0.5*U_x[4];
	double uR_xx =  0.5*U_x[2] -     U_x[3] + 0.5*U_x[4];

    // Smoothness Indicators
	double IS_L = uL_x*uL_x + 4.333333333333333*(uL_xx*uL_xx); // Left stencil
	double IS_C = uC_x*uC_x + 4.333333333333333*(uC_xx*uC_xx); // Central stencil
	double IS_R = uR_x*uR_x + 4.333333333333333*(uR_xx*uR_xx); // Right stencil

    double wL = 1.0/(std::pow((IS_L + epsilon), p));
	double wC = central_cell_wt/(std::pow((IS_C + epsilon), p));
	double wR = 1.0/(std::pow((IS_R + epsilon), p));

	double sum = wL + wC + wR;

	// Normalize the weights

	wL = wL/sum; wC = wC/sum; wR = wR/sum;

    u_x = wL*uL_x + wC*uC_x + wR*uR_x;

	//----------------------- get u_y and u_yy -----------------------

    // Bottom stencil
    double uB_y = -2.0*U_y[1] + 0.5*U_y[0] + 1.5*U_y[2];
	double uB_yy = 0.5*U_y[0] -     U_y[1] + 0.5*U_y[2];

	// Middle Stencil
	double uM_y  = 0.5*U_y[3] -          0.5*U_y[1];
	double uM_yy = 0.5*U_y[1] - U_y[2] + 0.5*U_y[3];

    // Top stencil
	double uT_y  = -1.5*U_y[2] + 2.0*U_y[3] - 0.5*U_y[4];
	double uT_yy =  0.5*U_y[2] -     U_y[3] + 0.5*U_y[4];

    // Smoothness Indicators
	double IS_B = uB_y*uB_y + 4.333333333333333*(uB_yy*uB_yy); // Bottom stencil
	double IS_M = uM_y*uM_y + 4.333333333333333*(uM_yy*uM_yy); // Middle stencil
	double IS_T = uT_y*uT_y + 4.333333333333333*(uT_yy*uT_yy); // Top stencil

    double wB = 1.0/(std::pow((IS_B + epsilon), p));
	double wM = central_cell_wt/(std::pow((IS_M + epsilon), p));
	double wT = 1.0/(std::pow((IS_T + epsilon), p));

    sum = wB + wM + wT;

	// Normalize the weights

	wB = wB/sum; wM = wM/sum; wT = wT/sum;

    u_y = wB*uB_y + wM*uM_y + wT*uT_y;
	
	*/
	
	double u_x = minmod(U_x[3] - U_x[2], U_x[2] - U_x[1]);
	double u_y = minmod(U_y[3] - U_y[2], U_y[2] - U_y[1]);
	
    coeffs[0] = U_xy[0];
    coeffs[1] = u_x;
    coeffs[2] = u_y;
}


// Third order WENO

void get_weno3_coefficients(const std::vector<double>& U_x,
		   const std::vector<double>& U_y,
		   const std::vector<double>& U_xy,
		   std::vector<double>& coeffs) {


	static const double epsilon = 1.0e-12;
    double central_cell_wt = 100.0;
    int p = 4;
    double u_x, u_xx, u_y, u_yy, u_xy;

    //----------------------- get u_x and u_xx -----------------------

    // Left stencil
    double uL_x = -2.0*U_x[1] + 0.5*U_x[0] + 1.5*U_x[2];
	double uL_xx = 0.5*U_x[0] -     U_x[1] + 0.5*U_x[2];

	// Central Stencil
	double uC_x  = 0.5*U_x[3] -          0.5*U_x[1];
	double uC_xx = 0.5*U_x[1] - U_x[2] + 0.5*U_x[3];

    // Right stencil
	double uR_x  = -1.5*U_x[2] + 2.0*U_x[3] - 0.5*U_x[4];
	double uR_xx =  0.5*U_x[2] -     U_x[3] + 0.5*U_x[4];

    // Smoothness Indicators
	double IS_L = uL_x*uL_x + 4.333333333333333*(uL_xx*uL_xx); // Left stencil
	double IS_C = uC_x*uC_x + 4.333333333333333*(uC_xx*uC_xx); // Central stencil
	double IS_R = uR_x*uR_x + 4.333333333333333*(uR_xx*uR_xx); // Right stencil

    double wL = 1.0/(std::pow((IS_L + epsilon), p));
	double wC = central_cell_wt/(std::pow((IS_C + epsilon), p));
	double wR = 1.0/(std::pow((IS_R + epsilon), p));


    double sum = wL + wC + wR;

	// Normalize the weights

	wL = wL/sum; wC = wC/sum; wR = wR/sum;

    u_x = wL*uL_x + wC*uC_x + wR*uR_x;
	u_xx = wL*uL_xx + wC*uC_xx + wR*uR_xx;

	//----------------------- get u_y and u_yy -----------------------

    // Bottom stencil
    double uB_y = -2.0*U_y[1] + 0.5*U_y[0] + 1.5*U_y[2];
	double uB_yy = 0.5*U_y[0] -     U_y[1] + 0.5*U_y[2];

	// Middle Stencil
	double uM_y  = 0.5*U_y[3] -          0.5*U_y[1];
	double uM_yy = 0.5*U_y[1] - U_y[2] + 0.5*U_y[3];

    // Top stencil
	double uT_y  = -1.5*U_y[2] + 2.0*U_y[3] - 0.5*U_y[4];
	double uT_yy =  0.5*U_y[2] -     U_y[3] + 0.5*U_y[4];

    // Smoothness Indicators
	double IS_B = uB_y*uB_y + 4.333333333333333*(uB_yy*uB_yy); // Bottom stencil
	double IS_M = uM_y*uM_y + 4.333333333333333*(uM_yy*uM_yy); // Middle stencil
	double IS_T = uT_y*uT_y + 4.333333333333333*(uT_yy*uT_yy); // Top stencil

    double wB = 1.0/(std::pow((IS_B + epsilon), p));
	double wM = central_cell_wt/(std::pow((IS_M + epsilon), p));
	double wT = 1.0/(std::pow((IS_T + epsilon), p));

    sum = wB + wM + wT;

	// Normalize the weights

	wB = wB/sum; wM = wM/sum; wT = wT/sum;

    u_y = wB*uB_y + wM*uM_y + wT*uT_y;
	u_yy = wB*uB_yy + wM*uM_yy + wT*uT_yy;

    // Get u_xy

    double u_xy1 =  U_xy[1] - U_xy[0] - u_x - u_y - u_xx - u_yy;
    double u_xy2 = -U_xy[2] + U_xy[0] + u_x - u_y + u_xx + u_yy;
    double u_xy3 = -U_xy[3] + U_xy[0] - u_x + u_y + u_xx + u_yy;
    double u_xy4 =  U_xy[4] - U_xy[0] + u_x + u_y - u_xx - u_yy;

    double IS1 = 4.0*(u_xx*u_xx + u_yy*u_yy) + u_xy1*u_xy1;
    double IS2 = 4.0*(u_xx*u_xx + u_yy*u_yy) + u_xy2*u_xy2;
    double IS3 = 4.0*(u_xx*u_xx + u_yy*u_yy) + u_xy3*u_xy3;
    double IS4 = 4.0*(u_xx*u_xx + u_yy*u_yy) + u_xy4*u_xy4;

    double w1 = 1.0/(std::pow((IS1 + epsilon), p));
    double w2 = 1.0/(std::pow((IS2 + epsilon), p));
    double w3 = 1.0/(std::pow((IS3 + epsilon), p));
    double w4 = 1.0/(std::pow((IS4 + epsilon), p));

    sum = w1 + w2 + w3 + w4;

	// Normalize the weights

	w1 = w1/sum; w2 = w2/sum; w3 = w3/sum; w4 = w4/sum;

    u_xy = w1*u_xy1 + w2*u_xy2 + w3*u_xy3 + w4*u_xy4;

    coeffs[0] = U_xy[0]; coeffs[1] = u_x;
    coeffs[2] = u_y;     coeffs[3] = u_xx;
    coeffs[4] = u_yy;    coeffs[5] = u_xy;
}

// Fourth order WENO

void get_weno4_coefficients(const std::vector<double>& U_x, const std::vector<double>& U_y, const std::vector<double>& U_xy, std::vector<double>& coeffs) {

	double u_s1 = U_x[0];
	double u_p1 = U_x[1];
	double U_0 = U_x[2];
	double u_p3 = U_x[3];
	double u_s5 = U_x[4];
	double u_s7 = U_y[0];
	double u_p4 = U_y[1];
	double u_p2 = U_y[3];
	double u_s3 = U_y[4];
	double u_s4 = U_xy[1];
	double u_s6 = U_xy[2];
	double u_s2 = U_xy[3];
	double u_s8 = U_xy[4];

    static const double epsilon = 1.0e-12;
	double gammaHi = 0.85; double gammaLo = 0.85;

	double gammaR4 = gammaHi; double gammaR3[5];
	double wR4; double wR3[5];

	gammaR3[0] = (1.0 - gammaHi)*gammaLo;
	gammaR3[1] = 0.25*(1.0 - gammaHi)*(1.0 - gammaLo);
	gammaR3[2] = gammaR3[1]; gammaR3[3] = gammaR3[1]; gammaR3[4] = gammaR3[1];

	// ----------------------- fourth order Stencil -----------------------

	double u_xR4   = (-82.*u_p1 + 82.*u_p3 + 11.*u_s1 - 11.*u_s5)/(120.);
	double u_yR4   = (82.*u_p2 - 82.*u_p4 - 11.*u_s3 + 11.*u_s7)/(120.);
	double u_xxR4  = (-2.*U_0 + u_p1 + u_p3)/(2.);
	double u_yyR4  = (-2.*U_0 + u_p2 + u_p4)/(2.);
	double u_xyR4  = (-u_s2 + u_s4 - u_s6 + u_s8)/(4.);
	double u_xxxR4 = (2.*u_p1 - 2.*u_p3 - u_s1 + u_s5)/(12.);
	double u_yyyR4 = (-2.*u_p2 + 2.*u_p4 + u_s3 - u_s7)/(12.);
	double u_xxyR4 = (-2.*u_p2 + 2.*u_p4 + u_s2 + u_s4 - u_s6 - u_s8)/(4.);
	double u_xyyR4 = (2.*u_p1 - 2.*u_p3 - u_s2 + u_s4 + u_s6 - u_s8)/(4.);

	// Smoothness Indicator

    double ISR4 = (u_xR4 + 0.1*u_xxxR4)*(u_xR4 + 0.1*u_xxxR4) +
    		      (u_yR4 + 0.1*u_yyyR4)*(u_yR4 + 0.1*u_yyyR4)+
		          4.333333333333333*(u_yyR4*u_yyR4 + u_xxR4*u_xxR4) +
		          39.05*(u_yyyR4*u_yyyR4 + u_xxxR4*u_xxxR4) +
		          1.1666666666666667*u_xyR4*u_xyR4 +
		          0.47*(u_xxyR4*u_xxyR4 + u_xyyR4*u_xyyR4);

    // ----------------------- third order stencils -----------------------

	double u_xR3[5]; double u_yR3[5];  double u_xxR3[5]; double u_yyR3[5]; double u_xyR3[5]; double ISR3[5];

	// Centered

	u_xR3[0]  = (-u_p1 + u_p3)/(2.);
	u_yR3[0]  = (u_p2 - u_p4)/(2.);
	u_xxR3[0] = (-2.*U_0 + u_p1 + u_p3)/(2.);
	u_yyR3[0] = (-2.*U_0 + u_p2 + u_p4)/(2.);
	u_xyR3[0] = (-u_s2 + u_s4 - u_s6 + u_s8)/(4.);

	// Directional 1

	u_xR3[1]  = (-1.5*U_0 + 2.*u_p3 - 0.5*u_s5);
	u_yR3[1]  = (u_p2 - u_p4)/(2.);
	u_xxR3[1] = (U_0 - 2.*u_p3 + u_s5)/(2.);
	u_yyR3[1] = (-2.*U_0 + u_p2 + u_p4)/(2.);
	u_xyR3[1] = (-u_p2 + u_p4 + u_s4 - u_s6)/(2.);

	// Directional 2

	u_xR3[2]  = (3.*U_0 - 4.*u_p1 + u_s1)/(2.);
	u_yR3[2]  = (u_p2 - u_p4)/(2.);
	u_xxR3[2] = (U_0 - 2*u_p1 + u_s1)/(2.);
	u_yyR3[2] = (-2.*U_0 + u_p2 + u_p4)/(2.);
	u_xyR3[2] = (u_p2 - u_p4 - u_s2 + u_s8)/(2.);

	// Directional 3

	u_xR3[3]  = (-u_p1 + u_p3)/(2.);
	u_yR3[3]  = (-3.*U_0 + 4.*u_p2 - u_s3)/(2.);
	u_xxR3[3] = (-2.*U_0 + u_p1 + u_p3)/(2.);
	u_yyR3[3] = (U_0 - 2.*u_p2 + u_s3)/(2.);
	u_xyR3[3] = (u_p1 - u_p3 - u_s2 + u_s4)/(2.);

	// Directional 4

	u_xR3[4]  = (-u_p1 + u_p3)/(2.);
	u_yR3[4]  = (3.*U_0 - 4.*u_p4 + u_s7)/(2.);
	u_xxR3[4] = (-2.*U_0 + u_p1 + u_p3)/(2.);
	u_yyR3[4] = (U_0 - 2.*u_p4 + u_s7)/(2.);
	u_xyR3[4] = (-u_p1 + u_p3 - u_s6 + u_s8)/(2.);

	// Find smoothness indicators

	for (unsigned int i = 0; i < 5; ++i) {
		ISR3[i] = u_xR3[i]*u_xR3[i] + u_yR3[i]*u_yR3[i] +
				  4.333333333333333*(u_xxR3[i]*u_xxR3[i] + u_yyR3[i]*u_yyR3[i]) +
				  1.1666666666666667*u_xyR3[i]*u_xyR3[i];
	}
	/*
	double tau =
	0.2*(std::abs(ISR4 - ISR3[0]) + std::abs(ISR4 - ISR3[1]) + std::abs(ISR4 - ISR3[2]) + std::abs(ISR4 - ISR3[3]) + std::abs(ISR4 - ISR3[4]));
	*/
	
	//wR4 = gammaR4*(1.0 + (tau*tau)/((ISR4 + epsilon)*(ISR4 + epsilon)) );

	wR4 = gammaR4*(1.0/((ISR4 + epsilon)*(ISR4 + epsilon)));
	
	double sum = wR4;

	for (int i = 0; i < 5; ++i) {
		//wR3[i] = gammaR3[i]*(1.0 + (tau*tau)/((ISR3[i] + epsilon)*(ISR3[i] + epsilon)) );
		wR3[i] = gammaR3[i]*( 1.0/((ISR3[i] + epsilon)*(ISR3[i] + epsilon)) );
		sum += wR3[i];
	}


    wR4 = wR4/sum;

	for (int i = 0; i < 5; ++i)
		wR3[i] = wR3[i]/sum;


	double gamma_x, w_x, gamma_y, w_y, gamma_xx, w_xx, gamma_yy, w_yy, gamma_xy, w_xy;

    gamma_x = 0.0; w_x = 0.0; gamma_y = 0.0; w_y = 0.0;  gamma_xx = 0.0;
    w_xx = 0.0; gamma_yy = 0.0; w_yy = 0.0; gamma_xy = 0.0; w_xy = 0.0;

    for (int i = 0; i < 5; i++) {
        gamma_x += gammaR3[i]*u_xR3[i];
        w_x += wR3[i]*u_xR3[i];

		gamma_y += gammaR3[i]*u_yR3[i];
        w_y += wR3[i]*u_yR3[i];

		gamma_xx += gammaR3[i]*u_xxR3[i];
        w_xx += wR3[i]*u_xxR3[i];

		gamma_yy += gammaR3[i]*u_yyR3[i];
		w_yy += wR3[i]*u_yyR3[i];

		gamma_xy += gammaR3[i]*u_xyR3[i];
        w_xy += wR3[i]*u_xyR3[i];
    }

    const double wt_ratio = (wR4/gammaR4);

    coeffs[0] = U_0;
    coeffs[1] = wt_ratio*(u_xR4 - gamma_x) + w_x;
    coeffs[2] = wt_ratio*(u_yR4 - gamma_y) + w_y;
    coeffs[3] = wt_ratio*(u_xxR4 - gamma_xx) + w_xx;
    coeffs[4] = wt_ratio*(u_yyR4 - gamma_yy) + w_yy;
    coeffs[5] = wt_ratio*(u_xyR4 - gamma_xy) + w_xy;
    coeffs[6] = wt_ratio*u_xxxR4;
    coeffs[7] = wt_ratio*u_yyyR4;
    coeffs[8] = wt_ratio*u_xxyR4;
    coeffs[9] = wt_ratio*u_xyyR4;
}


void weno_2d(const std::vector<double>& U_x,
		   const std::vector<double>& U_y,
		   const std::vector<double>& U_xy,
		   std::vector<double>& coeffs, int order) {


	switch (order) {
	case 2:
		get_weno2_coefficients(U_x, U_y, U_xy, coeffs);
		break;
	case 3:
		get_weno3_coefficients(U_x, U_y, U_xy, coeffs);
		break;
	case 4:
		get_weno4_coefficients(U_x, U_y, U_xy, coeffs);
	    break;
	default:
		Assert(order <= 4, ErrNotImplemented());
		break;

	}
}


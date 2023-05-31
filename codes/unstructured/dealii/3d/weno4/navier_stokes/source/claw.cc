/*
 * claw.cc
 *      Author: sunder
 */

#include "../include/claw.h"
#include "../include/Headers.h"
#include "../include/Exceptions.h"

void ErrNegativePressureDensity(double d, double p, double x, double y) {

	if (d < 0.0 || p < 0.0) {
	
		std::cerr << "Negative Pressure/Density" << std::endl;
        std::cerr << "Density = " << d << ", Pressure = " << p << std::endl;
        std::cerr << "at (" << x << ", " << y << ")" << std::endl;
		
		std::exit(1); 
	}
}
// Conservation law definitions

//----------------------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------------------


Conservation_Law::Conservation_Law() :
	specific_heat_cp(0.0),
	specific_heat_cv(0.0),
	prandtl_no(0.0),
	specific_heat_ratio(0.0)

{}

//----------------------------------------------------------------------------
// Constructor taking specific heat at constant pressure, specific heat at constant
// volume and Prandtl number as inputs
//----------------------------------------------------------------------------

Conservation_Law::Conservation_Law(double cp, double cv, double pr) :
	specific_heat_cp(cp),
	specific_heat_cv(cv),
	prandtl_no(pr)
{
	specific_heat_ratio = cp/cv;
}

//----------------------------------------------------------------------------
// Get specific heat at constant pressure
//----------------------------------------------------------------------------

double Conservation_Law::c_p() const {
	return specific_heat_cp;
}

//----------------------------------------------------------------------------
// Get specific heat at constant volume
//----------------------------------------------------------------------------

double Conservation_Law::c_v() const {
	return specific_heat_cv;
}

//----------------------------------------------------------------------------
// Get ratio of specific heats
//----------------------------------------------------------------------------

double Conservation_Law::gamma() const {
	return specific_heat_ratio;
}

//----------------------------------------------------------------------------
// Get gas constant
//----------------------------------------------------------------------------

double Conservation_Law::R() const {
	return specific_heat_cp - specific_heat_cv;
}

//----------------------------------------------------------------------------
// Get Prandtl number
//----------------------------------------------------------------------------

double Conservation_Law::Pr() const {
	return prandtl_no;
}

//----------------------------------------------------------------------------
// Find heat conductivity of the gas at the given state U
//----------------------------------------------------------------------------

double Conservation_Law::heat_conductivity(const Vector<double>& U) {

	return (specific_heat_cp*viscosity(U)/prandtl_no);
}

double Conservation_Law::heat_conductivity(const double mu) {

	return (specific_heat_cp*mu/prandtl_no);
}

//----------------------------------------------------------------------------
// Find the viscosity of the gas at the given state U
//----------------------------------------------------------------------------

double Conservation_Law::viscosity(const Vector<double>& U)  {

	// Use Sutherland's law
/*	
	double T = Temperature(U);
	static const double mu_0 = 0.1;
	static const double T_0 = 1.0;
	static const double beta = 1.5;
	static const double s = 1.0;

	double mu = mu_0*std::pow(T/T_0, beta)*((T_0 + s)/(T + s));
*/	
	double mu = 14.0;
	return mu;
}

double Conservation_Law::viscosity(const double T)  {

	// Use Sutherland's law
/*	
	static const double mu_0 = 0.1;
	static const double T_0 = 1.0;
	static const double beta = 1.5;
	static const double s = 1.0;

	double mu = mu_0*std::pow(T/T_0, beta)*((T_0 + s)/(T + s));
*/	
	double mu = 14.0;
	return mu;
}

//----------------------------------------------------------------------------
// Given a conserved variable, find the pressure of the gas
//----------------------------------------------------------------------------

double Conservation_Law::Pressure(const Vector<double>& U) {

	return ( (gamma()-1.0)*(U[4] - 0.5*((U[1]*U[1] + U[2]*U[2] + U[3]*U[3])/U[0])) );
}

//----------------------------------------------------------------------------
// Given a conserved variable, find temperature of the gas
//----------------------------------------------------------------------------

double Conservation_Law::Temperature(const Vector<double>& U) {


	return ( Pressure(U)/(U[0]) );
}


//----------------------------------------------------------------------------
// Given a conserved variable, find the speed of sound
//----------------------------------------------------------------------------

double Conservation_Law::speed_of_sound(const Vector<double>& U) {

	return std::sqrt(specific_heat_ratio*Pressure(U)/U[0]);
}

//----------------------------------------------------------------------------
// Get conserved variable from a primitive variable
//----------------------------------------------------------------------------

void Conservation_Law::primitive_to_conserved(const Vector<double>& W, Vector<double>& U) const {
/*
  Assert(W.size() == 5,
         ExcDimensionMismatch(W.size(), 5) );
  Assert(U.size() == 5,
         ExcDimensionMismatch(U.size(), 5) );
*/
	//double p = W[3]*W[0]*R(); // p = rho*R*T

    double e = W[4]/((gamma() - 1.0)*W[0]);
    double k = 0.5*(W[1]*W[1] + W[2]*W[2] + W[3]*W[3]);

    U[0] = W[0];
    U[1] = W[0]*W[1];
    U[2] = W[0]*W[2];
    U[3] = W[0]*W[3];
    U[4] = W[0]*(k + e);
}

//----------------------------------------------------------------------------
// Get primitive variable from a conserved variable
//----------------------------------------------------------------------------

void Conservation_Law::conserved_to_primitive(const Vector<double>& U, Vector<double>& W) const {
/*
  Assert(W.size() == 5,
         ExcDimensionMismatch(W.size(), 5) );
  Assert(U.size() == 5,
         ExcDimensionMismatch(U.size(), 5) );
*/
	W[0] = U[0];
	W[1] = U[1]/U[0];
	W[2] = U[2]/U[0];
	W[3] = U[3]/U[0];
	W[4] = (gamma()-1.0)*(U[4] - 0.5*U[0]*(W[1]*W[1] + W[2]*W[2] + W[3]*W[3]));
}


void Conservation_Law::Usual_To_Conservative(double*& Var_Q, double* Var_U, double gamma_l, double P_Infty_l) {
        Var_Q[0] = Var_U[0] ; // rho
        Var_Q[1] = Var_U[0]*Var_U[1] ; // rho*u
        Var_Q[2] = Var_U[0]*Var_U[2] ; // rho*v
        Var_Q[3] = Var_U[0]*Var_U[3] ; // rho*w
        Var_Q[4] = Var_U[4]/(gamma_l - 1.0) + gamma_l*P_Infty_l/(gamma_l - 1.0) + 0.5*Var_U[0]*(Var_U[1]*Var_U[1] + Var_U[2]*Var_U[2] + Var_U[3]*Var_U[3]) ; // E
}

void Conservation_Law::FluxFunction(Vector<double>& FV, double* Q, double nx, double ny, double nz) {
	double un; // normal velocity
		un = (Q[1]*nx + Q[2]*ny + Q[3]*nz) ;
      FV[0] = Q[0]*un ;
      FV[1] = Q[0]*Q[1]*un + Q[4]*nx ;
      FV[2] = Q[0]*Q[2]*un + Q[4]*ny ;
      FV[3] = Q[0]*Q[3]*un + Q[4]*nz ;
      FV[4] = ( gamma()*(Q[4])/(gamma() - 1.0) + 0.5*Q[0]*(Q[1]*Q[1] + Q[2]*Q[2] + Q[3]*Q[3]) )*un ;
}

/*
Vector<double> Conservation_Law::local_Lax_Friedrichs_riemann_solver(const Vector<double>& UL, const Vector<double>& UR, double nx, double ny, double nz, const Point<3>& P, bool boundary) {    
    if (UL[0] < 0.0 || UR[0] < 0.0 || UL[4] < 0.0 || UR[4] < 0.0) {
        ThrowNegativePressureDensityError(UL, UR, P, boundary); 
    }

	double U_L[] = {UL[0], UL[1], UL[2], UL[3], UL[4]};
   double U_R[] = {UR[0], UR[1], UR[2], UR[3], UR[4]};
   Vector<double> F(5);
	double flux[5];
 	double S_star;
 	double P_Infty_L = 0.0;
 	double P_Infty_R = 0.0;
 	int i ;
	double rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, P_L, P_R, c_L, c_R ;
	double un_L, un_R, S_L, S_R ;
	//double ut_L, ut_R;
	double *FL, *FR, *U ;

	U = new double[5] ; FL = new double[5] ; FR = new double[5] ;

	rho_L = U_L[0] ; u_L = U_L[1] ; v_L = U_L[2] ; w_L = U_L[3] ; P_L = U_L[4] ;
  rho_R = U_R[0] ; u_R = U_R[1] ; v_R = U_R[2] ; w_R = U_R[3] ; P_R = U_R[4] ;
	un_L = u_L*nx + v_L*ny + w_L*nz; //ut_L = -u_L*ny + v_L*nx ;
	un_R = u_R*nx + v_R*ny + w_R*nz; //ut_R = -u_R*ny + v_R*nx ;

  c_L = sqrt(gamma()*(P_L + P_Infty_L)/rho_L) ; c_R = sqrt(gamma()*(P_R + P_Infty_R)/rho_R) ;

	
	S_L = std::max(fabs(un_R - c_R), fabs(un_L - c_L)) ;
	S_R = std::max(fabs(un_L + c_L), fabs(un_R + c_R)) ;
	S_star = std::max(fabs(S_L),fabs(S_R)) ;

	// Now compute the left right and starred fluxes.
	Usual_To_Conservative(U,U_L,gamma(),P_Infty_L) ;
	FluxFunction(FL,U_L,gamma(),P_Infty_L,nx,ny,nz) ;
        for(i = 0 ; i < 5 ; i++) flux[i] = 0.5*FL[i] + 0.5*S_star*U[i] ;

	Usual_To_Conservative(U,U_R,gamma(),P_Infty_R) ;
	FluxFunction(FR,U_R,gamma(),P_Infty_R,nx,ny,nz);
        for(i = 0 ; i < 5 ; i++) flux[i] += 0.5*FR[i] - 0.5*S_star*U[i] ;

    F[0] = flux[0]; F[1] = flux[1]; F[2] = flux[2]; F[3] = flux[3]; F[4] = flux[4];

	delete [] U ; delete [] FL ; delete [] FR ;
	return F;
}
*/

Vector<double> Conservation_Law::HLLC_riemann_solver(const Vector<double>& WL, const Vector<double>& WR, double nx, double ny, double nz, const Point<3>& P, bool boundary) {

    if (WL[0] < 0.0 || WR[0] < 0.0 || WL[4] < 0.0 || WR[4] < 0.0) {
//        ThrowNegativePressureDensityError(WL, WR, P, boundary); 
		std::cout<<"WL: "<<WL<<"WR: "<<WR<<"Point: "<<P<<std::endl;
    } 

     	int i ;
		Vector<double> F(5);
	    double flux[5];
		double U_L[] = {WL[0], WL[1], WL[2], WL[3], WL[4]};
	    double U_R[] = {WR[0], WR[1], WR[2], WR[3], WR[4]};
        double rho_L, rho_R, u_L, u_R, v_L, v_R, w_L, w_R, P_L, P_R, c_L, c_R, P_Star, E_L, E_R;
        double un_L, un_R, ut1_L, ut1_R, ut2_L, ut2_R ;
        double t1x, t1y, t1z, t2x, t2y, t2z ;
        double un, ut1, ut2 ;
        double S_L, S_R, S_star ;
        double *UL_star, *UR_star, *FL_star, *FR_star, *U, Local, **Rot_Mat, **InvRot_Mat ;
	    double P_Infty_L = 0.0, P_Infty_R = 0.0;
		double gamma_L = gamma(), gamma_R = gamma();
		Vector<double> FL(5), FR(5);
        UL_star = new double[5] ; UR_star = new double[5] ; U = new double[5] ;
       	FL_star = new double[5] ; FR_star = new double[5] ;
        Allocate_2D_R(Rot_Mat,3,3) ; Allocate_2D_R(InvRot_Mat,3,3) ;
               
        rho_L = U_L[0] ; u_L = U_L[1] ; v_L = U_L[2] ; w_L = U_L[3] ; P_L = U_L[4] ; 
        rho_R = U_R[0] ; u_R = U_R[1] ; v_R = U_R[2] ; w_R = U_R[3] ; P_R = U_R[4] ; 

        // Define the two tangent directions following Miller and Colella's JCP paper
        if( fabs(ny + nz) >  fabs(ny - nz) ) {

                Local = sqrt( 2.0*(1.0 + nz*(ny - nx) + nx*ny) ) ;
                t1x = (ny + nz)/Local ; t1y = (nz - nx)/Local ; t1z = -(nx + ny)/Local ;
                t2x = ( nx*(nz - ny) - ny*ny - nz*nz )/Local ; 
                t2y = ( ny*(nx + nz) + nx*nx + nz*nz )/Local ; 
                t2z = ( nz*(nx - ny) - nx*nx - ny*ny )/Local ; 

        } else {

                Local = sqrt( 2.0*(1.0 - nx*ny - nx*nz - ny*nz) ) ;
                t1x = (ny - nz)/Local ; t1y = (nz - nx)/Local ; t1z = (nx - ny)/Local ;
                t2x = ( nx*(ny + nz) - ny*ny - nz*nz )/Local ; 
                t2y = ( ny*(nx + nz) - nx*nx - nz*nz )/Local ;
                t2z = ( nz*(nx + ny) - nx*nx - ny*ny )/Local ;

        }
/*      cout << nx << "\t" << ny << "\t" << nz << endl ;
        cout << t1x << "\t" << t1y << "\t" << t1z << endl ;
        cout << t2x << "\t" << t2y << "\t" << t2z << endl ;*/
        Rot_Mat[0][0] = nx ; Rot_Mat[0][1] = ny ; Rot_Mat[0][2] = nz ;
        Rot_Mat[1][0] = t1x ; Rot_Mat[1][1] = t1y ; Rot_Mat[1][2] = t1z ;
        Rot_Mat[2][0] = t2x ; Rot_Mat[2][1] = t2y ; Rot_Mat[2][2] = t2z ;

        LUInverse(InvRot_Mat, Rot_Mat, 3) ;

        un_L = u_L*nx + v_L*ny + w_L*nz ; ut1_L = u_L*t1x + v_L*t1y + w_L*t1z ; ut2_L = u_L*t2x + v_L*t2y + w_L*t2z ;
        un_R = u_R*nx + v_R*ny + w_R*nz ; ut1_R = u_R*t1x + v_R*t1y + w_R*t1z ; ut2_R = u_R*t2x + v_R*t2y + w_R*t2z ;


        c_L = sqrt(gamma_L*(P_L + P_Infty_L)/rho_L) ; c_R = sqrt(gamma_R*(P_R + P_Infty_R)/rho_R) ;

        E_L = (P_L + gamma_L*P_Infty_L)/(gamma_L - 1.0) + 0.5*rho_L*u_L*u_L + 0.5*rho_L*v_L*v_L + 0.5*rho_L*w_L*w_L ;
        E_R = (P_R + gamma_R*P_Infty_R)/(gamma_R - 1.0) + 0.5*rho_R*u_R*u_R + 0.5*rho_R*v_R*v_R + 0.5*rho_R*w_R*w_R ;

        S_L = std::min((un_R - c_R), (un_L - c_L)) ;
        S_R = std::max((un_L + c_L), (un_R + c_R)) ;
        S_star = ( P_R - P_L + rho_L*un_L*(S_L-un_L) - rho_R*un_R*(S_R - un_R) )/(rho_L*(S_L - un_L) - rho_R*(S_R - un_R)) ;
        P_Star = P_L + rho_L*(S_star - un_L)*(S_L - un_L) ;

        // Now compute the left right and starred fluxes.
        FluxFunction(FL,U_L,nx,ny,nz) ; FluxFunction(FR,U_R,nx,ny,nz) ;
//		FluxFunction(FL,U_L,gamma(),P_Infty_,nx,ny,nz) ; FluxFunction(FR,U_R,gamma(),P_Infty_,nx,ny,nz);

        UL_star[0] = rho_L*(S_L - un_L)/(S_L - S_star)  ;
        un = S_star ; ut1 = ut1_L ; ut2 = ut2_L ;
        UL_star[1] = UL_star[0]*(un*InvRot_Mat[0][0] + ut1*InvRot_Mat[0][1] + ut2*InvRot_Mat[0][2]) ;
        UL_star[2] = UL_star[0]*(un*InvRot_Mat[1][0] + ut1*InvRot_Mat[1][1] + ut2*InvRot_Mat[1][2]) ;
        UL_star[3] = UL_star[0]*(un*InvRot_Mat[2][0] + ut1*InvRot_Mat[2][1] + ut2*InvRot_Mat[2][2]) ;
        UL_star[4] = UL_star[0]*( (E_L/rho_L) + (S_star - un_L)*(S_star + P_L/(rho_L*(S_L - un_L)) ) )  ;

        UR_star[0] = rho_R*(S_R - un_R)/(S_R - S_star) ;
        un = S_star ; ut1 = ut1_R ; ut2 = ut2_R ;
        UR_star[1] = UR_star[0]*(un*InvRot_Mat[0][0] + ut1*InvRot_Mat[0][1] + ut2*InvRot_Mat[0][2]) ;
        UR_star[2] = UR_star[0]*(un*InvRot_Mat[1][0] + ut1*InvRot_Mat[1][1] + ut2*InvRot_Mat[1][2]) ;
        UR_star[3] = UR_star[0]*(un*InvRot_Mat[2][0] + ut1*InvRot_Mat[2][1] + ut2*InvRot_Mat[2][2]) ;
        UR_star[4] = UR_star[0]*( (E_R/rho_R) + (S_star - un_R)*(S_star + P_R/(rho_R*(S_R - un_R)) ) ) ;

        Usual_To_Conservative(U,U_L,gamma_L,P_Infty_L) ;
        for(i = 0 ; i < 5 ; i++) FL_star[i] = FL[i] + S_L*(UL_star[i] - U[i]) ;

        Usual_To_Conservative(U,U_R,gamma_R,P_Infty_R) ;
        for(i = 0 ; i < 5 ; i++) FR_star[i] = FR[i] + S_R*(UR_star[i] - U[i]) ;

        if( S_L > 0.0 ) {
                for(i = 0 ; i < 5 ; i++) flux[i] = FL[i] ;
        } else if((S_star >= 0.0) && (S_L < 0.0)) {
                for(i = 0 ; i < 5 ; i++) flux[i] = FL_star[i] ;
        } else if((S_star < 0.0) && (S_R >= 0.0)) {
                for(i = 0 ; i < 5 ; i++) flux[i] = FR_star[i] ;
        } else if(S_R < 0.0) {
                for(i = 0 ; i < 5 ; i++) flux[i] = FR[i] ;
        }

	    F[0] = flux[0]; F[1] = flux[1]; F[2] = flux[2]; F[3] = flux[3]; F[4] = flux[4];

        for(i = 0 ; i < 3 ; i++) { delete [] Rot_Mat[i] ; delete [] InvRot_Mat[i] ;}
        delete [] UL_star ; delete [] UR_star ;  delete [] U ;
        delete [] FL_star ; delete [] FR_star ;
        delete [] Rot_Mat ; delete [] InvRot_Mat ;

    return F;

}

//Vector<double> rotated_HLLC_riemann_solver(Vector<double> UL, Vector<double> UR, double nx, double ny, Point<2> P, bool boundary) {
Vector<double> Conservation_Law::rotated_HLLC_riemann_solver(const Vector<double>& UL, const Vector<double>& UR, double nx, double ny, double nz, const Point<3>& P, bool boundary) {

    Vector<double> flux1(5);
    Vector<double> flux2(5);
    Vector<double> flux3(5);
    Vector<double> Flux(5), WL(5), WR(5);

    int i ;
	double alpha1, alpha2, alpha3, n1x, n1y, n1z, t1x, t1y, t1z, t2x, t2y, t2z, u_L, u_R, v_L, v_R, w_L, w_R, du, dv, dw, dq, Local ;

	conserved_to_primitive(UL,WL);
	conserved_to_primitive(UR,WR);

	u_L = WL[1] ; v_L = WL[2] ; w_L = WL[3] ;
    u_R = WR[1] ; v_R = WR[2] ; w_R = WR[3] ;
	du = u_R - u_L ; dv = v_R - v_L ; dw = w_R - w_L ;
	dq = std::sqrt(du*du + dv*dv + dw*dw) ;

	if(dq < 1.0E-10) { n1x = nx ; n1y = ny ; n1z = nz;}
	else { n1x = du/dq ; n1y = dv/dq ; n1z = dw/dq ;}

	alpha1 = (n1x*nx + n1y*ny + n1z*nz) ;
	if(alpha1 < 0) { n1x = -n1x ; n1y = -n1y ; n1z = -n1z ; alpha1 = -alpha1 ; }

    if( std::fabs(n1y + n1z) >  std::fabs(n1y - n1z) ) {
	    Local = std::sqrt( 2.0*(1.0 + n1z*(n1y - n1x) + n1x*n1y) ) ;
        t1x = (n1y + n1z)/Local ; t1y = (n1z - n1x)/Local ; t1z = -(n1x + n1y)/Local ;
        t2x = ( n1x*(n1z - n1y) - n1y*n1y - n1z*n1z )/Local ; 
        t2y = ( n1y*(n1x + n1z) + n1x*n1x + n1z*n1z )/Local ; 
        t2z = ( n1z*(n1x - n1y) - n1x*n1x - n1y*n1y )/Local ; 

    } else {

		Local = std::sqrt( 2.0*(1.0 - n1x*n1y - n1x*n1z - n1y*n1z) ) ;
        t1x = (n1y - n1z)/Local ; t1y = (n1z - n1x)/Local ; t1z = (n1x - n1y)/Local ;
        t2x = ( n1x*(n1y + n1z) - n1y*n1y - n1z*n1z )/Local ; 
        t2y = ( n1y*(n1x + n1z) - n1x*n1x - n1z*n1z )/Local ;
        t2z = ( n1z*(n1x + n1y) - n1x*n1x - n1y*n1y )/Local ;

    }

	alpha2 = (t1x*nx + t1y*ny + t1z*nz) ; 
	if(alpha2 < 0) { t1x = -t1x ; t1y = -t1y ; t1z = -t1z ; alpha2 = -alpha2 ; }

	alpha3 = (t2x*nx + t2y*ny + t2z*nz) ; 
	if(alpha3 < 0) { t2x = -t2x ; t2y = -t2y ; t2z = -t2z ; alpha3 = -alpha3 ; }

	flux1 = HLLC_riemann_solver(WL, WR, n1x, n1y, n1z, P, boundary);
	flux2 = HLLC_riemann_solver(WL, WR, t1x, t1y, t1z, P, boundary);
	flux3 = HLLC_riemann_solver(WL, WR, t2x, t2y, t2z, P, boundary);

	for(i = 0 ; i < 5 ; i++) Flux[i] = alpha1*flux1[i] + alpha2*flux2[i] + alpha3*flux3[i] ;

	return Flux;

}


void Conservation_Law::compute_normal_viscous_flux(const Vector<double>& U,
	const std::vector<Vector<double> >& gradU,
	double nx, double ny, double nz, 
	Vector<double>& Flux) {

	Vector<double> W(5);
	conserved_to_primitive(U,W);
	double rho = W(0), u = W(1), v = W(2), w = W(3), p = W(4); 
	double T = p/(U[0]);

//	double mu = viscosity_from_primitive_variable(W);
//	double k = heat_conductivity_from_primitive_variable(W);
	double mu = viscosity(T);
//	double k = get_heat_conductivity(T);
	double k = heat_conductivity(mu);
	double rhou = U(1), rhov = U(2), rhow = U(3), e = U(4);
    double irho = 1.0/rho; 
	double rho2 = rho*rho;

	double rho_x  = gradU[0](0), rho_y  = gradU[0](1), rho_z  = gradU[0](2);
	double rhou_x = gradU[1](0), rhou_y = gradU[1](1), rhou_z = gradU[1](2);
	double rhov_x = gradU[2](0), rhov_y = gradU[2](1), rhov_z = gradU[2](2);
	double rhow_x = gradU[3](0), rhow_y = gradU[3](1), rhow_z = gradU[3](2);
	double e_x    = gradU[4](0),    e_y = gradU[4](1),    e_z = gradU[4](2);

	double u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z, T_x, T_y, T_z;
/*
	u = M_x/rho  
	v = M_y/rho
	w = M_z/rho
	p = (gamma-1)*(E - Rational(1,2)*(M_x**2 + M_y**2 + M_z**2 )/rho)         
	T = p/(rho) 

	div_v = u_x + v_y + w_z
	E = p/(GAMMA - 1) + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho
	p = (GAMMA - 1)*(E - 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho)
	p = rho*T
	T = ((GAMMA - 1)/rho)*(E - 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho)
*/
	u_x = irho*(rhou_x - u*rho_x);
	v_x = irho*(rhov_x - v*rho_x);
	w_x = irho*(rhow_x - w*rho_x);
	T_x = (-(gamma() - 1.0)*rho_x/rho2) * (e -0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho)
		  + ((gamma() - 1.0)/rho)*(e_x -(rhou*rhou_x + rhov*rhov_x + rhow*rhow_x)/rho + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)*rho_x/rho2);

	u_y = irho*(rhou_y - u*rho_y);
	v_y = irho*(rhov_y - v*rho_y);
	w_y = irho*(rhow_y - w*rho_y);
	T_y = (-(gamma() - 1.0)*rho_y/rho2) * (e -0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho)
		  + ((gamma() - 1.0)/rho)*(e_y -(rhou*rhou_y + rhov*rhov_y + rhow*rhow_y)/rho + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)*rho_y/rho2);

	u_z = irho*(rhou_z - u*rho_z);
	v_z = irho*(rhov_z - v*rho_z);
	w_z = irho*(rhow_z - w*rho_z);
	T_z = (-(gamma() - 1.0)*rho_z/rho2) * (e -0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho)
		  + ((gamma() - 1.0)/rho)*(e_z -(rhou*rhou_z + rhov*rhov_z + rhow*rhow_z)/rho + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)*rho_z/rho2);

	double tau_xx, tau_xy, tau_xz, tau_yz, tau_yy, tau_zz;

	tau_xx = (2.0*mu/3.0)*( 2.0*u_x - v_y - w_z) ;
	tau_xy = mu*(u_y + v_x) ;
	tau_xz = mu*(u_z + w_x) ;
	tau_yy = (2.0*mu/3.0)*( 2.0*v_y - u_x - w_z) ;
	tau_yz = mu*(v_z + w_y) ;
	tau_zz = (2.0*mu/3.0)*( 2.0*w_z - u_x - v_y) ;

	double q_x = -k*T_x;
	double q_y = -k*T_y;
	double q_z = -k*T_z;

	Vector<double> F(5), G(5), H(5);  // x, y and z viscous flux
//	double mu_2_3 = mu*2.0/3.0;

	// x viscous flux
	F(0) = 0.0;
	F(1) = tau_xx;
	F(2) = tau_xy;
	F(3) = tau_xz;
	F(4) = u*tau_xx + v*tau_xy + w*tau_xz - q_x;

	// y viscous flux
	G(0) = 0.0;
	G(1) = tau_xy;
	G(2) = tau_yy;
	G(3) = tau_yz;
	G(4) = u*tau_xy + v*tau_yy + w*tau_yz - q_y;

	// z viscous flux
	H(0) = 0.0;
	H(1) = tau_xz;
	H(2) = tau_yz;
	H(3) = tau_zz;
	H(4) = u*tau_xz + v*tau_yz + w*tau_zz - q_z;

	for(unsigned int i = 0; i < 5; ++i)
		Flux(i) = nx*F(i) + ny*G(i) + nz*H(i);
}

Vector<double> Conservation_Law::viscous_reiman_solver(const Vector<double>& UL, const Vector<double>& UR, const std::vector<Vector<double> >& gradUL, const std::vector<Vector<double> >& gradUR, double nx, double ny, double nz) {

	Vector<double> Flux(5), FL(5), FR(5);

	compute_normal_viscous_flux(UL, gradUL, nx, ny, nz, FL);
	compute_normal_viscous_flux(UR, gradUR, nx, ny, nz, FR);

	for(unsigned int i = 0; i < 5; ++i)
		Flux(i) = 0.5*(FL(i) + FR(i));

	return Flux;

}

Vector<double> Conservation_Law::NS_riemann_solver(const Vector<double>& UL, const Vector<double>& UR, const std::vector<Vector<double> >& gradUL, const std::vector<Vector<double> >& gradUR, double nx, double ny, double nz, Point<3> P, bool boundary) {

	Vector<double> viscous_flux(5), convective_flux(5), Flux(5);
	convective_flux = rotated_HLLC_riemann_solver(UL, UR, nx, ny, nz, P, boundary);
	viscous_flux = viscous_reiman_solver(UL, UR, gradUL, gradUR, nx, ny, nz);

	for(unsigned int i = 0; i < 5; ++i)
		Flux(i) = convective_flux(i) - viscous_flux(i);

	return Flux;
}

void Conservation_Law::additional_data(
		const Vector<double>& U,
		const std::vector< Vector<double> >& gradU,
		double& p, std::vector< Vector<double> >& stress) {

	// Convective flux

	Vector<double> W(5);
	conserved_to_primitive(U,W);
	double rho = W(0), u = W(1), v = W(2), w = W(3), pres = W(4); 
	double T = p/(U[0]);
	p = pres;
	double mu = viscosity(T);
	double k = heat_conductivity(mu);

	double rhou = U(1), rhov = U(2), rhow = U(3), e = U(4);
	double rho2 = rho*rho;

	double rho_x  = gradU[0](0), rho_y  = gradU[0](1), rho_z  = gradU[0](2);
	double rhou_x = gradU[1](0), rhou_y = gradU[1](1), rhou_z = gradU[1](2);
	double rhov_x = gradU[2](0), rhov_y = gradU[2](1), rhov_z = gradU[2](2);
	double rhow_x = gradU[3](0), rhow_y = gradU[3](1), rhow_z = gradU[3](2);
    double irho = 1.0/rho; 
	double u_x, u_y, u_z, v_x, v_y, v_z, w_x, w_y, w_z;
/*
	u = M_x/rho  
	v = M_y/rho
	w = M_z/rho
	p = (gamma-1)*(E - Rational(1,2)*(M_x**2 + M_y**2 + M_z**2 )/rho)         
	T = p/(rho) 

	div_v = u_x + v_y + w_z
	E = p/(GAMMA - 1) + 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho
	p = (GAMMA - 1)*(E - 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho)
	p = rho*T
	T = ((GAMMA - 1)/rho)*(E - 0.5*(rhou*rhou + rhov*rhov + rhow*rhow)/rho)
*/
	u_x = irho*(rhou_x - u*rho_x);
	v_x = irho*(rhov_x - v*rho_x);
	w_x = irho*(rhow_x - w*rho_x);

	u_y = irho*(rhou_y - u*rho_y);
	v_y = irho*(rhov_y - v*rho_y);
	w_y = irho*(rhow_y - w*rho_y);

	u_z = irho*(rhou_z - u*rho_z);
	v_z = irho*(rhov_z - v*rho_z);
	w_z = irho*(rhow_z - w*rho_z);

	double tau_xx, tau_xy, tau_xz, tau_yz, tau_yy, tau_zz;

	tau_xx = (2.0*mu/3.0)*( 2.0*u_x - v_y - w_z) ;
	tau_xy = mu*(u_y + v_x) ;
	tau_xz = mu*(u_z + w_x) ;
	tau_yy = (2.0*mu/3.0)*( 2.0*v_y - u_x - w_z) ;
	tau_yz = mu*(v_z + w_y) ;
	tau_zz = (2.0*mu/3.0)*( 2.0*w_z - u_x - v_y) ;

	stress[0](0) = tau_xx;
	stress[0](1) = tau_xy;
	stress[0](2) = tau_xz;
	stress[1](0) = tau_xy;
	stress[1](1) = tau_yy;
	stress[1](2) = tau_yz;
	stress[2](0) = tau_xz;
	stress[2](1) = tau_yz;
	stress[2](2) = tau_zz;

	// Viscous flux (for derivations see the sympy notebook)

}

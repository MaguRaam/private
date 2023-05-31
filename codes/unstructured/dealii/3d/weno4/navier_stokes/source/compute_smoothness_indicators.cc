#include "../include/Weno432.h"

// Find smoothness indicators 

// Second order polynomial 

double compute_second_order_smoothness_indicator(Vector<double> coeffs, Vector<double> IS, double h) {

    double u_x = coeffs(0); 
	double u_y = coeffs(1);
	double u_z = coeffs(2); 
    
    double const_term = IS(0)*(1./h)*(u_x*u_x + u_y*u_y + u_z*u_z);
    
    return (const_term);
	
}

// Third Order polynomial 

double compute_third_order_smoothness_indicator(Vector<double> coeffs, Vector<double> IS, double h) {
    
    double u_x  = coeffs(0); 	double u_y  = coeffs(1); 	double u_z  = coeffs(2);
	double u_xx = coeffs(3); 	double u_yy = coeffs(4); 	double u_zz = coeffs(5);
    double u_xy = coeffs(6); 	double u_xz = coeffs(7); 	double u_yz = coeffs(8); 
    
	double h_inv = 1.0/h;
	
    double const_term = IS(0)*(h_inv)*(u_x*u_x + 4.0*u_xx*u_xx + u_xy*u_xy + u_xz*u_xz + u_y*u_y + 4.0*u_yy*u_yy + u_yz*u_yz + u_z*u_z + 4.0*u_zz*u_zz);
//    double const_term = IS(0)*(h_inv)*(4.0*u_xx*u_xx + u_xy*u_xy + u_xz*u_xz + 4.0*u_yy*u_yy + u_yz*u_yz + 4.0*u_zz*u_zz + u_x*u_x + u_y*u_y + u_z*u_z);


    double x_term = IS(1)*(h_inv)*(4.0*u_x*u_xx + 2.0*u_xy*u_y + 2.0*u_xz*u_z);
//    double x_term = IS(1)*(h_inv)*(4.0*u_x*u_xx + 2.0*u_y*u_xy + 2.0*u_z*u_xz);


    double y_term = IS(2)*(h_inv)*(2.0*u_x*u_xy + 4.0*u_y*u_yy + 2.0*u_yz*u_z);
//    double y_term = IS(2)*(h_inv)*(2.0*u_x*u_xy + 4.0*u_y*u_yy + 2.0*u_z*u_yz); 

    double z_term = IS(3)*(h_inv)*(2.0*u_x*u_xz + 2.0*u_y*u_yz + 4.0*u_z*u_zz);
//    double z_term = IS(3)*(h_inv)*(2.0*u_x*u_xy + 2.0*u_y*u_yz + 4.0*u_z*u_zz); 


    double x2_term = IS(4)*(h_inv)*(4.0*u_xx*u_xx + u_xy*u_xy + u_xz*u_xz);
//    double x2_term = IS(4)*(h_inv)*(4.0*u_xx*u_xx + u_xy*u_xy + u_xz*u_xz);

    double y2_term = IS(5)*(h_inv)*(u_xy*u_xy + 4.0*u_yy*u_yy + u_yz*u_yz);
//    double y2_term = IS(5)*(h_inv)*(u_xy*u_xy + 4.0*u_yy*u_yy + u_yz*u_yz); 

    double z2_term = IS(6)*(h_inv)*(u_xz*u_xz + u_yz*u_yz + 4.0*u_zz*u_zz);
//    double z2_term = IS(6)*(h_inv)*(u_xz*u_xz + u_yz*u_yz + 4.0*u_zz*u_zz); 

    double xy_term = IS(7)*(h_inv)*(4.0*u_xx*u_xy + 4.0*u_xy*u_yy + 2.0*u_xz*u_yz);
//    double xy_term = IS(7)*(h_inv)*(4.0*u_xx*u_xy + 4.0*u_xy*u_yy + 2.0*u_xz*u_yz);

    double xz_term = IS(8)*(h_inv)*(4.0*u_xx*u_xz + 2.0*u_xy*u_yz + 4.0*u_xz*u_zz);
//    double xz_term = IS(8)*(h_inv)*(4.0*u_xx*u_xz + 4.0*u_xz*u_zz + 2.0*u_xy*u_yz);

    double yz_term = IS(9)*(h_inv)*(2.0*u_xy*u_xz + 4.0*u_yy*u_yz + 4.0*u_yz*u_zz);
//    double yz_term = IS(9)*(h_inv)*(4.0*u_yy*u_yz + 4.0*u_yz*u_zz + 2.0*u_xy*u_xz);
    
    return (const_term + x_term + y_term + z_term + x2_term + y2_term + z2_term + xy_term + xz_term + yz_term); 
 
}
// Fourth order polynomial 

double compute_fourth_order_smoothness_indicator(Vector<double> coeffs, Vector<double> IS, double h) {

    double u_x   = coeffs(0); 	double u_y   = coeffs(1); 	double u_z   = coeffs(2);
	double u_xx  = coeffs(3); 	double u_yy  = coeffs(4); 	double u_zz  = coeffs(5);
    double u_xy  = coeffs(6); 	double u_xz  = coeffs(7); 	double u_yz  = coeffs(8);
	double u_xxx = coeffs(9); 	double u_yyy = coeffs(10); 	double u_zzz = coeffs(11);
    double u_xxy = coeffs(12);	double u_xxz = coeffs(13);	double u_yyx = coeffs(14);
	double u_yyz = coeffs(15);	double u_zzx = coeffs(16);	double u_zzy = coeffs(17);
	double u_xyz = coeffs(18);
	
	double h_inv = 1.0/h;
	double const_term = IS(0)*(h_inv)*(u_x*u_x + 4.0*u_xx*u_xx + 36.0*u_xxx*u_xxx + 4.0*u_xxy*u_xxy + 4.0*u_xxz*u_xxz + u_xy*u_xy + u_xyz*u_xyz + u_xz*u_xz + u_y*u_y + 4.0*u_yy*u_yy + 4.0*u_yyx*u_yyx + 36.0*u_yyy*u_yyy + 4.0*u_yyz*u_yyz + u_yz*u_yz + u_z*u_z + 4.0*u_zz*u_zz + 4.0*u_zzx*u_zzx + 4.0*u_zzy*u_zzy + 36.0*u_zzz*u_zzz);
//	double const_term = IS(0)*(h_inv)*(36.0*u_xxx*u_xxx + 4.0*u_xxy*u_xxy + 4.0*u_xxz*u_xxz + u_xyz*u_xyz + 4.0*u_yyx*u_yyx + 36.0*u_yyy*u_yyy + 4.0*u_yyz*u_yyz + 4.0*u_zzx*u_zzx + 4.0*u_zzy*u_zzy + 36.0*u_zzz*u_zzz + 4.0*u_xx*u_xx + u_xy*u_xy + u_xz*u_xz + 4.0*u_yy*u_yy + u_yz*u_yz + 4.0*u_zz*u_zz + u_x*u_x + u_y*u_y + u_z*u_z);


    double x_term = IS(1)*(h_inv)*(4.0*u_x*u_xx + 24.0*u_xx*u_xxx + 4.0*u_xxy*u_xy + 4.0*u_xxz*u_xz + 2.0*u_xy*u_y + 2.0*u_xyz*u_yz + 2.0*u_xz*u_z + 8.0*u_yy*u_yyx + 8.0*u_zz*u_zzx);
//    double x_term = IS(1)*(h_inv)*(24.0*u_xx*u_xxx + 4.0*u_xxy*u_xy + 4.0*u_xxz*u_xz + 2.0*u_xyz*u_yz + 8.0*u_yyx*u_yy + 8.0*u_zzx*u_zz + 4.0*u_x*u_xx + 2.0*u_xy*u_y + 2.0*u_xz*u_z);


    double y_term = IS(2)*(h_inv)*(2.0*u_x*u_xy + 8.0*u_xx*u_xxy + 4.0*u_xy*u_yyx + 2.0*u_xyz*u_xz + 4.0*u_y*u_yy + 24.0*u_yy*u_yyy + 4.0*u_yyz*u_yz + 2.0*u_yz*u_z + 8.0*u_zz*u_zzy);
//    double y_term = IS(2)*(h_inv)*(24.0*u_yy*u_yyy + 4.0*u_yyx*u_xy + 4.0*u_yyz*u_yz + 2.0*u_xyz*u_xz + 8.0*u_xxy*u_xx + 8.0*u_zzy*u_zz + 4.0*u_y*u_yy + 2.0*u_xy*u_x + 2.0*u_yz*u_z);


    double z_term = IS(3)*(h_inv)*(2.0*u_x*u_xz + 8.0*u_xx*u_xxz + 2.0*u_xy*u_xyz + 4.0*u_xz*u_zzx + 2.0*u_y*u_yz + 8.0*u_yy*u_yyz + 4.0*u_yz*u_zzy + 4.0*u_z*u_zz + 24.0*u_zz*u_zzz);
//    double z_term = IS(3)*(h_inv)*(24.0*u_zz*u_zzz + 4.0*u_zzy*u_yz + 4.0*u_zzx*u_xz + 2.0*u_xyz*u_xy + 8.0*u_yyz*u_yy + 8.0*u_xxz*u_xx + 4.0*u_z*u_zz + 2.0*u_yz*u_y + 2.0*u_xz*u_x);


    double x2_term = IS(4)*(h_inv)*(6.0*u_x*u_xxx + 4.0*u_xx*u_xx + 36.0*u_xxx*u_xxx + 4.0*u_xxy*u_xxy + 2.0*u_xxy*u_y + 4.0*u_xxz*u_xxz + 2.0*u_xxz*u_z + u_xy*u_xy + u_xyz*u_xyz + u_xz*u_xz + 4.0*u_yyx*u_yyx + 4.0*u_zzx*u_zzx);
//    double x2_term = IS(4)*(h_inv)*(36.0*u_xxx*u_xxx + 4.0*u_xxy*u_xxy + 4.0*u_xxz*u_xxz + u_xyz*u_xyz + 4.0*u_yyx*u_yyx + 4.0*u_zzx*u_zzx + 6.0*u_x*u_xxx + 4.0*u_xx*u_xx + 2.0*u_xxy*u_y + 2.0*u_xxz*u_z + u_xy*u_xy + u_xz*u_xz);


    double y2_term = IS(5)*(h_inv)*(2.0*u_x*u_yyx + 4.0*u_xxy*u_xxy + u_xy*u_xy + u_xyz*u_xyz + 6.0*u_y*u_yyy + 4.0*u_yy*u_yy + 4.0*u_yyx*u_yyx + 36.0*u_yyy*u_yyy + 4.0*u_yyz*u_yyz + 2.0*u_yyz*u_z + u_yz*u_yz + 4.0*u_zzy*u_zzy);
//    double y2_term = IS(5)*(h_inv)*(36.0*u_yyy*u_yyy + 4.0*u_yyx*u_yyx + 4.0*u_yyz*u_yyz + u_xyz*u_xyz + 4.0*u_xxy*u_xxy + 4.0*u_zzy*u_zzy + 6.0*u_y*u_yyy + 4.0*u_yy*u_yy + 2.0*u_yyx*u_x + 2.0*u_yyz*u_z + u_xy*u_xy + u_yz*u_yz);


    double z2_term = IS(6)*(h_inv)*(2.0*u_x*u_zzx + 4.0*u_xxz*u_xxz + u_xyz*u_xyz + u_xz*u_xz + 2.0*u_y*u_zzy + 4.0*u_yyz*u_yyz + u_yz*u_yz + 6.0*u_z*u_zzz + 4.0*u_zz*u_zz + 4.0*u_zzx*u_zzx + 4.0*u_zzy*u_zzy + 36.0*u_zzz*u_zzz);
//    double z2_term = IS(6)*(h_inv)*(36.0*u_zzz*u_zzz + 4.0*u_zzy*u_zzy + 4.0*u_zzx*u_zzx + u_xyz*u_xyz + 4.0*u_yyz*u_yyz + 4.0*u_xxz*u_xxz + 6.0*u_z*u_zzz + 4.0*u_zz*u_zz + 2.0*u_zzy*u_y + 2.0*u_zzx*u_x + u_yz*u_yz + u_xz*u_xz);


    double xy_term = IS(7)*(h_inv)*(4.0*u_x*u_xxy + 4.0*u_xx*u_xy + 24.0*u_xxx*u_xxy + 8.0*u_xxy*u_yyx + 4.0*u_xxz*u_xyz + 4.0*u_xy*u_yy + 4.0*u_xyz*u_yyz + 2.0*u_xyz*u_z + 2.0*u_xz*u_yz + 4.0*u_y*u_yyx + 24.0*u_yyx*u_yyy + 8.0*u_zzx*u_zzy);
//    double xy_term = IS(7)*(h_inv)*(24.0*u_xxx*u_xxy + 8.0*u_xxy*u_yyx + 4.0*u_xxz*u_xyz + 4.0*u_xyz*u_yyz + 24.0*u_yyx*u_yyy + 8.0*u_zzx*u_zzy + 4.0*u_x*u_xxy + 4.0*u_xx*u_xy + 4.0*u_xy*u_yy + 2.0*u_xyz*u_z + 2.0*u_xz*u_yz + 4.0*u_yyx*u_y);


    double xz_term = IS(8)*(h_inv)*(4.0*u_x*u_xxz + 4.0*u_xx*u_xz + 24.0*u_xxx*u_xxz + 4.0*u_xxy*u_xyz + 8.0*u_xxz*u_zzx + 2.0*u_xy*u_yz + 2.0*u_xyz*u_y + 4.0*u_xyz*u_zzy + 4.0*u_xz*u_zz + 8.0*u_yyx*u_yyz + 4.0*u_z*u_zzx + 24.0*u_zzx*u_zzz);
//    double xz_term = IS(8)*(h_inv)*(24.0*u_xxx*u_xxz + 8.0*u_xxz*u_zzx + 4.0*u_xxy*u_xyz + 4.0*u_xyz*u_zzy + 24.0*u_zzx*u_zzz + 8.0*u_yyx*u_yyz + 4.0*u_x*u_xxz + 4.0*u_xx*u_xz + 4.0*u_xz*u_zz + 2.0*u_xyz*u_y + 2.0*u_xy*u_yz + 4.0*u_zzx*u_z);


    double yz_term = IS(9)*(h_inv)*(2.0*u_x*u_xyz + 8.0*u_xxy*u_xxz + 2.0*u_xy*u_xz + 4.0*u_xyz*u_yyx + 4.0*u_xyz*u_zzx + 4.0*u_y*u_yyz + 4.0*u_yy*u_yz + 24.0*u_yyy*u_yyz + 8.0*u_yyz*u_zzy + 4.0*u_yz*u_zz + 4.0*u_z*u_zzy + 24.0*u_zzy*u_zzz);
//    double yz_term = IS(9)*(h_inv)*(24.0*u_yyy*u_yyz + 8.0*u_yyz*u_zzy + 4.0*u_yyx*u_xyz + 4.0*u_xyz*u_zzx + 24.0*u_zzy*u_zzz + 8.0*u_xxy*u_xxz + 4.0*u_y*u_yyz + 4.0*u_yy*u_yz + 4.0*u_yz*u_zz + 2.0*u_xyz*u_x + 2.0*u_xy*u_xz + 4.0*u_zzy*u_z);


    double x3_term = IS(10)*(h_inv)*(12.0*u_xx*u_xxx + 2.0*u_xxy*u_xy + 2.0*u_xxz*u_xz);
//    double x3_term = IS(10)*(h_inv)*(12.0*u_xxx*u_xx + 2.0*u_xy*u_xxy + 2.0*u_xxz*u_xz); 


    double y3_term = IS(11)*(h_inv)*(2.0*u_xy*u_yyx + 12.0*u_yy*u_yyy + 2.0*u_yyz*u_yz);
//    double y3_term = IS(11)*(h_inv)*(12.0*u_yyy*u_yy + 2.0*u_xy*u_yyx + 2.0*u_yyz*u_yz); 


    double z3_term = IS(12)*(h_inv)*(2.0*u_xz*u_zzx + 2.0*u_yz*u_zzy + 12.0*u_zz*u_zzz);
//    double z3_term = IS(12)*(h_inv)*(12.0*u_zzz*u_zz + 2.0*u_yz*u_zzy + 2.0*u_zzx*u_xz);


    double x2y_term = IS(13)*(h_inv)*(8.0*u_xx*u_xxy + 6.0*u_xxx*u_xy + 4.0*u_xxy*u_yy + 2.0*u_xxz*u_yz + 4.0*u_xy*u_yyx + 2.0*u_xyz*u_xz);
//    double x2y_term = IS(13)*(h_inv)*(8.0*u_xx*u_xxy + 6.0*u_xxx*u_xy + 4.0*u_xxy*u_yy + 2.0*u_xxz*u_yz + 4.0*u_xy*u_yyx + 2.0*u_xyz*u_xz);


    double x2z_term = IS(14)*(h_inv)*(8.0*u_xx*u_xxz + 6.0*u_xxx*u_xz + 2.0*u_xxy*u_yz + 4.0*u_xxz*u_zz + 2.0*u_xy*u_xyz + 4.0*u_xz*u_zzx);
//    double x2z_term = IS(14)*(h_inv)*(8.0*u_xx*u_xxz + 6.0*u_xxx*u_xz + 4.0*u_xxz*u_zz + 2.0*u_xxy*u_yz + 4.0*u_xz*u_zzx + 2.0*u_xyz*u_xy);


    double y2x_term = IS(15)*(h_inv)*(4.0*u_xx*u_yyx + 4.0*u_xxy*u_xy + 6.0*u_xy*u_yyy + 2.0*u_xyz*u_yz + 2.0*u_xz*u_yyz + 8.0*u_yy*u_yyx);
//    double y2x_term = IS(15)*(h_inv)*(8.0*u_yy*u_yyx + 6.0*u_yyy*u_xy + 4.0*u_yyx*u_xx + 2.0*u_yyz*u_xz + 4.0*u_xy*u_xxy + 2.0*u_xyz*u_yz);


    double y2z_term = IS(16)*(h_inv)*(2.0*u_xy*u_xyz + 2.0*u_xz*u_yyx + 8.0*u_yy*u_yyz + 6.0*u_yyy*u_yz + 4.0*u_yyz*u_zz + 4.0*u_yz*u_zzy);
//    double y2z_term = IS(16)*(h_inv)*(8.0*u_yy*u_yyz + 6.0*u_yyy*u_yz + 4.0*u_yyz*u_zz + 2.0*u_yyz*u_xz + 4.0*u_yz*u_zzy + 2.0*u_xyz*u_xy);


    double z2x_term = IS(17)*(h_inv)*(4.0*u_xx*u_zzx + 4.0*u_xxz*u_xz + 2.0*u_xy*u_zzy + 2.0*u_xyz*u_yz + 6.0*u_xz*u_zzz + 8.0*u_zz*u_zzx);
//    double z2x_term = IS(17)*(h_inv)*(8.0*u_zz*u_zzx + 6.0*u_zzz*u_xz + 4.0*u_zzx*u_xx + 2.0*u_zzy*u_xy + 4.0*u_xz*u_xxz + 2.0*u_xyz*u_yz);


    double z2y_term = IS(18)*(h_inv)*(2.0*u_xy*u_zzx + 2.0*u_xyz*u_xz + 4.0*u_yy*u_zzy + 4.0*u_yyz*u_yz + 6.0*u_yz*u_zzz + 8.0*u_zz*u_zzy);
//    double z2y_term = IS(18)*(h_inv)*(8.0*u_zz*u_zzy + 6.0*u_zzz*u_yz + 4.0*u_zzy*u_yy + 2.0*u_zzx*u_xy + 4.0*u_yz*u_yyz + 2.0*u_xyz*u_xz);


	double xyz_term = IS(19)*(h_inv)*(4.0*u_xx*u_xyz + 4.0*u_xxy*u_xz + 4.0*u_xxz*u_xy + 4.0*u_xy*u_yyz + 4.0*u_xyz*u_yy + 4.0*u_xyz*u_zz + 4.0*u_xz*u_zzy + 4.0*u_yyx*u_yz + 4.0*u_yz*u_zzx);
//	double xyz_term = IS(19)*(h_inv)*(4.0*u_xx*u_xyz + 4.0*u_xxy*u_xz + 4.0*u_xxz*u_xy + 4.0*u_xy*u_yyz + 4.0*u_xyz*u_yy + 4.0*u_xyz*u_zz + 4.0*u_xz*u_zzy + 4.0*u_yyx*u_yz + 4.0*u_yz*u_zzx);


    double x4_term = IS(20)*(h_inv)*(9.0*u_xxx*u_xxx + u_xxy*u_xxy + u_xxz*u_xxz);
//    double x4_term = IS(20)*(h_inv)*(9.0*u_xxx*u_xxx + u_xxy*u_xxy + u_xxz*u_xxz);


    double y4_term = IS(21)*(h_inv)*(u_yyx*u_yyx + 9.0*u_yyy*u_yyy + u_yyz*u_yyz);
//    double y4_term = IS(21)*(h_inv)*(9.0*u_yyy*u_yyy + u_yyx*u_yyx + u_yyz*u_yyz);


    double z4_term = IS(22)*(h_inv)*(u_zzx*u_zzx + u_zzy*u_zzy + 9.0*u_zzz*u_zzz);
//    double z4_term = IS(22)*(h_inv)*(9.0*u_zzz*u_zzz + u_zzx*u_zzx + u_zzy*u_zzy);


    double x3y_term = IS(23)*(h_inv)*(12.0*u_xxx*u_xxy + 4.0*u_xxy*u_yyx + 2.0*u_xxz*u_xyz);
//    double x3y_term = IS(23)*(h_inv)*(12.0*u_xxx*u_xxy + 4.0*u_xxy*u_yyx + 2.0*u_xxz*u_xyz);


    double x3z_term = IS(24)*(h_inv)*(12.0*u_xxx*u_xxz + 2.0*u_xxy*u_xyz + 4.0*u_xxz*u_zzx);
//    double x3z_term = IS(24)*(h_inv)*(12.0*u_xxx*u_xxz + 4.0*u_xxz*u_zzx + 2.0*u_xxy*u_xyz);


    double y3x_term = IS(25)*(h_inv)*(4.0*u_xxy*u_yyx + 2.0*u_xyz*u_yyz + 12.0*u_yyx*u_yyy);
//    double y3x_term = IS(25)*(h_inv)*(12.0*u_yyy*u_yyx + 4.0*u_yyx*u_xxy + 2.0*u_yyz*u_xyz);


    double y3z_term = IS(26)*(h_inv)*(2.0*u_xyz*u_yyx + 12.0*u_yyy*u_yyz + 4.0*u_yyz*u_zzy);
//    double y3z_term = IS(26)*(h_inv)*(12.0*u_yyy*u_yyz + 4.0*u_yyz*u_zzy + 2.0*u_yyx*u_xyz);


    double z3x_term = IS(27)*(h_inv)*(4.0*u_xxz*u_zzx + 2.0*u_xyz*u_zzy + 12.0*u_zzx*u_zzz);
//    double z3x_term = IS(27)*(h_inv)*(12.0*u_zzz*u_zzx + 4.0*u_zzx*u_xxz + 2.0*u_zzy*u_xyz);


    double z3y_term = IS(28)*(h_inv)*(2.0*u_xyz*u_zzx + 4.0*u_yyz*u_zzy + 12.0*u_zzy*u_zzz);
//    double z3y_term = IS(28)*(h_inv)*(12.0*u_zzz*u_zzy + 4.0*u_zzy*u_yyz + 2.0*u_zzx*u_xyz);


    double x2y2_term = IS(29)*(h_inv)*(6.0*u_xxx*u_yyx + 4.0*u_xxy*u_xxy + 6.0*u_xxy*u_yyy + 2.0*u_xxz*u_yyz + u_xyz*u_xyz + 4.0*u_yyx*u_yyx);
//    double x2y2_term = IS(29)*(h_inv)*(6.0*u_xxx*u_yyx + 4.0*u_xxy*u_xxy + 6.0*u_xxy*u_yyy + 2.0*u_xxz*u_yyz + u_xyz*u_xyz + 4.0*u_yyx*u_yyx); 


    double x2z2_term = IS(30)*(h_inv)*(6.0*u_xxx*u_zzx + 2.0*u_xxy*u_zzy + 4.0*u_xxz*u_xxz + 6.0*u_xxz*u_zzz + u_xyz*u_xyz + 4.0*u_zzx*u_zzx);
//    double x2z2_term = IS(30)*(h_inv)*(6.0*u_xxx*u_zzx + 4.0*u_xxz*u_xxz + 6.0*u_xxz*u_zzz + 2.0*u_xxy*u_zzy + u_xyz*u_xyz + 4.0*u_zzx*u_zzx); 


    double y2z2_term = IS(31)*(h_inv)*(u_xyz*u_xyz + 2.0*u_yyx*u_zzx + 6.0*u_yyy*u_zzy + 4.0*u_yyz*u_yyz + 6.0*u_yyz*u_zzz + 4.0*u_zzy*u_zzy);
//    double y2z2_term = IS(31)*(h_inv)*(6.0*u_yyy*u_zzy + 4.0*u_yyz*u_yyz + 6.0*u_yyz*u_zzz + 2.0*u_yyz*u_zzx + u_xyz*u_xyz + 4.0*u_zzy*u_zzy); 

    
    return (const_term + x_term + y_term + z_term + x2_term + y2_term + z2_term + xy_term + xz_term + yz_term 
			+ x3_term + y3_term + z3_term + x2y_term + x2z_term + y2x_term + y2z_term + z2x_term + z2y_term + xyz_term
			+ x4_term + y4_term + z4_term + x3y_term + x3z_term + y3x_term + y3z_term + z3x_term + z3y_term + x2y2_term + x2z2_term + y2z2_term); 
	
} 

double compute_smoothness_indicator(const Vector<double>& coeffs) {

	unsigned int n_terms = coeffs.size();
	double IS = 0.0;

	for(unsigned int i = 0; i < n_terms; ++i)
		IS += coeffs(i)*coeffs(i);

    return IS;

} 

/*
 * main.cc
 *      Author: sunder
 */ 

#include "path_conservative_fv.hh"

std::vector<double> air_helium(double x, double y) {

	std::vector<double> V(5); 
	
	// Air part 
	
	if (x > 2.8314606741573036) {
	
		V[0] = 1.92691;
		V[1] = -0.3336067590079454;
		V[2] = 0.0; 
		V[3] = 1.5698;
		V[4] = 1.0 - 1.0e-6; 
	}
	
	else {
		
		V[0] = 1.4;
		V[1] = 0.0;
		V[2] = 0.0; 
		V[3] = 1.0;
		V[4] = 1.0 - 1.0e-6; 
	}
	
	// Helium 
	
	double x_0 = 2.05247191011236; double y_0 = 0.5;
	double r = 0.2808988764044944; 
	
	if ((x-x_0)*(x-x_0) + (y-y_0)*(y-y_0) <= r*r ) {
	
		V[0] = 0.2546;
		V[1] = 0.0;
		V[2] = 0.0; 
		V[3] = 1.0;
		V[4] = 1.0e-6; 
		
	}
	
	return V;
	
}


int main() {
	
	
	AppCtx Params; 
	
	Params.N_x = 630;
	Params.N_y = 210;
	Params.order = 4;
	Params.CFL = 0.3;
	Params.x_min = 0.0;
	Params.x_max = 3.0;
	Params.y_min = 0.0;
	Params.y_max = 1.0;
	Params.final_time = 2.0;
	Params.write_interval = 100;
	Params.g1 = 1.4;
	Params.g2 = 1.648; 
	Params.p1 = 0.0;
	Params.p2 = 0.0;

	Params.left_boundary   = transmissive;
	Params.right_boundary  = transmissive;
	Params.top_boundary    = reflective;
	Params.bottom_boundary = reflective;
	
	Params.reconstruct_primitive_variables = true;
	
	Path_Conservative_FV Problem(air_helium, Params);
	Problem.run(); 
	

	return 0;
}

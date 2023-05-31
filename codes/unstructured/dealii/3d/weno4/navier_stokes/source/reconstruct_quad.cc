#include "../include/Weno432.h"


// Perform the actual reconstruction 

void Weno4_3D::reconstruct_quad() {
		
	DoFHandler<3>::active_cell_iterator cell;
	unsigned int g_i;
	double h;
	Point<3> face_quadrature_point;
	double u_slip, v_slip, w_slip, u_slip_max = 0, v_slip_max = 0, w_slip_max = 0, global_u_slip_max, global_v_slip_max, global_w_slip_max;
	double gradRho, gradRho_max = 0.0, gradE, gradE_max = 0.0;
	double h_max = 0, global_h_max;
	unsigned int first_order_cell = 0, max_g_i, max_c;
	bool no_slip_cell = false;
	Vector<double> grad_rho(3), grad_e(3);

	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];
		h = Cell[c].h();

        for (unsigned int f = 0; f < 6; ++f) {

			if(cell->face(f)->at_boundary()) {

                if (cell->face(f)->boundary_id() == 2) {
					no_slip_cell = true;
					if(h > h_max)
						h_max = h;

					if(first_order_cell == 0)
						if(is_1st_order[c]) {
							first_order_cell = 1;
							std::cout<<"g_i: "<<g_i<<"\t"<<cell->center()<<std::endl;
						}

    	            for (unsigned int q = 0; q < 4; q++) {
						face_quadrature_point = Cell[c].face_quadrature_point(f, q);
			            u_slip = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point, h);
			            v_slip = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point, h);
			            w_slip = evaluate_weno_polynomial(coeffs_RHO_W[c], WENO_poly_consts[c], face_quadrature_point, h);
//						std::cout<<"g_i: "<<g_i<<"\tu: "<<u_slip<<"\tcoeff: "<<coeffs_RHO_U[c]<<"upwind: "<<upwind_coeffs_RHO_U[c]<<std::endl;					
//			            u_slip = evaluate_weno_polynomial(upwind_coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point, h);
//			            v_slip = evaluate_weno_polynomial(upwind_coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point, h);
//			            w_slip = evaluate_weno_polynomial(upwind_coeffs_RHO_W[c], WENO_poly_consts[c], face_quadrature_point, h);


						if(std::fabs(u_slip) > std::fabs(u_slip_max)) {
							u_slip_max = u_slip;
							max_c = c;
							max_g_i = g_i;
						}

						if(std::fabs(v_slip) > std::fabs(v_slip_max))
							v_slip_max = v_slip;
						if(std::fabs(w_slip) > std::fabs(w_slip_max))
							w_slip_max = w_slip;
						double nx = Cell[c].nx(f,q);
						double ny = Cell[c].ny(f,q);
						double nz = Cell[c].nz(f,q);

						grad_rho = evaluate_conservative_gradient(coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point, h);
//						grad_rho = evaluate_conservative_gradient(upwind_coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point, h);
						gradRho = grad_rho(0)*nx + grad_rho(1)*ny + grad_rho(2)*nz;
						if(std::fabs(gradRho) > std::fabs(gradRho_max))
							gradRho_max = gradRho;

						grad_e = evaluate_conservative_gradient(coeffs_E[c], WENO_poly_consts[c], face_quadrature_point, h);
//						grad_e = evaluate_conservative_gradient(upwind_coeffs_E[c], WENO_poly_consts[c], face_quadrature_point, h);
						gradE = grad_e(0)*nx + grad_e(1)*ny + grad_e(2)*nz;

						if(std::fabs(gradE) > std::fabs(gradE_max))
							gradE_max = gradE;

					}

				}

			}
		}
       
    } // End of cell loop 

	global_h_max = Utilities::MPI::max (h_max, MPI_COMM_WORLD);
	global_u_slip_max = Utilities::MPI::max (std::fabs(u_slip_max), MPI_COMM_WORLD);
	global_v_slip_max = Utilities::MPI::max (std::fabs(v_slip_max), MPI_COMM_WORLD);
	global_w_slip_max = Utilities::MPI::max (std::fabs(w_slip_max), MPI_COMM_WORLD);

	double global_gradRho_max = Utilities::MPI::max (std::fabs(gradRho_max), MPI_COMM_WORLD);
	double global_gradE_max = Utilities::MPI::max (std::fabs(gradE_max), MPI_COMM_WORLD);

	unsigned int global_first_order_cell = Utilities::MPI::max (first_order_cell, MPI_COMM_WORLD);

//	pcout<<"h_max_cylinder: "<<global_h_max<<std::endl;
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
		std::ofstream fout_convergence ; 
		fout_convergence.flags( std::ios::dec | std::ios::scientific ) ; 
		fout_convergence.precision(8) ;

    	const std::string filename = "recon_quad.dat";
		fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);
	
	   	fout_convergence << time 
   					<<"\t"<<global_first_order_cell
   					<<"\t"<<global_u_slip_max
   					<<"\t"<<global_v_slip_max
   					<<"\t"<<global_w_slip_max
   					<<"\t"<<global_gradRho_max
   					<<"\t"<<global_gradE_max
   					<< std::endl;
		fout_convergence.close();
	} 

} // End of function 

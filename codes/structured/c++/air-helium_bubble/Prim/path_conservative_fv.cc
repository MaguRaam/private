#include "path_conservative_fv.hh" 


Path_Conservative_FV::Path_Conservative_FV(std::vector<double> (*f)(double, double), AppCtx params):
	init_func(f),
	Params(params),
	FV(params.order),
	riemann(params.p1, params.p2, params.g1, params.g2),
	U(boost::extents[params.N_x+2*N_ph][params.N_y+2*N_ph][Riemann_Solver::n_comp]),
	W(boost::extents[params.N_x+2*N_ph][params.N_y+2*N_ph][Riemann_Solver::n_comp]),
	F(boost::extents[params.N_x+1][params.N_y][Riemann_Solver::n_comp]),
	G(boost::extents[params.N_x][params.N_y+1][Riemann_Solver::n_comp]),
	D(boost::extents[params.N_x+1][params.N_y][Riemann_Solver::n_comp]),
	E(boost::extents[params.N_x][params.N_y+1][Riemann_Solver::n_comp]),
	RHS(boost::extents[params.N_x][params.N_y][Riemann_Solver::n_comp]),
	x(boost::extents[params.N_x+2*N_ph]),
	y(boost::extents[params.N_y+2*N_ph]),
	dt(0.0),
	time(params.initial_time),
	time_step(0),
	rk_stage(1)
{
	
	dx = (Params.x_max - Params.x_min)/static_cast<double>(Params.N_x);
	dy = (Params.y_max - Params.y_min)/static_cast<double>(Params.N_y);
	
	boost::array<array_type::index, 2> bases_1d = {{-N_ph}};

	x.reindex(bases_1d);
	y.reindex(bases_1d);
	
	for (int i = -N_ph; i < Params.N_x + N_ph; i++)
		x[i] = Params.x_min + ((i+1)-0.5)*dx;


	for (int j = -N_ph; j < Params.N_y + N_ph; j++)
		y[j] = Params.y_min + ((j+1)-0.5)*dy;

	int N_gp_face;

	if (Params.order == 2) 
		N_gp_face = 1; 
	else 
		N_gp_face = 2; 
	
	// Arrays to store values of conserved variables at the quadrature points in space and time

	U_L.resize(boost::extents[Params.N_x+2][Params.N_y+2][N_gp_face][Riemann_Solver::n_comp]);
	U_R.resize(boost::extents[Params.N_x+2][Params.N_y+2][N_gp_face][Riemann_Solver::n_comp]);
	U_B.resize(boost::extents[Params.N_x+2][Params.N_y+2][N_gp_face][Riemann_Solver::n_comp]);
	U_T.resize(boost::extents[Params.N_x+2][Params.N_y+2][N_gp_face][Riemann_Solver::n_comp]);
	
	boost::array<array_type::index, 3> bases_3d = {{-N_ph, -N_ph,0}};
	
	U.reindex(bases_3d);
	W.reindex(bases_3d);
	
	boost::array<array_type::index, 4> bases_2d_face = {{-1, -1, 0, 0}};
	
	U_L.reindex(bases_2d_face); U_R.reindex(bases_2d_face); U_B.reindex(bases_2d_face); U_T.reindex(bases_2d_face);
	
	for (int i = -N_ph; i < Params.N_x + N_ph; ++i) {
		for (int j = -N_ph; j < Params.N_y + N_ph; ++j) {
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				U[i][j][c].resize(FV.no_dofs(), 0.0);
				W[i][j][c].resize(FV.no_dofs(), 0.0);
			}
		}
	}
}

void Path_Conservative_FV::initialize() {
	
	std::cout << "Initializing the solution" << std::endl;

	int N_gp = 5;
	std::vector<double> xi(4); std::vector<double> yi(4);

	std::vector<double> ti(N_gp);
	std::vector<double> wi(N_gp);

	QGauss<double> (ti, wi, N_gp);

	std::vector<double> integral(Riemann_Solver::n_comp);
	double xp, yp;

	std::vector<double> Q(Riemann_Solver::n_comp); 
	std::vector<double> V(Riemann_Solver::n_comp);

	// Loop through all the cells

	for (int i = 0; i < Params.N_x; ++i) {
		for (int j = 0; j < Params.N_y; ++j) {

			// Get the vertex coordinates of the cell

			xi[0] = x[i] - dx/2.0; xi[1] = x[i] + dx/2.0;
			xi[2] = x[i] - dx/2.0; xi[3] = x[i] + dx/2.0;

			yi[0] = y[j] - dy/2.0; yi[1] = y[j] - dy/2.0;
			yi[2] = y[j] + dy/2.0; yi[3] = y[j] + dy/2.0;

			// For cell (i,j) initialize with cell averages
			
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) 
				integral[c] = 0.0; 

			for(int l = 0; l < N_gp; l++) {
				for (int m = 0; m < N_gp; m++) {

					isoparametric_transform(ti[l], ti[m], xi, yi, xp, yp);

					V = init_func(xp, yp);

					riemann.primitive_to_conserved(V, Q);
					
					for (int c = 0; c < Riemann_Solver::n_comp; ++c) 
						integral[c] += wi[l]*wi[m]*Q[c];
        		}
        	}
			
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) 
				U[i][j][c][0] = integral[c];
        }
    }

    std::cout << "Done!" << std::endl;
}


void Path_Conservative_FV::compute_primitive_variables() {
	
	std::vector<double> Q(Riemann_Solver::n_comp); 
	std::vector<double> V(Riemann_Solver::n_comp);
	
	// Loop through all the cells

    for (index i = 0; i < Params.N_x; i++) {
        for (index j = 0; j < Params.N_y; j++) {
			
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) 
				Q[c] = U[i][j][c][0];
			
			riemann.conserved_to_primitive(Q, V);
			
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) 
				W[i][j][c][0] = V[c];
		}
	}
}

void Path_Conservative_FV::apply_boundary_conditions() {

	int oned_begin, oned_end, ilhs, irhs;

	// ---------------------- Left boundary ----------------------

	for (index j = -N_ph; j < Params.N_y + N_ph; ++j) {

		oned_begin = 0; oned_end = Params.N_x-1;

		for (index i = 0; i < N_ph; ++i) {

			// Outflow/Transmissive boundary

			if (Params.left_boundary == transmissive) {

				ilhs = oned_begin - i - 1;
				irhs = oned_begin;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[ilhs][j][c][0] = U[irhs][j][c][0];
					W[ilhs][j][c][0] = W[irhs][j][c][0];
					
				}
			}
			
			// Reflective boundary 
			
			if (Params.left_boundary == reflective) {
				ilhs = oned_begin - i - 1;
				irhs = oned_begin + i;

				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[ilhs][j][c][0] = U[irhs][j][c][0];
					W[ilhs][j][c][0] = W[irhs][j][c][0];
				}
				
				U[ilhs][j][1][0] = -U[ilhs][j][1][0];
				W[ilhs][j][1][0] = -W[ilhs][j][1][0];
			}
			
			// Periodic boundary 

			if (Params.left_boundary == periodic) {
				
				ilhs = oned_begin - i - 1;
				irhs = oned_end - i;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[ilhs][j][c][0] = U[irhs][j][c][0];
					W[ilhs][j][c][0] = W[irhs][j][c][0];
					
				}
			}
		}
	}

	// ---------------------- Right boundary ----------------------

	for (index j = -N_ph; j < Params.N_y + N_ph; ++j) {

		oned_begin = 0; oned_end = Params.N_x-1;

		for (index i = 0; i < N_ph; ++i) {

			// Outflow/Transmissive boundary

			if (Params.right_boundary == transmissive) {

				ilhs = oned_end + i + 1;
				irhs = oned_end;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[ilhs][j][c][0] = U[irhs][j][c][0];
					W[ilhs][j][c][0] = W[irhs][j][c][0];
				}
				
			}
			
			// Reflective boundary
			
			if (Params.right_boundary == reflective) {

				ilhs = oned_end + i + 1;
				irhs = oned_end - i;

				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[ilhs][j][c][0] = U[irhs][j][c][0];
					W[ilhs][j][c][0] = W[irhs][j][c][0];
				}
				
				U[ilhs][j][1][0] = -U[ilhs][j][1][0];
				W[ilhs][j][1][0] = -W[ilhs][j][1][0];
			}
			
			// Periodic boundary 

			if (Params.right_boundary == periodic) {

				ilhs = oned_end + i + 1;
				irhs = oned_begin + i;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[ilhs][j][c][0] = U[irhs][j][c][0];
					W[ilhs][j][c][0] = W[irhs][j][c][0];
					
				}
			}
		}
	}

	// ---------------------- Bottom boundary ----------------------

	for (index i = -N_ph; i < Params.N_x + N_ph; ++i) {

		oned_begin = 0; oned_end = Params.N_y-1;

		for (index j = 0; j < N_ph; ++j) {

			// Outflow/Transmissive boundary

			if (Params.bottom_boundary == transmissive) {

				ilhs = oned_begin - j - 1;
				irhs = oned_begin;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[i][ilhs][c][0] = U[i][irhs][c][0];
					W[i][ilhs][c][0] = W[i][irhs][c][0];
				}
			}
			
			// Reflective boundary
			
			if (Params.bottom_boundary == reflective) {
				ilhs = oned_begin - j - 1;
				irhs = oned_begin + j;

				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[i][ilhs][c][0] = U[i][irhs][c][0];
					W[i][ilhs][c][0] = W[i][irhs][c][0];
				}

                                 
                                U[i][ilhs][2][0] = -U[i][ilhs][2][0];
                                W[i][ilhs][2][0] = -W[i][ilhs][2][0];

			}
			
			
			// Periodic boundary 

			if (Params.bottom_boundary == periodic) {
				ilhs = oned_begin - j - 1;
				irhs = oned_end - j;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[i][ilhs][c][0] = U[i][irhs][c][0];
					W[i][ilhs][c][0] = W[i][irhs][c][0];
				}
			}
			
		}
	}

	// ---------------------- Top boundary ----------------------

	for (index i = -N_ph; i < Params.N_x + N_ph; ++i) {

		oned_begin = 0; oned_end = Params.N_y-1;

		for (index j = 0; j < N_ph; ++j) {

			// Outflow/Transmissive boundary

			if (Params.top_boundary == transmissive) {

				ilhs = oned_end + j + 1;
				irhs = oned_end;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[i][ilhs][c][0] = U[i][irhs][c][0];
					W[i][ilhs][c][0] = W[i][irhs][c][0];
				}
			}
			
			// Reflective boundary
			
			if (Params.top_boundary == reflective) {

				ilhs = oned_end + j + 1;
				irhs = oned_end;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[i][ilhs][c][0] = U[i][irhs][c][0];
					W[i][ilhs][c][0] = W[i][irhs][c][0];
				}
				
				U[i][ilhs][2][0] = -U[i][irhs][2][0];
				W[i][ilhs][2][0] = -W[i][irhs][2][0];
			}
			
			if (Params.top_boundary == reflective) {
				
				ilhs = oned_end + j + 1;
				irhs = oned_end - j;

				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[i][ilhs][c][0] = U[i][irhs][c][0];
					W[i][ilhs][c][0] = W[i][irhs][c][0];
				}
				
				U[i][ilhs][2][0] = -U[i][irhs][2][0];
				W[i][ilhs][2][0] = -W[i][irhs][2][0];
			}
			
			// Periodic boundary 

			if (Params.top_boundary == periodic) {
				ilhs = oned_end + j + 1;
				irhs = oned_begin + j;
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					U[i][ilhs][c][0] = U[i][irhs][c][0];
					W[i][ilhs][c][0] = W[i][irhs][c][0];
				}

			}
		}
	}
}

void Path_Conservative_FV::compute_rhs(double t) {

	t = t*1.0; 
	
	int N_gp_face;

	if (Params.order == 2)
		N_gp_face = 1; 
	else 
		N_gp_face = 2;
	
	double s, s_max = 0.0; 


	const double r1_dx = 1./dx;
	const double r1_dy = 1./dx;
	
	// Flux Vectors

    std::vector<double> cons_F(Riemann_Solver::n_comp); 
	std::vector<double> ncons_F(Riemann_Solver::n_comp);
    
	std::vector<double> QL(Riemann_Solver::n_comp); 
	std::vector<double> QR(Riemann_Solver::n_comp);
	std::vector<double> QB(Riemann_Solver::n_comp); 
	std::vector<double> QT(Riemann_Solver::n_comp);
	
	std::vector<double> VL(Riemann_Solver::n_comp); 
	std::vector<double> VR(Riemann_Solver::n_comp);
	std::vector<double> VT(Riemann_Solver::n_comp);
	std::vector<double> VB(Riemann_Solver::n_comp);

	// Vectors for Stencils

	std::vector<double> u_xloc(5); std::vector<double> u_yloc(5);
	std::vector<double> u_xyloc(5);

	// Quadrature points in space and time

	std::vector<double> face_q_points(N_gp_face);
	std::vector<double> face_q_weights(N_gp_face);
	QGauss<double> (face_q_points, face_q_weights, N_gp_face);
	
	apply_boundary_conditions();
	
	if (Params.reconstruct_primitive_variables)
		compute_primitive_variables(); 

	//----------------------- Do WENO reconstruction -----------------------

	for (int i = -1; i < Params.N_x+1; i++) {
		for (int j = -1; j < Params.N_y+1; j++) {
			
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
			
				// Stencil in x-direction

				u_xloc[0] = U[i-2][j][c][0]; u_xloc[1] = U[i-1][j][c][0]; u_xloc[2] = U[i][j][c][0];
				u_xloc[3] = U[i+1][j][c][0]; u_xloc[4] = U[i+2][j][c][0];
				
				// Stencil in y-direction

				u_yloc[0] = U[i][j-2][c][0]; u_yloc[1] = U[i][j-1][c][0]; u_yloc[2] = U[i][j][c][0];
				u_yloc[3] = U[i][j+1][c][0]; u_yloc[4] = U[i][j+2][c][0];
				
				// Stencil in xy direction

				u_xyloc[0] = U[i][j][c][0];     u_xyloc[1] = U[i+1][j+1][c][0]; u_xyloc[2] = U[i+1][j-1][c][0];
				u_xyloc[3] = U[i-1][j+1][c][0]; u_xyloc[4] = U[i-1][j-1][c][0];
				
				weno_2d(u_xloc, u_yloc, u_xyloc, U[i][j][c], Params.order);
			}
			
			if (Params.reconstruct_primitive_variables) {
			
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
			
					// Stencil in x-direction

					u_xloc[0] = W[i-2][j][c][0]; u_xloc[1] = W[i-1][j][c][0]; u_xloc[2] = W[i][j][c][0];
					u_xloc[3] = W[i+1][j][c][0]; u_xloc[4] = W[i+2][j][c][0];
					
					// Stencil in y-direction

					u_yloc[0] = W[i][j-2][c][0]; u_yloc[1] = W[i][j-1][c][0]; u_yloc[2] = W[i][j][c][0];
					u_yloc[3] = W[i][j+1][c][0]; u_yloc[4] = W[i][j+2][c][0];
					
					// Stencil in xy direction

					u_xyloc[0] = W[i][j][c][0];     u_xyloc[1] = W[i+1][j+1][c][0]; u_xyloc[2] = W[i+1][j-1][c][0];
					u_xyloc[3] = W[i-1][j+1][c][0]; u_xyloc[4] = W[i-1][j-1][c][0];
					
					weno_2d(u_xloc, u_yloc, u_xyloc, W[i][j][c], Params.order);
				}
			}
			
			// Find the values of conserved variables at quadrature points
			
			if (Params.reconstruct_primitive_variables) {
			
				for (int f = 0; f < N_gp_face; ++f) {
					
					for (int c = 0; c < Riemann_Solver::n_comp; ++c) {

						VL[c] = FV.evaluate(W[i][j][c], -0.5, face_q_points[f]);
						VR[c] = FV.evaluate(W[i][j][c],  0.5, face_q_points[f]);
						VB[c] = FV.evaluate(W[i][j][c], face_q_points[f], -0.5);
						VT[c] = FV.evaluate(W[i][j][c], face_q_points[f],  0.5);

					}
					
					riemann.primitive_to_conserved(VL, QL);
					riemann.primitive_to_conserved(VR, QR);
					riemann.primitive_to_conserved(VB, QB);
					riemann.primitive_to_conserved(VT, QT);
					
					for (int c = 0; c < Riemann_Solver::n_comp; ++c) {

						U_L[i][j][f][c] = QL[c];
						U_R[i][j][f][c] = QR[c];
						U_B[i][j][f][c] = QB[c];
						U_T[i][j][f][c] = QT[c];

					}
				}
			}
			
			else {

				for (int f = 0; f < N_gp_face; ++f) {
					for (int c = 0; c < Riemann_Solver::n_comp; ++c) {

						U_L[i][j][f][c] = FV.evaluate(U[i][j][c], -0.5, face_q_points[f]);
						U_R[i][j][f][c] = FV.evaluate(U[i][j][c],  0.5, face_q_points[f]);
						U_B[i][j][f][c] = FV.evaluate(U[i][j][c], face_q_points[f], -0.5);
						U_T[i][j][f][c] = FV.evaluate(U[i][j][c], face_q_points[f],  0.5);

					}
				}
			}
		}
	}
	
	// Find the x-component of the flux

	for (int i = 0; i < Params.N_x + 1; ++i) {
		for (int j = 0; j < Params.N_y; ++j) {
			
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				F[i][j][c] = 0.0; 
				D[i][j][c] = 0.0; 
			}
			
			for (int f = 0; f < N_gp_face; ++f) {
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					QL[c] = U_R[i-1][j][f][c];   QR[c] = U_L[i][j][f][c];
				}

				s = riemann.llf_riemann_solver(QL, QR, 1.0, 0.0, 0.5*(x[i-1] + x[i]), y[j], cons_F, ncons_F);
				
				if (s > s_max)
					s_max = s; 

				// Get average values
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					F[i][j][c] += face_q_weights[f]*cons_F[c];
					D[i][j][c] += face_q_weights[f]*ncons_F[c];
				}
			}
		
			
		}
	}
	
	// Find the y-component of the flux

	for (int i = 0; i < Params.N_x; ++i) {
		for (int j = 0; j < Params.N_y + 1; ++j) {
			
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				G[i][j][c] = 0.0; 
				E[i][j][c] = 0.0; 
			}
			
			for (int f = 0; f < N_gp_face; ++f) {
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					QL[c] = U_T[i][j-1][f][c]; QR[c] = U_B[i][j][f][c];
				}

				s = riemann.llf_riemann_solver(QL, QR, 0.0, 1.0, x[i], 0.5*(y[j-1] + y[j]), cons_F, ncons_F);

				if (s > s_max)
					s_max = s; 
				
				// Get average values
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
					G[i][j][c] += face_q_weights[f]*cons_F[c];
					E[i][j][c] += face_q_weights[f]*ncons_F[c];
				}
			}
		}
	}
	
	// Find RHS

	for (int i = 0; i < Params.N_x; ++i) {
		for (int j = 0; j < Params.N_y; ++j) {
			
			for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				RHS[i][j][c] =  -r1_dx*( F[i+1][j][c] - F[i][j][c] + 0.5*(D[i+1][j][c] + D[i][j][c]) ) 
								-r1_dy*( G[i][j+1][c] - G[i][j][c] + 0.5*(E[i][j+1][c] + E[i][j][c]) );
			}
		}
	}
	
	// Add the source term to the RHS
	
	int N_gp = 4;
	
	double x_gp[] = {-0.28867513459481287, -0.28867513459481287,  0.28867513459481287, 0.28867513459481287};
	double y_gp[] = {-0.28867513459481287,  0.28867513459481287, -0.28867513459481287, 0.28867513459481287};
	double w_gp[] = {0.25,  0.25, 0.25, 0.25};
	
	std::vector<double> grad_Qx(Riemann_Solver::n_comp);
	std::vector<double> grad_Qy(Riemann_Solver::n_comp);
	std::vector<double> Q(Riemann_Solver::n_comp);
	
	double* grad; grad = new double[2];
	
	for (int i = 0; i < Params.N_x; ++i) {
		for (int j = 0; j < Params.N_y; ++j) {
			
			for (int q = 0; q < N_gp; ++q) {
			
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				
					Q[c] = FV.evaluate(U[i][j][c], x_gp[q], y_gp[q]);
					
					FV.evaluate_grad(U[i][j][c], x_gp[q], y_gp[q], grad);
					
					grad_Qx[c] = grad[0]; grad_Qy[c] = grad[1];
					
				}
				
				riemann.non_conservative_flux(Q, grad_Qx, grad_Qy, ncons_F);
				
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) 
					RHS[i][j][c] += -r1_dx*w_gp[q]*ncons_F[c];
			}
		}
	}
	
	delete[] grad; 
	
	// Find the time step size (only if rk_stage = 1)
	
	if (rk_stage == 1) {
		dt = Params.CFL*dx/s_max;
		
		// If time step exceeds the final time, reduce it accordingly 
	
		if((time + dt)>Params.final_time)
			dt = Params.final_time - time;
	}
}


void Path_Conservative_FV::solve_ssprk22() {
	
	std::cout << "Solving using SSPRK (2,2) method" << std::endl;
	
	boost::multi_array<double, 3> U_old(boost::extents[Params.N_x][Params.N_y][Riemann_Solver::n_comp]);
	
	while (time < Params.final_time) {
		
		printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.final_time);
		
		if (Params.write_interval != 0) {
			if (time_step % Params.write_interval == 0) {
				compute_primitive_variables();
				plot_tecplot(time_step);
			}
		}
		
        //  Stage 1
		
		rk_stage = 1;
        compute_rhs(time);
		
		for (int i = 0; i < Params.N_x; ++i) {
			for (int j = 0; j < Params.N_y; ++j) {
			
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				
					U_old[i][j][c] =   U[i][j][c][0];
					
					U[i][j][c][0] +=   dt*RHS[i][j][c];
					
				}
			}
		}
			
		//  Stage 2
		
		rk_stage = 2;
		compute_rhs(time + dt);
		
		
		for (int i = 0; i < Params.N_x; ++i) {
			for (int j = 0; j < Params.N_y; ++j) {
			
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				
					U[i][j][c][0] =   0.5*( U_old[i][j][c] +   U[i][j][c][0] + dt*RHS[i][j][c]);
					
				}
			}
		}

		time += dt;
		time_step++;
	}
	
	printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.final_time);
	
	compute_primitive_variables();
	plot_tecplot(time_step);
	
}

void Path_Conservative_FV::solve_ssprk33() {
	
	std::cout << "Solving using SSPRK (3,3) method" << std::endl;
	
	const double r2_3 = 2./3.; 
	
	boost::multi_array<double, 3> U_old(boost::extents[Params.N_x][Params.N_y][Riemann_Solver::n_comp]);
	
	while (time < Params.final_time) {
		
		printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.final_time);
		
		if (Params.write_interval != 0) {
			if (time_step % Params.write_interval == 0) {
				compute_primitive_variables();
				plot_tecplot(time_step);
			}
		}
		
        //  Stage 1
		
		rk_stage = 1;
        compute_rhs(time);
		
		for (int i = 0; i < Params.N_x; ++i) {
			for (int j = 0; j < Params.N_y; ++j) {
			
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				
					U_old[i][j][c] =   U[i][j][c][0];
					
					U[i][j][c][0] +=   dt*RHS[i][j][c];
					
				}
			}
		}
			
		//  Stage 2
		
		rk_stage = 2;
		compute_rhs(time + dt);
		
		
		for (int i = 0; i < Params.N_x; ++i) {
			for (int j = 0; j < Params.N_y; ++j) {
			
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				
					U[i][j][c][0] =   0.25*(3.0*U_old[i][j][c] +   U[i][j][c][0] + dt*RHS[i][j][c]);
					
				}
			}
		}

		//  Stage 2
		
		rk_stage = 3;
		compute_rhs(time + 0.5*dt);
		
		
		for (int i = 0; i < Params.N_x; ++i) {
			for (int j = 0; j < Params.N_y; ++j) {
			
				for (int c = 0; c < Riemann_Solver::n_comp; ++c) {
				
					U[i][j][c][0] =   r2_3*( 0.5*U_old[i][j][c] +   U[i][j][c][0] + dt*RHS[i][j][c]);
					
				}
			}
		}
		
		time += dt;
		time_step++;
	}
	
	printf ("time = %4.3e, dt = %4.3e, final time = %4.3e\n", time, dt, Params.final_time);
	
	compute_primitive_variables();
	plot_tecplot(time_step);
	
}



void Path_Conservative_FV::plot_tecplot(int i, const int digits)  {

    std::ofstream tecplot;
    const std::string filename = "../plots/plot_" + Utilities::int_to_string (i, digits) + ".dat";
    tecplot.open (filename);
    tecplot.flags( std::ios::dec | std::ios::scientific );
	tecplot.precision(6);

	Assert(tecplot.is_open(), ErrIO(filename));

	// Plot everything

    tecplot << "TITLE = \"PathConservative-GAMMA_MODEL\" " << std::endl
	    << "VARIABLES = \"x\", \"y\", \"rho\", \"v_x\", \"v_y\", \"p\", \"phi\" " << std::endl;
    tecplot << "Zone I = " << Params.N_y << " J = " << Params.N_x << std::endl;

    for (int i=0; i < Params.N_x; i++) {
	    for (int j=0; j < Params.N_y; j++) {
            tecplot << x[i] << "\t" << y[j] << "\t"
            		<< W[i][j][0][0] << "\t"
			        << W[i][j][1][0] << "\t"
					<< W[i][j][2][0] << "\t"
					<< W[i][j][3][0] << "\t"
					<< W[i][j][4][0] << std::endl;
	    }
    }

    tecplot.close();
}

void Path_Conservative_FV::run() {
	
	initialize(); 
	
	solve_ssprk33();
	
	
}

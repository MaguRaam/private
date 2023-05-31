#include "../include/Weno32.h"

// Initialize the solution 

void Weno3_2D::initialize() {

    TimerOutput::Scope t(computing_timer, "initialize");

    pcout << "Initializing the solution" << std::endl;

    unsigned int N_gp = 4;               // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);

    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

    double V0;
	double h_old,h_min_local; 
	h_min_local = 1e6; 

	unsigned int c, g_i;
    

	Vector<double> U(4); Vector<double> W(4);

	Point<2> q_point;

	for (unsigned int c = 0; cell != endc; ++cell, ++c) {

		if (cell->is_locally_owned()) {
	
			cell->get_dof_indices(local_dof_indices);
			fv_values.reinit(cell);
			g_i = local_dof_indices[0];				


			V0 = 0.0; 
			for (unsigned int i=0; i<fv_values.n_quadrature_points; ++i)
				V0 += fv_values.JxW (i);
	
			h_old = std::sqrt(V0);

			if (h_old < h_min_local) {
				h_min_local = h_old; 
			}
	
			double rho = 0.0, rho_u = 0.0, rho_v = 0.0, e = 0.0;

			for (unsigned int i = 0; i < N_gp*N_gp; i++) {

				q_point = fv_values.quadrature_point(i);

				W = initial_condition(q_point); 
				U = primitive_to_conserved(W); 
	
				rho   += (1./V0)*fv_values.JxW (i)*U(0);
				rho_u += (1./V0)*fv_values.JxW (i)*U(1);
				rho_v += (1./V0)*fv_values.JxW (i)*U(2);
				e     += (1./V0)*fv_values.JxW (i)*U(3);
			}

			local_RHO(g_i)   = rho;
			local_RHO_U(g_i) = rho_u;
			local_RHO_V(g_i) = rho_v;
			local_E(g_i)     = e;

		}

	}
	
	

	local_RHO.compress(VectorOperation::insert);
	local_RHO_U.compress(VectorOperation::insert);
	local_RHO_V.compress(VectorOperation::insert);
	local_E.compress(VectorOperation::insert);

	RHO = local_RHO;
	RHO_U = local_RHO_U;
	RHO_V = local_RHO_V;
	E = local_E;

	h_min = Utilities::MPI::min (h_min_local, MPI_COMM_WORLD);
    
    pcout << "h_min: " <<h_min<< std::endl;  

} 

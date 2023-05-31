#include "../include/Weno432.h"

void Weno4_3D::L_norm() {

	double local_linfty_norm_RHO = 0.0, local_linfty_norm_RHO_U = 0.0, local_linfty_norm_RHO_V = 0.0, local_linfty_norm_RHO_W = 0.0, local_linfty_norm_E = 0.0;
	double linfty_norm_RHO, linfty_norm_RHO_U, linfty_norm_RHO_V, linfty_norm_RHO_W, linfty_norm_E;

	double local_l2_norm_RHO = 0.0, local_l2_norm_RHO_U = 0.0, local_l2_norm_RHO_V = 0.0, local_l2_norm_RHO_W = 0.0, local_l2_norm_E = 0.0;
	double l2_norm_RHO, l2_norm_RHO_U, l2_norm_RHO_V, l2_norm_RHO_W, l2_norm_E;

	Vector<double> local_difference_RHO(n_locally_cells), local_difference_RHO_U(n_locally_cells), local_difference_RHO_V(n_locally_cells), local_difference_RHO_W(n_locally_cells), local_difference_E(n_locally_cells);

	unsigned int n_cells = triangulation.n_active_cells();

	// Get iterators for cell
		
    Vector<double> U(5), U_exact(5); 

	double h, max_error = 0, volume = 0.0;

	unsigned int g_i, g_i_max;

	Point<3> error_location, quadrature_point;

	// Get iterators for cell
		
    DoFHandler<3>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();

	QGauss<3> quadrature_formula(3);

	const UpdateFlags update_flags =   update_values | update_quadrature_points | update_JxW_values;

	FEValues<3> fv_values (mapping, fv, quadrature_formula, update_flags);

	const unsigned int n_q_points	    	    = fv_values.n_quadrature_points;
	std::vector<Vector<double> >	     exact_solution(n_q_points, Vector<double>(5));
		double error_u, error_u_max = 0, error_v, error_v_max = 0, error_w, error_w_max = 0, error_rho, error_rho_max = 0, error_e, error_e_max = 0;		  		 		
	for (unsigned int c = 0; c < n_locally_cells; ++c) {

		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];
		h = Cell[c].h();
		fv_values.reinit(cell);

		true_solution.soln_list(fv_values.get_quadrature_points(),exact_solution,time);


		for (unsigned int i = 0; i < n_q_points; ++i) {
			quadrature_point = fv_values.quadrature_point(i);
/*
   	    	U(0) = evaluate_weno_polynomial(upwind_coeffs_RHO[c], WENO_poly_consts[c], quadrature_point, h);
   	    	U(1) = evaluate_weno_polynomial(upwind_coeffs_RHO_U[c], WENO_poly_consts[c], quadrature_point, h);
   	    	U(2) = evaluate_weno_polynomial(upwind_coeffs_RHO_V[c], WENO_poly_consts[c], quadrature_point, h);
   	    	U(3) = evaluate_weno_polynomial(upwind_coeffs_RHO_W[c], WENO_poly_consts[c], quadrature_point, h);
   	    	U(4) = evaluate_weno_polynomial(upwind_coeffs_E[c], WENO_poly_consts[c], quadrature_point, h);
*/
   	    	U(0) = evaluate_weno_polynomial(coeffs_RHO[c], WENO_poly_consts[c], quadrature_point, h);
   	    	U(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], quadrature_point, h);
   	    	U(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], quadrature_point, h);
   	    	U(3) = evaluate_weno_polynomial(coeffs_RHO_W[c], WENO_poly_consts[c], quadrature_point, h);
   	    	U(4) = evaluate_weno_polynomial(coeffs_E[c], WENO_poly_consts[c], quadrature_point, h);

			claw.primitive_to_conserved(exact_solution[i], U_exact); 
			error_rho = U_exact(0) - U(0);
			error_u = U_exact(1) - U(1);
			error_v = U_exact(2) - U(2);
			error_w = U_exact(3) - U(3);
			error_e = U_exact(4) - U(4);

			if(std::fabs(local_difference_RHO(c)) < std::fabs(error_rho))
				local_difference_RHO(c) = error_rho;
			if(std::fabs(local_difference_RHO_U(c)) < std::fabs(error_u))
				local_difference_RHO_U(c) = error_u;
			if(std::fabs(local_difference_RHO_V(c)) < std::fabs(error_v))
				local_difference_RHO_V(c) = error_v;
			if(std::fabs(local_difference_RHO_W(c)) < std::fabs(error_w))
				local_difference_RHO_W(c) = error_w;
			if(std::fabs(local_difference_E(c)) < std::fabs(error_e))
				local_difference_E(c) = error_e;

//			if(std::fabs(error_u) > std::fabs(error_u_max) ) { error_u_max = error_u; g_i_max = g_i;}

		}

		local_error_RHO(g_i) = local_difference_RHO(c);	
		local_error_RHO_U(g_i) = local_difference_RHO_U(c);	
		local_error_RHO_V(g_i) = local_difference_RHO_V(c);	
		local_error_RHO_W(g_i) = local_difference_RHO_W(c);	
		local_error_E(g_i) = local_difference_E(c);	

	}

//	pcout<<"g_i: "<<g_i_max<<"\nerror: "<<error_u_max<<"\nwe coeff: "<<coeffs_RHO_U[g_i_max]<<"\nup coeff: "<<upwind_coeffs_RHO_U[g_i_max]<<std::endl;

	local_error_RHO.compress(VectorOperation::insert);
	local_error_RHO_U.compress(VectorOperation::insert);
	local_error_RHO_V.compress(VectorOperation::insert);
	local_error_RHO_W.compress(VectorOperation::insert);
	local_error_E.compress(VectorOperation::insert);

	error_RHO = local_error_RHO;
	error_RHO_U = local_error_RHO_U;
	error_RHO_V = local_error_RHO_V;
	error_RHO_W = local_error_RHO_W;
	error_E = local_error_E;

//	std::cout<<"max_error: "<<max_error<<"\terror_location: "<<error_location<<std::endl;
	local_linfty_norm_RHO = local_difference_RHO.linfty_norm();
	local_linfty_norm_RHO_U = local_difference_RHO_U.linfty_norm();
	local_linfty_norm_RHO_V = local_difference_RHO_V.linfty_norm();
	local_linfty_norm_RHO_W = local_difference_RHO_W.linfty_norm();
	local_linfty_norm_E = local_difference_E.linfty_norm();

	local_l2_norm_RHO = local_difference_RHO.l2_norm();
	local_l2_norm_RHO_U = local_difference_RHO_U.l2_norm();
	local_l2_norm_RHO_V = local_difference_RHO_V.l2_norm();
	local_l2_norm_RHO_W = local_difference_RHO_W.l2_norm();
	local_l2_norm_E = local_difference_E.l2_norm();

	linfty_norm_RHO = Utilities::MPI::max (local_linfty_norm_RHO, MPI_COMM_WORLD);
	linfty_norm_RHO_U = Utilities::MPI::max (local_linfty_norm_RHO_U, MPI_COMM_WORLD);
	linfty_norm_RHO_V = Utilities::MPI::max (local_linfty_norm_RHO_V, MPI_COMM_WORLD);
	linfty_norm_RHO_W = Utilities::MPI::max (local_linfty_norm_RHO_W, MPI_COMM_WORLD);
	linfty_norm_E = Utilities::MPI::max (local_linfty_norm_E, MPI_COMM_WORLD);

	l2_norm_RHO = std::sqrt ( ( Utilities::MPI::sum (local_l2_norm_RHO * local_l2_norm_RHO, MPI_COMM_WORLD) )/n_cells);
	l2_norm_RHO_U = std::sqrt ( (Utilities::MPI::sum (local_l2_norm_RHO_U * local_l2_norm_RHO_U, MPI_COMM_WORLD) )/n_cells);
	l2_norm_RHO_V = std::sqrt ( (Utilities::MPI::sum (local_l2_norm_RHO_V * local_l2_norm_RHO_V, MPI_COMM_WORLD))/n_cells);
	l2_norm_RHO_W = std::sqrt ( (Utilities::MPI::sum (local_l2_norm_RHO_W * local_l2_norm_RHO_W, MPI_COMM_WORLD))/n_cells);
	l2_norm_E = std::sqrt ( (Utilities::MPI::sum (local_l2_norm_E * local_l2_norm_E, MPI_COMM_WORLD))/n_cells);

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
		std::ofstream fout_convergence ; 
		fout_convergence.flags( std::ios::dec | std::ios::scientific ) ; 
		fout_convergence.precision(8) ;

    	const std::string filename = "norm_" + Utilities::int_to_string(log10(weight_ls),1) + ".dat";
		fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);
	
	   	fout_convergence << n_cells 
					<<"\t"<<weight_ls
   					<<"\t"<<linfty_norm_RHO
   					<<"\t"<<l2_norm_RHO
   					<<"\t"<<linfty_norm_RHO_U
   					<<"\t"<<l2_norm_RHO_U
   					<<"\t"<<linfty_norm_RHO_V
   					<<"\t"<<l2_norm_RHO_V
   					<<"\t"<<linfty_norm_RHO_W
   					<<"\t"<<l2_norm_RHO_W
   					<<"\t"<<linfty_norm_E
   					<<"\t"<<l2_norm_E
   					<< std::endl;
		fout_convergence.close();
	}
}

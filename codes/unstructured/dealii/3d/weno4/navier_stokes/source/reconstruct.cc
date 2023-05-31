#include "../include/Weno432.h"


// Perform the actual reconstruction 

void Weno4_3D::reconstruct() {

    unsigned int n_faces = GeometryInfo<3>::faces_per_cell, global = 1001;
    Point<3> face_quadrature_point;

//	double gradRho, gradRho_max = 0, rho_slip_max = 0.0, u_slip_max = 0.0, rho_slip, u_slip;

    unsigned int no_stencils = 8;
    unsigned int p = 2; 
    double epsilon = 1.0e-12; 
    double h; // Measure of cell size 
    
    double fourth_order_wt = 0.86;
	double third_order_wt = 0.02; 
	double second_order_wt = 0.0; 
//	double second_order_wt = 0.0125; 
	double sum_gamma;

	Vector<double> gamma(no_stencils);
	
	gamma(0) = fourth_order_wt; 
	gamma(1) = third_order_wt; 
	gamma(2) = third_order_wt;
	gamma(3) = third_order_wt;
	gamma(4) = third_order_wt;
	gamma(5) = third_order_wt;
	gamma(6) = third_order_wt;
	gamma(7) = third_order_wt;
/*
	gamma(8) = second_order_wt;
	gamma(9) = second_order_wt;
	gamma(10) = second_order_wt;
	gamma(11) = second_order_wt;
	gamma(12) = second_order_wt;
	gamma(13) = second_order_wt;
	gamma(14) = second_order_wt;
	gamma(15) = second_order_wt;
*/
	Vector<double> U(5), W(5); 

    // Variables for reconstruction of RHO

	double rho0; 
    
	// Fourth order stencil
	Vector<double> d_rho_4; 
	Vector<double> b_rho_4; 
	Vector<double> rho_coeff_4(19);

	// Third order centered stencil // 
	Vector<double> d_rho_3; 
	Vector<double> b_rho_3;
	Vector<double> rho_coeff_3(9);

	// Third order directional stencil // 
	std::vector<Vector<double> > b_rho_3_d(6);         
	std::vector<Vector<double> > d_rho_3_d(6);
	std::vector<Vector<double> > rho_coeff_3_d(6, Vector<double>(9) );

	std::vector<Vector<double> > b_rho_2_d(8, Vector<double>(3));
	std::vector<Vector<double> > rho_coeff_2_d(8, Vector<double>(3));

    Vector<double> IS_RHO(no_stencils); Vector<double> w_RHO(no_stencils);  double sum_RHO;

    // Variables for reconstruction of RHO_U

	double rho_u0; 

	// Fourth order stencil
	Vector<double> d_rho_u_4; 
	Vector<double> b_rho_u_4; 
	Vector<double> rho_u_coeff_4(19);

	// Third order centered stencil // 
	Vector<double> d_rho_u_3; 
	Vector<double> b_rho_u_3;
	Vector<double> rho_u_coeff_3(9);

	// Third order directional stencil // 
	std::vector<Vector<double> > b_rho_u_3_d(6);         
	std::vector<Vector<double> > d_rho_u_3_d(6);
	std::vector<Vector<double> > rho_u_coeff_3_d(6, Vector<double>(9) );

	/* Second order directional stencil */ 
	std::vector<Vector<double> > b_rho_u_2_d(8, Vector<double>(3));
	std::vector<Vector<double> > rho_u_coeff_2_d(8, Vector<double>(3));

    Vector<double> IS_RHO_U(no_stencils); Vector<double> w_RHO_U(no_stencils);  double sum_RHO_U;

    // Variables for reconstruction of RHO_V

	double rho_v0; 

	// Fourth order stencil
	Vector<double> d_rho_v_4; 
	Vector<double> b_rho_v_4; 
	Vector<double> rho_v_coeff_4(19);

	// Third order centered stencil // 
	Vector<double> d_rho_v_3; 
	Vector<double> b_rho_v_3;
	Vector<double> rho_v_coeff_3(9);

	// Third order directional stencil // 
	std::vector<Vector<double> > b_rho_v_3_d(6);         
	std::vector<Vector<double> > d_rho_v_3_d(6);
	std::vector<Vector<double> > rho_v_coeff_3_d(6, Vector<double>(9) );

	/* Second order directional stencil */ 
	std::vector<Vector<double> > b_rho_v_2_d(8, Vector<double>(3));
	std::vector<Vector<double> > rho_v_coeff_2_d(8, Vector<double>(3));

    Vector<double> IS_RHO_V(no_stencils); Vector<double> w_RHO_V(no_stencils);  double sum_RHO_V;

    // Variables for reconstruction of RHO_W

	double rho_w0; 

	// Fourth order stencil
	Vector<double> d_rho_w_4; 
	Vector<double> b_rho_w_4; 
	Vector<double> rho_w_coeff_4(19);

	// Third order centered stencil // 
	Vector<double> d_rho_w_3; 
	Vector<double> b_rho_w_3;
	Vector<double> rho_w_coeff_3(9);

	// Third order directional stencil // 
	std::vector<Vector<double> > b_rho_w_3_d(6);         
	std::vector<Vector<double> > d_rho_w_3_d(6);
	std::vector<Vector<double> > rho_w_coeff_3_d(6, Vector<double>(9) );

	/* Second order directional stencil */ 
	std::vector<Vector<double> > b_rho_w_2_d(8, Vector<double>(3));
	std::vector<Vector<double> > rho_w_coeff_2_d(8, Vector<double>(3));

    Vector<double> IS_RHO_W(no_stencils); Vector<double> w_RHO_W(no_stencils);  double sum_RHO_W;

	double e0; 
	// Fourth order stencil
	Vector<double> d_e_4; 
	Vector<double> b_e_4; 
	Vector<double> e_coeff_4(19);

	// Third order centered stencil // 
	Vector<double> d_e_3; 
	Vector<double> b_e_3;
	Vector<double> e_coeff_3(9);

	// Third order directional stencil // 
	std::vector<Vector<double> > b_e_3_d(6);         
	std::vector<Vector<double> > d_e_3_d(6);
	std::vector<Vector<double> > e_coeff_3_d(6, Vector<double>(9) );

	/* Second order directional stencil */ 
	std::vector<Vector<double> > b_e_2_d(8, Vector<double>(3));
	std::vector<Vector<double> > e_coeff_2_d(8, Vector<double>(3));

    Vector<double> IS_E(no_stencils); Vector<double> w_E(no_stencils);  double sum_E;

    // Iterate over all the cells 

    DoFHandler<3>::active_cell_iterator cell, neighbor;

	std::vector<bool> n_n_face(6, false);
	std::vector<unsigned int> n_n_face_index(6), neighbor_index(6, numbers::invalid_unsigned_int);  
	unsigned int index_least;
    unsigned int ROWS, ROWS_slip, index, g_i, fff, n_boundary_quad_points = 2, n_boundary_quad_points_edge = 1;
	bool negative = false;
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];
		is_1st_order[c] = false;

      	for (unsigned int f = 0; f < n_faces; f++){	
			if (cell_neighbor_neighbor_index[c][f].size() > 0) {
				n_n_face[f] = true;
				n_n_face_index[f] = cell_neighbor_neighbor_index[c][f][0];
			}
			else {
				n_n_face[f] = false;
				n_n_face_index[f] = numbers::invalid_unsigned_int;
			}

			if(!cell->face(f)->at_boundary()) {
				neighbor = cell->neighbor(f);
			    neighbor->get_dof_indices(local_neighbor_dof_indices);
				neighbor_index[f] = local_neighbor_dof_indices[0];
			}
			else
				neighbor_index[f] = numbers::invalid_unsigned_int;
		}

        rho0   =   RHO(g_i);
    	rho_u0 =   RHO_U(g_i);
    	rho_v0 =   RHO_V(g_i);
    	rho_w0 =   RHO_W(g_i);
    	e0     =   E(g_i);

		h = Cell[c].h();
        
        if ( !(cell->at_boundary()) ) {
//			pcout<<"interior g_i: "<<g_i<<"\tcenter: "<<cell->center()<<std::endl;
            // =====================================================================
            // r = 4 stencil 
            // =====================================================================

			d_rho_4.reinit(6);			d_rho_u_4.reinit(6);
			d_rho_v_4.reinit(6);		d_rho_w_4.reinit(6);
			d_e_4.reinit(6);

           	for (unsigned int f = 0; f < n_faces; f++) {			
            	
				d_rho_4(f)  = (RHO(neighbor_index[f]) - rho0);
				d_rho_u_4(f)  = (RHO_U(neighbor_index[f]) - rho_u0);
				d_rho_v_4(f)  = (RHO_V(neighbor_index[f]) - rho_v0);
				d_rho_w_4(f)  = (RHO_W(neighbor_index[f]) - rho_w0);
				d_e_4(f)  = (E(neighbor_index[f]) - e0);
			} 

			ROWS = cell_diagonal_neighbor_index[c].size(); 
			
			// neighbor of neighbors 

	      	for (unsigned int f = 0; f < n_faces; f++)
				if (n_n_face[f])
					ROWS++; 
				
			index = 0; 
            
            // Least Squares Part  
            
            b_rho_4.reinit(ROWS);	b_rho_u_4.reinit(ROWS);	b_rho_v_4.reinit(ROWS);
			b_rho_w_4.reinit(ROWS);	b_e_4.reinit(ROWS);
			
			// vertex neighbors
			
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];
				b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
				b_rho_u_4(index)   = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0; 
				b_rho_v_4(index)   = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0; 
				b_rho_w_4(index)   = RHO_W(local_neighbor_dof_indices[0] ) - rho_w0; 
				b_e_4(index)   = E(local_neighbor_dof_indices[0] ) - e0; 
				index++; 
			}
			
			// neighbors of neighbors 

	      	for (unsigned int f = 0; f < n_faces; f++) {
				if (n_n_face[f]){
					local_neighbor_dof_indices[0] = n_n_face_index[f];
					b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0;
					b_rho_u_4(index)   = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_4(index)   = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_rho_w_4(index)   = RHO_W(local_neighbor_dof_indices[0] ) - rho_w0;
					b_e_4(index)   = E(local_neighbor_dof_indices[0] ) - e0;
					index++;
				}
            }

//			////pcout<<"cls_R4 interior: "<<c<<std::endl;
//			////pcout<<"b_rho: "<<b_rho_u_4<<std::endl;
//			////pcout<<"d_rho: "<<d_rho_u_4<<std::endl;
//           	if(ROWS == 19) std::cout<<"c: "<<c<<std::endl;
			CLS_R4[c].solve(b_rho_4, d_rho_4, rho_coeff_4); 
			CLS_R4[c].solve(b_rho_u_4, d_rho_u_4, rho_u_coeff_4); 
			CLS_R4[c].solve(b_rho_v_4, d_rho_v_4, rho_v_coeff_4); 
			CLS_R4[c].solve(b_rho_w_4, d_rho_w_4, rho_w_coeff_4); 
			CLS_R4[c].solve(b_e_4, d_e_4, e_coeff_4); 

//			////pcout<<"cls_R4 interior"<<std::endl;

//			if(c == 1) std::cout<<"b_e_4: "<<b_e_4<<std::endl<<"d_e_4: "<<d_e_4<<std::endl<<"e_coeff_4: "<<e_coeff_4<<std::endl;
//			if(c == 1) std::cout<<"b_e_4: "<<b_e_4(0)<<"\t"<<b_e_4(1)<<"\t"<<b_e_4(2)<<"\t"<<b_e_4(3)<<std::endl;
//			////pcout<<"cls_R4 interior end"<<std::endl;

			// =====================================================================
            // r = 3 stencil (Centered stencil)
            // =====================================================================

			d_rho_3.reinit(6);			d_rho_u_3.reinit(6);
			d_rho_v_3.reinit(6);		d_rho_w_3.reinit(6);
			d_e_3.reinit(6);

           	for (unsigned int f = 0; f < n_faces; f++) {			
				d_rho_3(f)  = (RHO(neighbor_index[f]) - rho0);
				d_rho_u_3(f)  = (RHO_U(neighbor_index[f]) - rho_u0);
				d_rho_v_3(f)  = (RHO_V(neighbor_index[f]) - rho_v0);
				d_rho_w_3(f)  = (RHO_W(neighbor_index[f]) - rho_w0);
				d_e_3(f)  = (E(neighbor_index[f]) - e0);
			}

			ROWS = cell_diagonal_neighbor_index[c].size(); 
//			////pcout<<"270: "<<ROWS<<std::endl;
            b_rho_3.reinit(ROWS);	b_rho_u_3.reinit(ROWS);	b_rho_v_3.reinit(ROWS);
			b_rho_w_3.reinit(ROWS);	b_e_3.reinit(ROWS);
			index = 0;
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];
				b_rho_3(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
//				if(c == 764) ////pcout<<"solve: b_rho_3\nindex: "<<index<<"\tneighbor: "<<local_neighbor_dof_indices[0]<<"\trho: "<<RHO(local_neighbor_dof_indices[0] )<<std::endl;
				b_rho_u_3(index)   = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0; 
				b_rho_v_3(index)   = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0; 
				b_rho_w_3(index)   = RHO_W(local_neighbor_dof_indices[0] ) - rho_w0; 
				b_e_3(index)   = E(local_neighbor_dof_indices[0] ) - e0; 
				index++; 
			}

			CLS_R3[c].solve(b_rho_3, d_rho_3, rho_coeff_3); 
			CLS_R3[c].solve(b_rho_u_3, d_rho_u_3, rho_u_coeff_3); 
			CLS_R3[c].solve(b_rho_v_3, d_rho_v_3, rho_v_coeff_3); 
			CLS_R3[c].solve(b_rho_w_3, d_rho_w_3, rho_w_coeff_3); 
			CLS_R3[c].solve(b_e_3, d_e_3, e_coeff_3); 

//			if(g_i == global) pcout<<"solve: cell: "<<g_i<<"\nb_rho_3"<<b_rho_3<<"d_rho_3: "<<d_rho_3<<"rho_coeff_3: "<<rho_coeff_3<<std::endl;
			// =====================================================================
            // r = 3 stencil (directional stencil)
            // =====================================================================

	      	for (unsigned int f = 0; f < n_faces; ++f) {

            	if (is_admissible_R3_d[c][f]) {

					d_rho_3_d[f].reinit(5);					d_rho_u_3_d[f].reinit(5);
					d_rho_v_3_d[f].reinit(5);				d_rho_w_3_d[f].reinit(5);
					d_e_3_d[f].reinit(5);

					index = 0; 

					if(f%2 == 0)
						fff = f+1;            	
					else 
						fff = f-1;

	           		for (unsigned int ff = 0; ff < n_faces; ff++) {
						if(fff != ff) {
							d_rho_3_d[f](index)  = (RHO(neighbor_index[ff]) - rho0);
							d_rho_u_3_d[f](index)  = (RHO_U(neighbor_index[ff]) - rho_u0);
							d_rho_v_3_d[f](index)  = (RHO_V(neighbor_index[ff]) - rho_v0);
							d_rho_w_3_d[f](index)  = (RHO_W(neighbor_index[ff]) - rho_w0);
							d_e_3_d[f](index)  = (E(neighbor_index[ff]) - e0);
//							if(g_i == global && f == 0) pcout<<"solve: cell: "<<g_i<<"\tpossible: "<<is_admissible_R3_d[c][f]<<"\tf: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\nindex: "<<index<<"\tneig: "<<neighbor_index[ff]<<"\trho: "<<RHO(neighbor_index[ff])<<"\td_rho_3_d: "<<d_rho_3_d[f](index)<<std::endl;
							index++;
						}
					}

//					if(g_i == global && f == 0) pcout<<"solve: cell: "<<g_i<<"\tpossible: "<<is_admissible_R3_d[c][f]<<"\tf: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\nindex: "<<index<<"\tneig: "<<n_n_face_index[f]<<"\trho: "<<RHO(n_n_face_index[f])<<"\td_rho_3_d: "<<d_rho_3_d[f](index)<<std::endl;

    	            ROWS = cell_neighbor_index[c][f].size() + 1; 
//					////pcout<<"328: "<<ROWS<<std::endl;
					b_rho_3_d[f].reinit(ROWS); b_rho_u_3_d[f].reinit(ROWS); b_rho_v_3_d[f].reinit(ROWS); 
					b_rho_w_3_d[f].reinit(ROWS); 
					b_e_3_d[f].reinit(ROWS);
					index = 0;

					for (unsigned int d = 0; d < ROWS - 1; ++d) {			
						local_neighbor_dof_indices[0] = cell_neighbor_index[c][f][d]; 	//////pcout<<"334: d: "<<d<<"\tb_rho_3_d: "<<b_rho_3_d[f].size()<<std::endl;
						b_rho_3_d[f](index)  = (RHO(local_neighbor_dof_indices[0]) - rho0);
						b_rho_u_3_d[f](index)  = (RHO_U(local_neighbor_dof_indices[0]) - rho_u0);
						b_rho_v_3_d[f](index)  = (RHO_V(local_neighbor_dof_indices[0]) - rho_v0);
						b_rho_w_3_d[f](index)  = (RHO_W(local_neighbor_dof_indices[0]) - rho_w0);
						b_e_3_d[f](index)  = (E(local_neighbor_dof_indices[0]) - e0);
//						if(g_i == global && f == 0) pcout<<"solve: cell: "<<g_i<<"\tpossible: "<<is_admissible_R3_d[c][f]<<"\tf: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\nindex: "<<index<<"\tneig: "<<local_neighbor_dof_indices[0]<<"\trho: "<<RHO(local_neighbor_dof_indices[0])<<"\tb_rho_3_d: "<<b_rho_3_d[f](index)<<std::endl;
						index++;

					} //////pcout<<"341: "<<std::endl;

					b_rho_3_d[f](index)  = (RHO(n_n_face_index[f]) - rho0);
					b_rho_u_3_d[f](index)  = (RHO_U(n_n_face_index[f]) - rho_u0);
					b_rho_v_3_d[f](index)  = (RHO_V(n_n_face_index[f]) - rho_v0);
					b_rho_w_3_d[f](index)  = (RHO_W(n_n_face_index[f]) - rho_w0);
					b_e_3_d[f](index)  = (E(n_n_face_index[f]) - e0);

//					if(g_i == global && f == 2) pcout<<"solve: cell: "<<g_i<<"\tpossible: "<<is_admissible_R3_d[c][f]<<"\tf: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\nb_rho_3_d: "<<b_rho_3_d[f]<<"d_rho_3_d: "<<d_rho_3_d[f]<<"rho_coeff_3_d: "<<rho_coeff_3_d[f]<<std::endl;

					CLS_R3_d[c][f].solve(b_rho_3_d[f], d_rho_3_d[f], rho_coeff_3_d[f]); 
					CLS_R3_d[c][f].solve(b_rho_u_3_d[f], d_rho_u_3_d[f], rho_u_coeff_3_d[f]); 
					CLS_R3_d[c][f].solve(b_rho_v_3_d[f], d_rho_v_3_d[f], rho_v_coeff_3_d[f]); 
					CLS_R3_d[c][f].solve(b_rho_w_3_d[f], d_rho_w_3_d[f], rho_w_coeff_3_d[f]); 
					CLS_R3_d[c][f].solve(b_e_3_d[f], d_e_3_d[f], e_coeff_3_d[f]); //	////pcout<<"347: c: "<<c<<"\tf: "<<f<<std::endl;

				}

	            else {
					rho_coeff_3_d[f] = rho_coeff_3; 
					rho_u_coeff_3_d[f] = rho_u_coeff_3;
					rho_v_coeff_3_d[f] = rho_v_coeff_3;
					rho_w_coeff_3_d[f] = rho_w_coeff_3;
					e_coeff_3_d[f] = e_coeff_3;
				}
//				if(g_i == global && f == 0) pcout<<"solve: cell: "<<g_i<<"\tpossible: "<<is_admissible_R3_d[c][f]<<"\tf: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\nb_rho_3_d: "<<b_rho_3_d[f]<<"d_rho_3_d: "<<d_rho_3_d[f]<<"rho_coeff_3_d: "<<rho_coeff_3_d[f]<<std::endl;
			} // End of face loop for 3rd order directional stencils

/*
			// =====================================================================
			// r = 2 stencil 0 (f0 f2 & f4 )		
			// =====================================================================
			
			b_rho_2_d[0](0) = d_rho_3(0);			b_rho_2_d[0](1) = d_rho_3(2);			b_rho_2_d[0](2) = d_rho_3(4);
			b_rho_u_2_d[0](0) = d_rho_u_3(0);		b_rho_u_2_d[0](1) = d_rho_u_3(2);		b_rho_u_2_d[0](2) = d_rho_u_3(4);
			b_rho_v_2_d[0](0) = d_rho_v_3(0);		b_rho_v_2_d[0](1) = d_rho_v_3(2);		b_rho_v_2_d[0](2) = d_rho_v_3(4);
			b_rho_w_2_d[0](0) = d_rho_w_3(0);		b_rho_w_2_d[0](1) = d_rho_w_3(2);		b_rho_w_2_d[0](2) = d_rho_w_3(4);
			b_e_2_d[0](0) 	  = d_e_3(0);			b_e_2_d[0](1) 	  = d_e_3(2);			b_e_2_d[0](2) 	  = d_e_3(4);

			// =====================================================================
			// r = 2 stencil 1 (f1 f2 & f4 )		
			// =====================================================================
			
			b_rho_2_d[1](0) = d_rho_3(1);			b_rho_2_d[1](1) = d_rho_3(2);			b_rho_2_d[1](2) = d_rho_3(4);
			b_rho_u_2_d[1](0) = d_rho_u_3(1);		b_rho_u_2_d[1](1) = d_rho_u_3(2);		b_rho_u_2_d[1](2) = d_rho_u_3(4);
			b_rho_v_2_d[1](0) = d_rho_v_3(1);		b_rho_v_2_d[1](1) = d_rho_v_3(2);		b_rho_v_2_d[1](2) = d_rho_v_3(4);
			b_rho_w_2_d[1](0) = d_rho_w_3(1);		b_rho_w_2_d[1](1) = d_rho_w_3(2);		b_rho_w_2_d[1](2) = d_rho_w_3(4);
			b_e_2_d[1](0) 	  = d_e_3(1);			b_e_2_d[1](1) 	  = d_e_3(2);			b_e_2_d[1](2) 	  = d_e_3(4);

			// =====================================================================
			// r = 2 stencil 2 (f0 f2 & f5 )		
			// =====================================================================
			
			b_rho_2_d[2](0) = d_rho_3(0);			b_rho_2_d[2](1) = d_rho_3(2);			b_rho_2_d[2](2) = d_rho_3(5);
			b_rho_u_2_d[2](0) = d_rho_u_3(0);		b_rho_u_2_d[2](1) = d_rho_u_3(2);		b_rho_u_2_d[2](2) = d_rho_u_3(5);
			b_rho_v_2_d[2](0) = d_rho_v_3(0);		b_rho_v_2_d[2](1) = d_rho_v_3(2);		b_rho_v_2_d[2](2) = d_rho_v_3(5);
			b_rho_w_2_d[2](0) = d_rho_w_3(0);		b_rho_w_2_d[2](1) = d_rho_w_3(2);		b_rho_w_2_d[2](2) = d_rho_w_3(5);
			b_e_2_d[2](0) 	  = d_e_3(0);			b_e_2_d[2](1) 	  = d_e_3(2);			b_e_2_d[2](2) 	  = d_e_3(5);

			// =====================================================================
			// r = 2 stencil 3 (f1 f2 & f5 )		
			// =====================================================================
			
			b_rho_2_d[3](0) = d_rho_3(1);			b_rho_2_d[3](1) = d_rho_3(2);			b_rho_2_d[3](2) = d_rho_3(5);
			b_rho_u_2_d[3](0) = d_rho_u_3(1);		b_rho_u_2_d[3](1) = d_rho_u_3(2);		b_rho_u_2_d[3](2) = d_rho_u_3(5);
			b_rho_v_2_d[3](0) = d_rho_v_3(1);		b_rho_v_2_d[3](1) = d_rho_v_3(2);		b_rho_v_2_d[3](2) = d_rho_v_3(5);
			b_rho_w_2_d[3](0) = d_rho_w_3(1);		b_rho_w_2_d[3](1) = d_rho_w_3(2);		b_rho_w_2_d[3](2) = d_rho_w_3(5);
			b_e_2_d[3](0) 	  = d_e_3(1);			b_e_2_d[3](1) 	  = d_e_3(2);			b_e_2_d[3](2) 	  = d_e_3(5);

			// =====================================================================
			// r = 2 stencil 4 (f0 f3 & f4 )		
			// =====================================================================
			
			b_rho_2_d[4](0) = d_rho_3(0);			b_rho_2_d[4](1) = d_rho_3(3);			b_rho_2_d[4](2) = d_rho_3(4);
			b_rho_u_2_d[4](0) = d_rho_u_3(0);		b_rho_u_2_d[4](1) = d_rho_u_3(3);		b_rho_u_2_d[4](2) = d_rho_u_3(4);
			b_rho_v_2_d[4](0) = d_rho_v_3(0);		b_rho_v_2_d[4](1) = d_rho_v_3(3);		b_rho_v_2_d[4](2) = d_rho_v_3(4);
			b_rho_w_2_d[4](0) = d_rho_w_3(0);		b_rho_w_2_d[4](1) = d_rho_w_3(3);		b_rho_w_2_d[4](2) = d_rho_w_3(4);
			b_e_2_d[4](0) 	  = d_e_3(0);			b_e_2_d[4](1) 	  = d_e_3(3);			b_e_2_d[4](2) 	  = d_e_3(4);

			// =====================================================================
			// r = 2 stencil 5 (f1 f3 & f4 )		
			// =====================================================================
			
			b_rho_2_d[5](0) = d_rho_3(1);			b_rho_2_d[5](1) = d_rho_3(3);			b_rho_2_d[5](2) = d_rho_3(4);
			b_rho_u_2_d[5](0) = d_rho_u_3(1);		b_rho_u_2_d[5](1) = d_rho_u_3(3);		b_rho_u_2_d[5](2) = d_rho_u_3(4);
			b_rho_v_2_d[5](0) = d_rho_v_3(1);		b_rho_v_2_d[5](1) = d_rho_v_3(3);		b_rho_v_2_d[5](2) = d_rho_v_3(4);
			b_rho_w_2_d[5](0) = d_rho_w_3(1);		b_rho_w_2_d[5](1) = d_rho_w_3(3);		b_rho_w_2_d[5](2) = d_rho_w_3(4);
			b_e_2_d[5](0) 	  = d_e_3(1);			b_e_2_d[5](1) 	  = d_e_3(3);			b_e_2_d[5](2) 	  = d_e_3(4);

			// =====================================================================
			// r = 2 stencil 6 (f0 f3 & f5 )		
			// =====================================================================
			
			b_rho_2_d[6](0) = d_rho_3(0);			b_rho_2_d[6](1) = d_rho_3(3);			b_rho_2_d[6](2) = d_rho_3(5);
			b_rho_u_2_d[6](0) = d_rho_u_3(0);		b_rho_u_2_d[6](1) = d_rho_u_3(3);		b_rho_u_2_d[6](2) = d_rho_u_3(5);
			b_rho_v_2_d[6](0) = d_rho_v_3(0);		b_rho_v_2_d[6](1) = d_rho_v_3(3);		b_rho_v_2_d[6](2) = d_rho_v_3(5);
			b_rho_w_2_d[6](0) = d_rho_w_3(0);		b_rho_w_2_d[6](1) = d_rho_w_3(3);		b_rho_w_2_d[6](2) = d_rho_w_3(5);
			b_e_2_d[6](0) 	  = d_e_3(0);			b_e_2_d[6](1) 	  = d_e_3(3);			b_e_2_d[6](2) 	  = d_e_3(5);

			// =====================================================================
			// r = 2 stencil 7 (f1 f3 & f5 )		
			// =====================================================================
			
			b_rho_2_d[7](0) = d_rho_3(1);			b_rho_2_d[7](1) = d_rho_3(3);			b_rho_2_d[7](2) = d_rho_3(5);
			b_rho_u_2_d[7](0) = d_rho_u_3(1);		b_rho_u_2_d[7](1) = d_rho_u_3(3);		b_rho_u_2_d[7](2) = d_rho_u_3(5);
			b_rho_v_2_d[7](0) = d_rho_v_3(1);		b_rho_v_2_d[7](1) = d_rho_v_3(3);		b_rho_v_2_d[7](2) = d_rho_v_3(5);
			b_rho_w_2_d[7](0) = d_rho_w_3(1);		b_rho_w_2_d[7](1) = d_rho_w_3(3);		b_rho_w_2_d[7](2) = d_rho_w_3(5);
			b_e_2_d[7](0) 	  = d_e_3(1);			b_e_2_d[7](1) 	  = d_e_3(3);			b_e_2_d[7](2) 	  = d_e_3(5);

	   	  	for (unsigned int v = 0; v < 8; v++) {

				LU_R2_d[c][v].solve(b_rho_2_d[v], rho_coeff_2_d[v]);
				LU_R2_d[c][v].solve(b_rho_u_2_d[v], rho_u_coeff_2_d[v]);
				LU_R2_d[c][v].solve(b_rho_v_2_d[v], rho_v_coeff_2_d[v]);
				LU_R2_d[c][v].solve(b_rho_w_2_d[v], rho_w_coeff_2_d[v]);
				LU_R2_d[c][v].solve(b_e_2_d[v], e_coeff_2_d[v]);

			}
*/
/*
			rho_coeff_4 = 0.0;
			rho_u_coeff_4 = 0.0;
			rho_v_coeff_4 = 0.0;
			rho_w_coeff_4 = 0.0;
			e_coeff_4 = 0.0;

           	for (unsigned int i = 0; i < 9; i++) {						

                rho_coeff_4(i) = rho_coeff_3(i);
                rho_u_coeff_4(i) = rho_u_coeff_3(i);
                rho_v_coeff_4(i) = rho_v_coeff_3(i);
                rho_w_coeff_4(i) = rho_w_coeff_3(i);
                e_coeff_4(i) = e_coeff_3(i);

			}
*/


//			if(g_i == 960) pcout<<"interior out g_i: "<<g_i<<"\n";
		} // End of interior cells loop

		else {
//			if(g_i == 960) pcout<<"boundary g_i: "<<g_i<<"\tcenter: "<<cell->center()<<std::endl;
				////pcout<<"corner in\n";
/*
			    fourth_order_wt = 0.76;
				third_order_wt = 0.02; 
				second_order_wt = 0.0; 

				gamma(0) = fourth_order_wt; 
				gamma(1) = third_order_wt; 
				gamma(2) = third_order_wt;
				gamma(3) = third_order_wt;
				gamma(4) = third_order_wt;
				gamma(5) = third_order_wt;
				gamma(6) = third_order_wt;
				gamma(7) = third_order_wt;

				gamma(8) = second_order_wt;
				gamma(9) = second_order_wt;
				gamma(10) = second_order_wt;
				gamma(11) = second_order_wt;
				gamma(12) = second_order_wt;
				gamma(13) = second_order_wt;
				gamma(14) = second_order_wt;
				gamma(15) = second_order_wt;
*/
				unsigned int n_boundary_faces = 0, n_no_slip_boundary = 0;
			    for (unsigned int f=0; f<n_faces; ++f) {
					if(cell->face(f)->at_boundary()) {
						if(cell->face(f)->boundary_id() == 2) {	
							n_no_slip_boundary += 1;
						}
						n_boundary_faces += 1;
					}
				}

				// =====================================================================
	            // r = 4 stencil 
	            // =====================================================================
//			if (n_boundary_faces < 3) {

				index_least = 0;

				ROWS = cell_diagonal_neighbor_index[c].size(); 
			
				// neighbor of neighbors 

		      	for (unsigned int f = 0; f < n_faces; f++)
					if (n_n_face[f])
						ROWS++; 

				ROWS_slip = ROWS + (n_boundary_faces - n_no_slip_boundary)*4;
				ROWS += n_boundary_faces*4;
    	        // Least Squares Part  
//            	if(ROWS == 19) std::cout<<"boundary c: "<<c<<std::endl;
    	        b_rho_4.reinit(ROWS);	
				b_e_4.reinit(ROWS);

				b_rho_u_4.reinit(ROWS_slip);	b_rho_v_4.reinit(ROWS_slip);
				b_rho_w_4.reinit(ROWS_slip);	

			
				// vertex neighbors
			
				for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

					local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];
					b_rho_4(index_least)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_4(index_least)   = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0; 
					b_rho_v_4(index_least)   = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0; 
					b_rho_w_4(index_least)   = RHO_W(local_neighbor_dof_indices[0] ) - rho_w0; 
					b_e_4(index_least)   = E(local_neighbor_dof_indices[0] ) - e0; 
//					if(g_i == 960) pcout<<"neigh g_i: "<<local_neighbor_dof_indices[0]<<"\trho: "<<RHO(local_neighbor_dof_indices[0] )<<"\trho0: "<<rho0<<std::endl;
					index_least++; 
				}
			
				// neighbors of neighbors 

		      	for (unsigned int f = 0; f < n_faces; f++) {
					if (n_n_face[f]){
						local_neighbor_dof_indices[0] = n_n_face_index[f];
						b_rho_4(index_least)   = RHO(local_neighbor_dof_indices[0] ) - rho0;
						b_rho_u_4(index_least)   = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
						b_rho_v_4(index_least)   = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
						b_rho_w_4(index_least)   = RHO_W(local_neighbor_dof_indices[0] ) - rho_w0;
						b_e_4(index_least)   = E(local_neighbor_dof_indices[0] ) - e0;
						index_least++;
					}
    	        }

				index = 0; 
				ROWS_slip = 6 - n_boundary_faces + n_no_slip_boundary*4;
				ROWS = 6 - n_boundary_faces;
				d_rho_4.reinit(ROWS);		
				d_rho_u_4.reinit(ROWS_slip);
				d_rho_v_4.reinit(ROWS_slip);
				d_rho_w_4.reinit(ROWS_slip);
				d_e_4.reinit(ROWS);
				unsigned int index_least_slip = index_least;
	           	for (unsigned int f = 0; f < n_faces; f++) {			

					if(!cell->face(f)->at_boundary()) {	
						neighbor = cell->neighbor(f);
					    neighbor->get_dof_indices(local_neighbor_dof_indices);
						neighbor_index[f] = local_neighbor_dof_indices[0];
						d_rho_4(index)  = (RHO(neighbor_index[f]) - rho0);
						d_e_4(index)  = (E(neighbor_index[f]) - e0);
						index++;
					}
					else {

			                for (unsigned int i = 0; i < 4; i++) {
								b_rho_4(index_least) = 0.0;
								b_e_4(index_least) = 0.0;
								index_least++;
							}
					}
				}

				index = 0; 
				ROWS_slip = 6 - n_boundary_faces + n_no_slip_boundary*4;
				d_rho_u_4.reinit(ROWS_slip);
				d_rho_v_4.reinit(ROWS_slip);
				d_rho_w_4.reinit(ROWS_slip);

	           	for (unsigned int f = 0; f < n_faces; f++) {			

					if(!cell->face(f)->at_boundary()) {	
						neighbor = cell->neighbor(f);
					    neighbor->get_dof_indices(local_neighbor_dof_indices);
						neighbor_index[f] = local_neighbor_dof_indices[0];
						d_rho_u_4(index)  = (RHO_U(neighbor_index[f]) - rho_u0);
						d_rho_v_4(index)  = (RHO_V(neighbor_index[f]) - rho_v0);
						d_rho_w_4(index)  = (RHO_W(neighbor_index[f]) - rho_w0);
						index++;
					}
					else {
						if (cell->face(f)->boundary_id() == 2) {
			                for (unsigned int i = 0; i < 4; i++) {
								d_rho_u_4(index) = -rho_u0;
								d_rho_v_4(index) = -rho_v0;
								d_rho_w_4(index) = -rho_w0;
								index++;
							}
						}
						else {
			                for (unsigned int i = 0; i < 4; i++) {
								b_rho_u_4(index_least_slip) = 0.0;
								b_rho_v_4(index_least_slip) = 0.0;
								b_rho_w_4(index_least_slip) = 0.0;
								index_least_slip++;
							}
						}
					}
				}

//				if(g_i == 960) std::cout<<"upwind nuemann \n"<<"b_rho_4: "<<b_rho_4<<"d_rho_4: "<<d_rho_4<<std::endl;
//				std::cout<<"den c: "<<c<<"\n";
                CLS_R4[c].solve(b_rho_4, d_rho_4, rho_coeff_4);
//				std::cout<<"vel c: "<<c<<"\n";
                CLS_R4_slip[c].solve(b_rho_u_4, d_rho_u_4, rho_u_coeff_4);
                CLS_R4_slip[c].solve(b_rho_v_4, d_rho_v_4, rho_v_coeff_4);
                CLS_R4_slip[c].solve(b_rho_w_4, d_rho_w_4, rho_w_coeff_4);
                CLS_R4[c].solve(b_e_4, d_e_4, e_coeff_4);
//				if(g_i == 960) std::cout<<"upwind\n";
//					if(g_i == 2) pcout<<"solve: cell: "<<g_i<<"\nb_rho_u_4: "<<b_rho_u_4<<"d_rho_u_4: "<<d_rho_u_4<<"rho_u_coeff_4: "<<rho_u_coeff_4<<std::endl;
/*
				if(std::fabs(cell->center()[0] - 1.1) < 0.06 && std::fabs(cell->center()[2] - 0.45) < 0.01) {

					for (unsigned int d = 0; d < cell_all_neighbor_index[c].size(); ++d) {
						local_neighbor_dof_indices[0] = cell_all_neighbor_index[c][d];
						unsigned int local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];			
						neighbor = local_index_to_iterator[local_index];
						std::cout<<"cell center: "<<cell->center()<<"\tRHO : "<<RHO(local_neighbor_dof_indices[0])<<std::endl;
					}	
				}
*/
				// =====================================================================
	            // r = 3 stencil 
	            // =====================================================================
//					//pcout<<"472: 3rd c: "<<c<<std::endl;
					index = 0;
					index_least = 0;
/*
					d_rho_3.reinit(3);			d_rho_u_3.reinit(3);
					d_rho_v_3.reinit(3);		d_rho_w_3.reinit(3);
					d_e_3.reinit(3);
*/
					ROWS = 6 - n_boundary_faces;
					d_rho_3.reinit(ROWS);			d_rho_u_3.reinit(ROWS);
					d_rho_v_3.reinit(ROWS);		d_rho_w_3.reinit(ROWS);
					d_e_3.reinit(ROWS);

					ROWS = cell_diagonal_neighbor_index[c].size() + n_boundary_faces*4;
//					ROWS_slip = ROWS + n_boundary_faces*4;
//					////pcout<<"477: "<<ROWS<<std::endl;
	    	        b_rho_3.reinit(ROWS);	b_rho_u_3.reinit(ROWS);	b_rho_v_3.reinit(ROWS);
					b_rho_w_3.reinit(ROWS);	b_e_3.reinit(ROWS);

//	    	        b_rho_3.reinit(ROWS + 9);	b_rho_u_3.reinit(ROWS + 9);	b_rho_v_3.reinit(ROWS + 9);
//					b_rho_w_3.reinit(ROWS + 9);	b_e_3.reinit(ROWS + 9);

//					if(c == 5376) //pcout<<"c: "<<c<<"\tg_i: "<<g_i<<std::endl;	
		           	for (unsigned int f = 0; f < n_faces; f++) {			
	
						if(!cell->face(f)->at_boundary()) {	
							neighbor = cell->neighbor(f);
						    neighbor->get_dof_indices(local_neighbor_dof_indices);
							neighbor_index[f] = local_neighbor_dof_indices[0];

							d_rho_3(index)  = (RHO(neighbor_index[f]) - rho0);
							d_rho_u_3(index)  = (RHO_U(neighbor_index[f]) - rho_u0);
							d_rho_v_3(index)  = (RHO_V(neighbor_index[f]) - rho_v0);
							d_rho_w_3(index)  = (RHO_W(neighbor_index[f]) - rho_w0);
							d_e_3(index)  = (E(neighbor_index[f]) - e0);
//							if(c == 5376) //pcout<<"index: "<<index<<"\tneighbor_index: "<<neighbor_index[f]<<"\tb_rho_3: "<<b_rho_3(index)<<"\tb_rho_u_3: "<<b_rho_u_3(index)<<std::endl;	
//							if(c == 5376 && neighbor->face(0)->at_boundary()) //pcout<<"index: "<<neighbor_index[f]<<std::endl;
							index++;

						}

						else {

							if (cell->face(f)->boundary_id() == 2) {

					            for (unsigned int i = 0; i < 4; i++) {
									b_rho_3(index_least) = 0.0;
									b_rho_u_3(index_least) = -rho_u0*weight_ls;
									b_rho_v_3(index_least) = -rho_v0*weight_ls;
									b_rho_w_3(index_least) = -rho_w0*weight_ls;
									b_e_3(index_least) = 0.0;
									index_least++;
								}
							}
							else {

					            for (unsigned int i = 0; i < 4; i++) {
									b_rho_3(index_least) = 0.0;
									b_rho_u_3(index_least) = 0.0;
									b_rho_v_3(index_least) = 0.0;
									b_rho_w_3(index_least) = 0.0;
									b_e_3(index_least) = 0.0;
									index_least++;
								}
							}
						}
					} //////pcout<<"500\n";
					for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

						local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];
	
						b_rho_3(index_least)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
						b_rho_u_3(index_least)   = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0; 
						b_rho_v_3(index_least)   = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0; 
						b_rho_w_3(index_least)   = RHO_W(local_neighbor_dof_indices[0] ) - rho_w0; 
						b_e_3(index_least)   = E(local_neighbor_dof_indices[0] ) - e0; 
//						if(c == 5376) //pcout<<"index_least: "<<index_least<<"\tneighbor: "<<local_neighbor_dof_indices[0]<<"\tb_rho_3: "<<b_rho_3(index_least)<<"\tb_rho_u_3: "<<b_rho_u_3(index_least)<<std::endl;	
						index_least++; 
						index++;
					}

					CLS_R3[c].solve(b_rho_3, d_rho_3, rho_coeff_3);
					CLS_R3_slip[c].solve(b_rho_u_3, d_rho_u_3, rho_u_coeff_3); 
					CLS_R3_slip[c].solve(b_rho_v_3, d_rho_v_3, rho_v_coeff_3); 
					CLS_R3_slip[c].solve(b_rho_w_3, d_rho_w_3, rho_w_coeff_3); 
					CLS_R3[c].solve(b_e_3, d_e_3, e_coeff_3);

//					if(g_i == 2) pcout<<"solve: cell: "<<g_i<<"\nb_rho_u_3: "<<b_rho_u_3<<"d_rho_u_3: "<<d_rho_u_3<<"rho_u_coeff_3: "<<rho_u_coeff_3<<std::endl; 
//					if(g_i == 960) std::cout<<"3 upwind\n";
//					if(g_i == global) pcout<<"solve: cell: "<<g_i<<"\nb_rho_3"<<b_rho_3<<"d_rho_3: "<<d_rho_3<<"rho_coeff_3: "<<rho_coeff_3<<std::endl;
//					if(g_i == 1290)
//						pcout<<"rho_coeff_3: "<<rho_coeff_3<<std::endl;

/*
					CLS_R3[c].solve(b_rho_3, rho_coeff_3); 
					CLS_R3[c].solve(b_rho_u_3, rho_u_coeff_3); 
					CLS_R3[c].solve(b_rho_v_3, rho_v_coeff_3); 
					CLS_R3[c].solve(b_rho_w_3, rho_w_coeff_3); 
					CLS_R3[c].solve(b_e_3, e_coeff_3); 					//////pcout<<"525 end: "<<ROWS<<std::endl;
*/
//					////pcout<<"528: 3rd over c: "<<c<<std::endl;
//					if(g_i == 21376) {
/*
					if(g_i == 21376) {
						//pcout<<"solve: b_rho_3: ";
						for(unsigned int i = 0; i < 17; ++i)
							//pcout<<b_rho_3(i)<<"\t";
						//pcout<<std::endl;
						//pcout<<"d_rho_3: "<<d_rho_3<<"rho_coeff_3: "<<rho_coeff_3<<std::endl;
					}
*/
//					if(g_i == 21408) std::cout<<"solve: b_rho_3: "<<b_rho_3<<"d_rho_3: "<<d_rho_3<<"rho_coeff_3: "<<rho_coeff_3<<"rhou: "<<rho_u_coeff_3<<std::endl;

				//pcout<<"715 end: "<<std::endl;
				// =====================================================================
	            // r = 3 directiomal stencil 
	            // =====================================================================

	           	for (unsigned int f = 0; f < n_faces; f++) {	
		
					unsigned int f_e = numbers::invalid_unsigned_int, n_boundary_face_d = 0, n_no_slip_boundary_face_d = 0;;
					bool no_slip_boundary_d = false, transmissive_boundary_d = false;

					if(f%2 == 0)
						f_e = f + 1;
					else
						f_e = f - 1;

					if(is_admissible_R3_d[c][f]) {

		           		for (unsigned int ff = 0; ff < n_faces; ff++)
							if(ff != f_e) 
								if(cell->face(ff)->boundary_id() == 2) {
									no_slip_boundary_d = true;
									n_boundary_face_d += 1;
									n_no_slip_boundary_face_d += 1;
								}
								else if(cell->face(ff)->boundary_id() == 0 || cell->face(ff)->boundary_id() == 1) {
									transmissive_boundary_d = true;
									n_boundary_face_d += 1;
								}


							//pcout<<"732: 3rd direc c: "<<c<<"\tf: "<<f<<std::endl;
								ROWS = 5 - n_boundary_face_d;						
								index = 0;
								index_least = 0;
								d_rho_3_d[f].reinit(ROWS);			d_rho_u_3_d[f].reinit(ROWS);
								d_rho_v_3_d[f].reinit(ROWS);		d_rho_w_3_d[f].reinit(ROWS);
								d_e_3_d[f].reinit(ROWS);

			    	            ROWS = cell_neighbor_index[c][f].size() + 1 + n_boundary_face_d*4; // 4 face neighbor and 1 neighbor of face neighbor in f direction
//								ROWS_slip = ROWS + n_boundary_face_d*4;
								//////pcout<<"615: "<<ROWS<<std::endl;		
					            b_rho_3_d[f].reinit(ROWS);	b_rho_u_3_d[f].reinit(ROWS);	b_rho_v_3_d[f].reinit(ROWS);
								b_rho_w_3_d[f].reinit(ROWS);	b_e_3_d[f].reinit(ROWS);

				           		for (unsigned int ff = 0; ff < n_faces; ff++) {

									if(ff != f_e) {

										if(!cell->face(ff)->at_boundary()) {
											d_rho_3_d[f](index)  = (RHO(neighbor_index[ff]) - rho0);
											d_rho_u_3_d[f](index)  = (RHO_U(neighbor_index[ff]) - rho_u0);
											d_rho_v_3_d[f](index)  = (RHO_V(neighbor_index[ff]) - rho_v0);
											d_rho_w_3_d[f](index)  = (RHO_W(neighbor_index[ff]) - rho_w0);
											d_e_3_d[f](index)  = (E(neighbor_index[ff]) - e0);
//											if(g_i == global && f == 0) pcout<<"solve: cell: "<<g_i<<"\tpossible: "<<is_admissible_R3_d[c][f]<<"\tf: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\nindex: "<<index<<"\tneig: "<<neighbor_index[ff]<<"\trho: "<<RHO(neighbor_index[ff])<<"\td_rho_3_d: "<<d_rho_3_d[f](index)<<std::endl;
											index++;
										}	

										else {

											if (cell->face(ff)->boundary_id() == 2) {
/*
									        	for (unsigned int i = 0; i < 1; i++) {
													d_rho_u_3_d[f](index) = -rho_u0;
													d_rho_v_3_d[f](index) = -rho_v0;
													d_rho_w_3_d[f](index) = -rho_w0;
													index++;
												}
*/
									            for (unsigned int i = 0; i < 4; i++) {
													b_rho_3_d[f](index_least)  = 0.0;
													b_rho_u_3_d[f](index_least) = -rho_u0*weight_ls;
													b_rho_v_3_d[f](index_least) = -rho_v0*weight_ls;
													b_rho_w_3_d[f](index_least) = -rho_w0*weight_ls;
													b_e_3_d[f](index_least)  = 0.0;
													index_least++;
												}
											}
											else {
									            for (unsigned int i = 0; i < 4; i++) {
													b_rho_3_d[f](index_least)  = 0.0;
													b_rho_u_3_d[f](index_least) = 0.0;
													b_rho_v_3_d[f](index_least) = 0.0;
													b_rho_w_3_d[f](index_least) = 0.0;
													b_e_3_d[f](index_least)  = 0.0;
													index_least++;
												}
											}

										}
						  			}
								}


								for (unsigned int d = 0; d < cell_neighbor_index[c][f].size(); ++d) {
	
									local_neighbor_dof_indices[0] = cell_neighbor_index[c][f][d];
									b_rho_3_d[f](index_least)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
									b_rho_u_3_d[f](index_least)   = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0; 
									b_rho_v_3_d[f](index_least)   = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0; 
									b_rho_w_3_d[f](index_least)   = RHO_W(local_neighbor_dof_indices[0] ) - rho_w0; 
									b_e_3_d[f](index_least)   = E(local_neighbor_dof_indices[0] ) - e0; 
									index++; 
									index_least++;
								}

								b_rho_3_d[f](index_least)  = (RHO(n_n_face_index[f]) - rho0);
								b_rho_u_3_d[f](index_least)  = (RHO_U(n_n_face_index[f]) - rho_u0);
								b_rho_v_3_d[f](index_least)  = (RHO_V(n_n_face_index[f]) - rho_v0);
								b_rho_w_3_d[f](index_least)  = (RHO_W(n_n_face_index[f]) - rho_w0);
								b_e_3_d[f](index_least)  = (E(n_n_face_index[f]) - e0);

								CLS_R3_d[c][f].solve(b_rho_3_d[f], d_rho_3_d[f], rho_coeff_3_d[f]); 
								CLS_R3_d_slip[c][f].solve(b_rho_u_3_d[f], d_rho_u_3_d[f], rho_u_coeff_3_d[f]); 
								CLS_R3_d_slip[c][f].solve(b_rho_v_3_d[f], d_rho_v_3_d[f], rho_v_coeff_3_d[f]); 
								CLS_R3_d_slip[c][f].solve(b_rho_w_3_d[f], d_rho_w_3_d[f], rho_w_coeff_3_d[f]); 
								CLS_R3_d[c][f].solve(b_e_3_d[f], d_e_3_d[f], e_coeff_3_d[f]);

/*
								if(c == 36 && f == 0) 
									pcout<<"3d0"<<"\nb_rho_u_3_d: "<<b_rho_u_3_d[f]
											<<"d_rho_u_3_d: "<<d_rho_u_3_d[f]
											<<"rho_u_coeff_3_d: "<<rho_u_coeff_3_d[f];
*/
/*
								if( c == 0) {
									pcout<<"3"<<f<<"\t"<<is_admissible_R3_d[c][f]<<"\nb_rho_u_3_d: "<<b_rho_u_3_d[f]
											<<"d_rho_u_3_d: "<<d_rho_u_3_d[f]
											<<"rho_u_coeff_3_d: "<<rho_u_coeff_3_d[f];
	
								}
*/
							//pcout<<"848: 3rd direc one over c: "<<c<<"\tf: "<<f<<std::endl;
					} // End for admissible stencil

					else {

						rho_coeff_3_d[f] = rho_coeff_3; 
						rho_u_coeff_3_d[f] = rho_u_coeff_3;
						rho_v_coeff_3_d[f] = rho_v_coeff_3;
						rho_w_coeff_3_d[f] = rho_w_coeff_3;
						e_coeff_3_d[f] = e_coeff_3;

					}
//								if(g_i == 960) std::cout<<"3 f: "<<f<<"\n";
//					if(g_i == global) pcout<<"solve: cell: "<<g_i<<"\tpossible: "<<is_admissible_R3_d[c][f]<<"\tf: "<<f<<"\tcenter: "<<cell->face(f)->center()<<"\nb_rho_3_d: "<<b_rho_3_d[f]<<"d_rho_3_d: "<<d_rho_3_d[f]<<"rho_coeff_3_d: "<<rho_coeff_3_d[f]<<std::endl;
//					if(g_i == 1290)
//						pcout<<"rho_coeff_3_d f: "<<f<<"\t"<<rho_coeff_3_d[f]<<std::endl;



				} // End of face loop for directional stencil 
				//pcout<<"999 end: "<<std::endl;


/*
			// =====================================================================
			// r = 2 stencil 0 (f0 f2 & f4 )		
			// =====================================================================
			
			b_rho_2_d[0](0) = d_rho_3(0);			b_rho_2_d[0](1) = d_rho_3(2);			b_rho_2_d[0](2) = d_rho_3(4);
			b_rho_u_2_d[0](0) = d_rho_u_3(0);		b_rho_u_2_d[0](1) = d_rho_u_3(2);		b_rho_u_2_d[0](2) = d_rho_u_3(4);
			b_rho_v_2_d[0](0) = d_rho_v_3(0);		b_rho_v_2_d[0](1) = d_rho_v_3(2);		b_rho_v_2_d[0](2) = d_rho_v_3(4);
			b_rho_w_2_d[0](0) = d_rho_w_3(0);		b_rho_w_2_d[0](1) = d_rho_w_3(2);		b_rho_w_2_d[0](2) = d_rho_w_3(4);
			b_e_2_d[0](0) 	  = d_e_3(0);			b_e_2_d[0](1) 	  = d_e_3(2);			b_e_2_d[0](2) 	  = d_e_3(4);

			// =====================================================================
			// r = 2 stencil 1 (f1 f2 & f4 )		
			// =====================================================================
			
			b_rho_2_d[1](0) = d_rho_3(1);			b_rho_2_d[1](1) = d_rho_3(2);			b_rho_2_d[1](2) = d_rho_3(4);
			b_rho_u_2_d[1](0) = d_rho_u_3(1);		b_rho_u_2_d[1](1) = d_rho_u_3(2);		b_rho_u_2_d[1](2) = d_rho_u_3(4);
			b_rho_v_2_d[1](0) = d_rho_v_3(1);		b_rho_v_2_d[1](1) = d_rho_v_3(2);		b_rho_v_2_d[1](2) = d_rho_v_3(4);
			b_rho_w_2_d[1](0) = d_rho_w_3(1);		b_rho_w_2_d[1](1) = d_rho_w_3(2);		b_rho_w_2_d[1](2) = d_rho_w_3(4);
			b_e_2_d[1](0) 	  = d_e_3(1);			b_e_2_d[1](1) 	  = d_e_3(2);			b_e_2_d[1](2) 	  = d_e_3(4);

			// =====================================================================
			// r = 2 stencil 2 (f0 f2 & f5 )		
			// =====================================================================
			
			b_rho_2_d[2](0) = d_rho_3(0);			b_rho_2_d[2](1) = d_rho_3(2);			b_rho_2_d[2](2) = d_rho_3(5);
			b_rho_u_2_d[2](0) = d_rho_u_3(0);		b_rho_u_2_d[2](1) = d_rho_u_3(2);		b_rho_u_2_d[2](2) = d_rho_u_3(5);
			b_rho_v_2_d[2](0) = d_rho_v_3(0);		b_rho_v_2_d[2](1) = d_rho_v_3(2);		b_rho_v_2_d[2](2) = d_rho_v_3(5);
			b_rho_w_2_d[2](0) = d_rho_w_3(0);		b_rho_w_2_d[2](1) = d_rho_w_3(2);		b_rho_w_2_d[2](2) = d_rho_w_3(5);
			b_e_2_d[2](0) 	  = d_e_3(0);			b_e_2_d[2](1) 	  = d_e_3(2);			b_e_2_d[2](2) 	  = d_e_3(5);

			// =====================================================================
			// r = 2 stencil 3 (f1 f2 & f5 )		
			// =====================================================================
			
			b_rho_2_d[3](0) = d_rho_3(1);			b_rho_2_d[3](1) = d_rho_3(2);			b_rho_2_d[3](2) = d_rho_3(5);
			b_rho_u_2_d[3](0) = d_rho_u_3(1);		b_rho_u_2_d[3](1) = d_rho_u_3(2);		b_rho_u_2_d[3](2) = d_rho_u_3(5);
			b_rho_v_2_d[3](0) = d_rho_v_3(1);		b_rho_v_2_d[3](1) = d_rho_v_3(2);		b_rho_v_2_d[3](2) = d_rho_v_3(5);
			b_rho_w_2_d[3](0) = d_rho_w_3(1);		b_rho_w_2_d[3](1) = d_rho_w_3(2);		b_rho_w_2_d[3](2) = d_rho_w_3(5);
			b_e_2_d[3](0) 	  = d_e_3(1);			b_e_2_d[3](1) 	  = d_e_3(2);			b_e_2_d[3](2) 	  = d_e_3(5);

			// =====================================================================
			// r = 2 stencil 4 (f0 f3 & f4 )		
			// =====================================================================
			
			b_rho_2_d[4](0) = d_rho_3(0);			b_rho_2_d[4](1) = d_rho_3(3);			b_rho_2_d[4](2) = d_rho_3(4);
			b_rho_u_2_d[4](0) = d_rho_u_3(0);		b_rho_u_2_d[4](1) = d_rho_u_3(3);		b_rho_u_2_d[4](2) = d_rho_u_3(4);
			b_rho_v_2_d[4](0) = d_rho_v_3(0);		b_rho_v_2_d[4](1) = d_rho_v_3(3);		b_rho_v_2_d[4](2) = d_rho_v_3(4);
			b_rho_w_2_d[4](0) = d_rho_w_3(0);		b_rho_w_2_d[4](1) = d_rho_w_3(3);		b_rho_w_2_d[4](2) = d_rho_w_3(4);
			b_e_2_d[4](0) 	  = d_e_3(0);			b_e_2_d[4](1) 	  = d_e_3(3);			b_e_2_d[4](2) 	  = d_e_3(4);

			// =====================================================================
			// r = 2 stencil 5 (f1 f3 & f4 )		
			// =====================================================================
			
			b_rho_2_d[5](0) = d_rho_3(1);			b_rho_2_d[5](1) = d_rho_3(3);			b_rho_2_d[5](2) = d_rho_3(4);
			b_rho_u_2_d[5](0) = d_rho_u_3(1);		b_rho_u_2_d[5](1) = d_rho_u_3(3);		b_rho_u_2_d[5](2) = d_rho_u_3(4);
			b_rho_v_2_d[5](0) = d_rho_v_3(1);		b_rho_v_2_d[5](1) = d_rho_v_3(3);		b_rho_v_2_d[5](2) = d_rho_v_3(4);
			b_rho_w_2_d[5](0) = d_rho_w_3(1);		b_rho_w_2_d[5](1) = d_rho_w_3(3);		b_rho_w_2_d[5](2) = d_rho_w_3(4);
			b_e_2_d[5](0) 	  = d_e_3(1);			b_e_2_d[5](1) 	  = d_e_3(3);			b_e_2_d[5](2) 	  = d_e_3(4);

			// =====================================================================
			// r = 2 stencil 6 (f0 f3 & f5 )		
			// =====================================================================
			
			b_rho_2_d[6](0) = d_rho_3(0);			b_rho_2_d[6](1) = d_rho_3(3);			b_rho_2_d[6](2) = d_rho_3(5);
			b_rho_u_2_d[6](0) = d_rho_u_3(0);		b_rho_u_2_d[6](1) = d_rho_u_3(3);		b_rho_u_2_d[6](2) = d_rho_u_3(5);
			b_rho_v_2_d[6](0) = d_rho_v_3(0);		b_rho_v_2_d[6](1) = d_rho_v_3(3);		b_rho_v_2_d[6](2) = d_rho_v_3(5);
			b_rho_w_2_d[6](0) = d_rho_w_3(0);		b_rho_w_2_d[6](1) = d_rho_w_3(3);		b_rho_w_2_d[6](2) = d_rho_w_3(5);
			b_e_2_d[6](0) 	  = d_e_3(0);			b_e_2_d[6](1) 	  = d_e_3(3);			b_e_2_d[6](2) 	  = d_e_3(5);

			// =====================================================================
			// r = 2 stencil 7 (f1 f3 & f5 )		
			// =====================================================================
			
			b_rho_2_d[7](0) = d_rho_3(1);			b_rho_2_d[7](1) = d_rho_3(3);			b_rho_2_d[7](2) = d_rho_3(5);
			b_rho_u_2_d[7](0) = d_rho_u_3(1);		b_rho_u_2_d[7](1) = d_rho_u_3(3);		b_rho_u_2_d[7](2) = d_rho_u_3(5);
			b_rho_v_2_d[7](0) = d_rho_v_3(1);		b_rho_v_2_d[7](1) = d_rho_v_3(3);		b_rho_v_2_d[7](2) = d_rho_v_3(5);
			b_rho_w_2_d[7](0) = d_rho_w_3(1);		b_rho_w_2_d[7](1) = d_rho_w_3(3);		b_rho_w_2_d[7](2) = d_rho_w_3(5);
			b_e_2_d[7](0) 	  = d_e_3(1);			b_e_2_d[7](1) 	  = d_e_3(3);			b_e_2_d[7](2) 	  = d_e_3(5);

	   	  	for (unsigned int v = 0; v < 8; v++) {

				LU_R2_d[c][v].solve(b_rho_2_d[v], rho_coeff_2_d[v]);
				LU_R2_d_slip[c][v].solve(b_rho_u_2_d[v], rho_u_coeff_2_d[v]);
				LU_R2_d_slip[c][v].solve(b_rho_v_2_d[v], rho_v_coeff_2_d[v]);
				LU_R2_d_slip[c][v].solve(b_rho_w_2_d[v], rho_w_coeff_2_d[v]);
				LU_R2_d[c][v].solve(b_e_2_d[v], e_coeff_2_d[v]);
//				if(g_i == 1290)
//					pcout<<"rho_coeff_2_d v: "<<v<<"\t"<<rho_coeff_2_d[v]<<std::endl;

			}
*/
/*
			if (n_boundary_faces > 2) {
				rho_coeff_4 = 0.0;
				rho_u_coeff_4 = 0.0;
				rho_v_coeff_4 = 0.0;
				rho_w_coeff_4 = 0.0;
				e_coeff_4 = 0.0;

    	       	for (unsigned int i = 0; i < 9; i++) {						

    	            rho_coeff_4(i) = rho_coeff_3(i);
    	            rho_u_coeff_4(i) = rho_u_coeff_3(i);
    	            rho_v_coeff_4(i) = rho_v_coeff_3(i);
    	            rho_w_coeff_4(i) = rho_w_coeff_3(i);
    	            e_coeff_4(i) = e_coeff_3(i);

				}
			}
*/
/*
				rho_u_coeff_4 = 0.0;
    	       	for (unsigned int i = 0; i < 9; i++) {						

    	            rho_coeff_4(i) = rho_coeff_3(i);
    	            rho_u_coeff_4(i) = rho_u_coeff_3(i);
    	            rho_v_coeff_4(i) = rho_v_coeff_3(i);
    	            rho_w_coeff_4(i) = rho_w_coeff_3(i);
    	            e_coeff_4(i) = e_coeff_3(i);

				}
*/
//			if(g_i == 960) pcout<<" boundary g_i: "<<g_i<<"\n";
		} // End of boundary cells loop

		// Compute smoothness indicator 
//		////pcout<<"843:\n";
/*
				rho_coeff_4 = 0.0;
				rho_u_coeff_4 = 0.0;
				rho_v_coeff_4 = 0.0;
				rho_w_coeff_4 = 0.0;
				e_coeff_4 = 0.0;

    	       	for (unsigned int i = 0; i < 9; i++) {						

    	            rho_coeff_4(i) = rho_coeff_3(i);
    	            rho_u_coeff_4(i) = rho_u_coeff_3(i);
    	            rho_v_coeff_4(i) = rho_v_coeff_3(i);
    	            rho_w_coeff_4(i) = rho_w_coeff_3(i);
    	            e_coeff_4(i) = e_coeff_3(i);

				}
*/
		// =====================================================================
		// r = 4 Stencil 
		// =====================================================================

		IS_RHO(0)   = compute_smoothness_indicator(rho_coeff_4);
		IS_RHO_U(0)   = compute_smoothness_indicator(rho_u_coeff_4);
		IS_RHO_V(0)   = compute_smoothness_indicator(rho_v_coeff_4);
		IS_RHO_W(0)   = compute_smoothness_indicator(rho_w_coeff_4);
		IS_E(0)   = compute_smoothness_indicator(e_coeff_4);

		// =====================================================================
		// r = 3 centred Stencil 
		// =====================================================================

		IS_RHO(1)   = compute_smoothness_indicator(rho_coeff_3);
		IS_RHO_U(1)   = compute_smoothness_indicator(rho_u_coeff_3);
		IS_RHO_V(1)   = compute_smoothness_indicator(rho_v_coeff_3);
		IS_RHO_W(1)   = compute_smoothness_indicator(rho_w_coeff_3);
		IS_E(1)   = compute_smoothness_indicator(e_coeff_3);

		// =====================================================================
		// r = 3 directional Stencil 
		// =====================================================================

		for(unsigned int f = 2; f < 8; ++f) {

			IS_RHO(f)   = compute_smoothness_indicator(rho_coeff_3_d[f-2]); 
			IS_RHO_U(f)   = compute_smoothness_indicator(rho_u_coeff_3_d[f-2]); 
			IS_RHO_V(f)   = compute_smoothness_indicator(rho_v_coeff_3_d[f-2]); 
			IS_RHO_W(f)   = compute_smoothness_indicator(rho_w_coeff_3_d[f-2]); 
			IS_E(f)   = compute_smoothness_indicator(e_coeff_3_d[f-2]); 

		}
/*
		// =====================================================================
		// r = 2 directional Stencil 
		// =====================================================================

		for(unsigned int v = 8; v < 16; ++v) {

			IS_RHO(v)   = compute_smoothness_indicator(rho_coeff_2_d[v-8]); 
			IS_RHO_U(v) = compute_smoothness_indicator(rho_u_coeff_2_d[v-8]); 
			IS_RHO_V(v) = compute_smoothness_indicator(rho_v_coeff_2_d[v-8]); 
			IS_RHO_W(v) = compute_smoothness_indicator(rho_w_coeff_2_d[v-8]);
			IS_E(v)     = compute_smoothness_indicator(e_coeff_2_d[v-8]); 

		}
*/
		// Combine the polynomials 

		sum_RHO = 0.0; sum_RHO_U = 0.0; sum_RHO_V = 0.0; sum_RHO_W = 0.0; sum_E = 0.0; 
		sum_gamma = 0.0; 		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			
			w_RHO(j)   = gamma(j)/(std::pow((IS_RHO(j) + epsilon), p));
			w_RHO_U(j) = gamma(j)/(std::pow((IS_RHO_U(j) + epsilon), p));
			w_RHO_V(j) = gamma(j)/(std::pow((IS_RHO_V(j) + epsilon), p));
			w_RHO_W(j) = gamma(j)/(std::pow((IS_RHO_W(j) + epsilon), p));
			w_E(j)    	 = gamma(j)/(std::pow((IS_E(j) + epsilon), p));
			
			sum_RHO += w_RHO(j); sum_RHO_U += w_RHO_U(j); sum_RHO_V += w_RHO_V(j); sum_RHO_W += w_RHO_W(j); sum_E += w_E(j); 
			sum_gamma += gamma(j); 
		}
		
		// Normalize the weight_lss 
		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			w_RHO(j) = w_RHO(j)/sum_RHO; 
			w_RHO_U(j) = w_RHO_U(j)/sum_RHO_U;
			w_RHO_V(j) = w_RHO_V(j)/sum_RHO_V;
			w_RHO_W(j) = w_RHO_W(j)/sum_RHO_W;
			w_E(j) = w_E(j)/sum_E;
			gamma(j) = gamma(j)/sum_gamma; 
		}

		coeffs_RHO[c](0) = rho0;
		coeffs_RHO_U[c](0) = rho_u0;
		coeffs_RHO_V[c](0) = rho_v0;
		coeffs_RHO_W[c](0) = rho_w0;
		coeffs_E[c](0) = e0;

		upwind_coeffs_RHO[c](0) = rho0;
		upwind_coeffs_RHO_U[c](0) = rho_u0;
		upwind_coeffs_RHO_V[c](0) = rho_v0;
		upwind_coeffs_RHO_W[c](0) = rho_w0;
		upwind_coeffs_E[c](0) = e0;

		for (unsigned int i = 0; i < 19; ++i) {
			upwind_coeffs_RHO[c](i+1) = rho_coeff_4(i);
			upwind_coeffs_RHO_U[c](i+1) = rho_u_coeff_4(i);
			upwind_coeffs_RHO_V[c](i+1) = rho_v_coeff_4(i);
			upwind_coeffs_RHO_W[c](i+1) = rho_w_coeff_4(i);
			upwind_coeffs_E[c](i+1) = e_coeff_4(i); 
		}

		// 1st order derrivative
/*
		if(g_i == 127) pcout<<"all coeff:\n4th: "<<w_RHO_U(0)<<"\tp: "<<rho_u_coeff_4
							<<"\n3rd: "<<w_RHO_U(1)<<"\tp: "<<rho_u_coeff_3
							<<"\n3rd d1: "<<is_admissible_R3_d[c][0]<<"\t"<<w_RHO_U(2)<<"\tp: "<<rho_u_coeff_3_d[0]
							<<"\n3rd d2: "<<is_admissible_R3_d[c][1]<<"\t"<<w_RHO_U(3)<<"\tp: "<<rho_u_coeff_3_d[1]
							<<"\n3rd d3: "<<is_admissible_R3_d[c][2]<<"\t"<<w_RHO_U(4)<<"\tp: "<<rho_u_coeff_3_d[2]
							<<"\n3rd d4: "<<is_admissible_R3_d[c][3]<<"\t"<<w_RHO_U(5)<<"\tp: "<<rho_u_coeff_3_d[3]
							<<"\n3rd d5: "<<is_admissible_R3_d[c][4]<<"\t"<<w_RHO_U(6)<<"\tp: "<<rho_u_coeff_3_d[4]
							<<"\n3rd d6: "<<is_admissible_R3_d[c][5]<<"\t"<<w_RHO_U(7)<<"\tp: "<<rho_u_coeff_3_d[5]
							<<std::endl;
*/

		for(unsigned int co = 1; co < 4; ++co){

			unsigned int co_1 = co - 1;
			coeffs_RHO[c][co]   = (w_RHO(0)/gamma(0))*(rho_coeff_4(co_1) 	   - gamma(1)*rho_coeff_3(co_1))   + w_RHO(1)*rho_coeff_3(co_1);
			coeffs_RHO_U[c][co] = (w_RHO_U(0)/gamma(0)) * (rho_u_coeff_4(co_1) - gamma(1)*rho_u_coeff_3(co_1)) + w_RHO_U(1)*rho_u_coeff_3(co_1);
			coeffs_RHO_V[c][co] = (w_RHO_V(0)/gamma(0)) * (rho_v_coeff_4(co_1) - gamma(1)*rho_v_coeff_3(co_1)) + w_RHO_V(1)*rho_v_coeff_3(co_1);
			coeffs_RHO_W[c][co] = (w_RHO_W(0)/gamma(0)) * (rho_w_coeff_4(co_1) - gamma(1)*rho_w_coeff_3(co_1)) + w_RHO_W(1)*rho_w_coeff_3(co_1);
			coeffs_E[c][co]     = (w_E(0)/gamma(0)) * (e_coeff_4(co_1)		   - gamma(1)*e_coeff_3(co_1))     + w_E(1)*e_coeff_3(co_1);


			for(unsigned int f = 0; f < 6; ++f) {
				coeffs_RHO[c][co]   = coeffs_RHO[c][co]   - (w_RHO(0)/gamma(0))*gamma(f+2)*rho_coeff_3_d[f](co_1)     + w_RHO(f+2)*rho_coeff_3_d[f](co_1);
				coeffs_RHO_U[c][co] = coeffs_RHO_U[c][co] - (w_RHO_U(0)/gamma(0))*gamma(f+2)*rho_u_coeff_3_d[f](co_1) + w_RHO_U(f+2)*rho_u_coeff_3_d[f](co_1); 
				coeffs_RHO_V[c][co] = coeffs_RHO_V[c][co] - (w_RHO_V(0)/gamma(0))*gamma(f+2)*rho_v_coeff_3_d[f](co_1) + w_RHO_V(f+2)*rho_v_coeff_3_d[f](co_1); 
				coeffs_RHO_W[c][co] = coeffs_RHO_W[c][co] - (w_RHO_W(0)/gamma(0))*gamma(f+2)*rho_w_coeff_3_d[f](co_1) + w_RHO_W(f+2)*rho_w_coeff_3_d[f](co_1); 
				coeffs_E[c][co] 	= coeffs_E[c][co] 	  - (w_E(0)/gamma(0))*gamma(f+2)*e_coeff_3_d[f](co_1)         + w_E(f+2)*e_coeff_3_d[f](co_1);						
			}
/*
			for(unsigned int v = 0; v < 8; ++v) {
				coeffs_RHO[c][co] 	= coeffs_RHO[c][co]   - (w_RHO(0)/gamma(0))*gamma(v+8)*rho_coeff_2_d[v](co_1) 	  + w_RHO(v+8)*rho_coeff_2_d[v](co_1);
				coeffs_RHO_U[c][co] = coeffs_RHO_U[c][co] - (w_RHO_U(0)/gamma(0))*gamma(v+8)*rho_u_coeff_2_d[v](co_1) + w_RHO_U(v+8)*rho_u_coeff_2_d[v](co_1); 
				coeffs_RHO_V[c][co] = coeffs_RHO_V[c][co] - (w_RHO_V(0)/gamma(0))*gamma(v+8)*rho_v_coeff_2_d[v](co_1) + w_RHO_V(v+8)*rho_v_coeff_2_d[v](co_1); 
				coeffs_RHO_W[c][co] = coeffs_RHO_W[c][co] - (w_RHO_W(0)/gamma(0))*gamma(v+8)*rho_w_coeff_2_d[v](co_1) + w_RHO_W(v+8)*rho_w_coeff_2_d[v](co_1); 
				coeffs_E[c][co] 	= coeffs_E[c][co]	  - (w_E(0)/gamma(0))*gamma(v+8)*e_coeff_2_d[v](co_1)  		  + w_E(v+8)*e_coeff_2_d[v](co_1);

			}
*/
		}
//		////pcout<<"932:\n";
		// 2nd order derrivative

		for(unsigned int co = 4; co < 10; ++co){

			unsigned int co_1 = co - 1;

			coeffs_RHO[c][co]   = (w_RHO(0)/gamma(0))*(rho_coeff_4(co_1)	   - gamma(1)*rho_coeff_3(co_1))   + w_RHO(1)*rho_coeff_3(co_1);
			coeffs_RHO_U[c][co] = (w_RHO_U(0)/gamma(0)) * (rho_u_coeff_4(co_1) - gamma(1)*rho_u_coeff_3(co_1)) + w_RHO_U(1)*rho_u_coeff_3(co_1);
			coeffs_RHO_V[c][co] = (w_RHO_V(0)/gamma(0)) * (rho_v_coeff_4(co_1) - gamma(1)*rho_v_coeff_3(co_1)) + w_RHO_V(1)*rho_v_coeff_3(co_1);
			coeffs_RHO_W[c][co] = (w_RHO_W(0)/gamma(0)) * (rho_w_coeff_4(co_1) - gamma(1)*rho_w_coeff_3(co_1)) + w_RHO_W(1)*rho_w_coeff_3(co_1);
			coeffs_E[c][co]    	= (w_E(0)/gamma(0)) * (e_coeff_4(co_1) 		   - gamma(1)*e_coeff_3(co_1))     + w_E(1)*e_coeff_3(co_1);

			for(unsigned int f = 0; f < 6; ++f) {
				coeffs_RHO[c][co] = coeffs_RHO[c][co]	  - (w_RHO(0)/gamma(0))*gamma(f+2)*rho_coeff_3_d[f](co_1)	  + w_RHO(f+2)*rho_coeff_3_d[f](co_1);
				coeffs_RHO_U[c][co] = coeffs_RHO_U[c][co] - (w_RHO_U(0)/gamma(0))*gamma(f+2)*rho_u_coeff_3_d[f](co_1) + w_RHO_U(f+2)*rho_u_coeff_3_d[f](co_1); 
				coeffs_RHO_V[c][co] = coeffs_RHO_V[c][co] - (w_RHO_V(0)/gamma(0))*gamma(f+2)*rho_v_coeff_3_d[f](co_1) + w_RHO_V(f+2)*rho_v_coeff_3_d[f](co_1); 
				coeffs_RHO_W[c][co] = coeffs_RHO_W[c][co] - (w_RHO_W(0)/gamma(0))*gamma(f+2)*rho_w_coeff_3_d[f](co_1) + w_RHO_W(f+2)*rho_w_coeff_3_d[f](co_1); 
				coeffs_E[c][co] 	= coeffs_E[c][co]	  - (w_E(0)/gamma(0))*gamma(f+2)*e_coeff_3_d[f](co_1)  		  + w_E(f+2)*e_coeff_3_d[f](co_1);						

			}
		}

		// 3rd order derrivative

		for(unsigned int co = 10; co < 20; ++co){
			unsigned int co_1 = co - 1;
			coeffs_RHO[c][co] 	 = (w_RHO(0)/gamma(0)) * rho_coeff_4(co_1);
			coeffs_RHO_U[c][co]  = (w_RHO_U(0)/gamma(0)) * rho_u_coeff_4(co_1);
			coeffs_RHO_V[c][co]  = (w_RHO_V(0)/gamma(0)) * rho_v_coeff_4(co_1);
			coeffs_RHO_W[c][co]  = (w_RHO_W(0)/gamma(0)) * rho_w_coeff_4(co_1);
			coeffs_E[c][co]      = (w_E(0)/gamma(0)) * e_coeff_4(co_1);
      	}
/*
		if( c == 0) {
			pcout<<"4r: "<<w_RHO_U(0)<<"\tcoeff: "<<rho_u_coeff_4
				<<"3r: "<<w_RHO_U(1)<<"\tcoeff: "<<rho_u_coeff_3
				<<"31: "<<is_admissible_R3_d[c][0]<<"\t"<<w_RHO_U(2)<<"\tcoeff: "<<rho_u_coeff_3_d[0]
				<<"32: "<<is_admissible_R3_d[c][1]<<"\t"<<w_RHO_U(3)<<"\tcoeff: "<<rho_u_coeff_3_d[1]
				<<"33: "<<is_admissible_R3_d[c][2]<<"\t"<<w_RHO_U(4)<<"\tcoeff: "<<rho_u_coeff_3_d[2]
				<<"34: "<<is_admissible_R3_d[c][3]<<"\t"<<w_RHO_U(5)<<"\tcoeff: "<<rho_u_coeff_3_d[3]
				<<"35: "<<is_admissible_R3_d[c][4]<<"\t"<<w_RHO_U(6)<<"\tcoeff: "<<rho_u_coeff_3_d[4]
				<<"36: "<<is_admissible_R3_d[c][5]<<"\t"<<w_RHO_U(7)<<"\tcoeff: "<<rho_u_coeff_3_d[5]
				<<"final: "<<coeffs_RHO_U[c];

		}
*/
/*
		if(g_i == 1289) {
			std::cout<<"g_i: "<<g_i<<"\tcell center: "<<cell->center()
					<<"\n3rd c: "<<rho_coeff_3;
			for(unsigned int f = 0; f < 6; ++f){
				std::cout<<"f: "<<f<<"\tcenter: "<<cell->face(f)->center()
					<<"\n3rd d: "<<rho_coeff_3_d[f];
			}		
			std::cout<<"RHO coeff: "<<coeffs_RHO[c]<<std::endl;
		}
*/

//		if(g_i == global)
//			std::cout<<"g_i: "<<g_i<<"\nRHO coeff: "<<coeffs_RHO[c]<<std::endl;

/*
		for(unsigned int co = 1; co < 20; ++co){

			unsigned int co_1 = co - 1;
			coeffs_RHO[c][co] 	 = rho_coeff_4(co_1);
			coeffs_RHO_U[c][co]  = rho_u_coeff_4(co_1);
			coeffs_RHO_V[c][co]  = rho_v_coeff_4(co_1);
			coeffs_RHO_W[c][co]  = rho_w_coeff_4(co_1);
			coeffs_E[c][co]  = e_coeff_4(co_1);
      	}


		for(unsigned int co = 1; co < 10; ++co){

			unsigned int co_1 = co - 1;
			coeffs_RHO[c][co] 	 = rho_coeff_3(co_1);
			coeffs_RHO_U[c][co]  = rho_u_coeff_3(co_1);
			coeffs_RHO_V[c][co]  = rho_v_coeff_3(co_1);
			coeffs_RHO_W[c][co]  = rho_w_coeff_3(co_1);
			coeffs_E[c][co]  = e_coeff_3(co_1);
      	}


		for(unsigned int co = 1; co < 10; ++co){

			unsigned int co_1 = co - 1, f = 1;
			coeffs_RHO[c][co] 	 = rho_coeff_3_d[f](co_1);
			coeffs_RHO_U[c][co]  = rho_u_coeff_3_d[f](co_1);
			coeffs_RHO_V[c][co]  = rho_v_coeff_3_d[f](co_1);
			coeffs_RHO_W[c][co]  = rho_w_coeff_3_d[f](co_1);
			coeffs_E[c][co]  = e_coeff_3_d[f](co_1);
      	}
*/
/*
        for (unsigned int i = 10; i < 20; i++) {
	        coeffs_RHO[c](i) = 0.0; 
            coeffs_RHO_U[c](i) = 0.0;
            coeffs_RHO_V[c](i) = 0.0;
            coeffs_RHO_W[c](i) = 0.0;
            coeffs_E[c](i) = 0.0;
        }
*/
/*
		h = Cell[c].h();
		negative = false;
        for (unsigned int f = 0; f < n_faces; ++f) {

	    	for (unsigned int q = 0; q < 4; q++) {

                face_quadrature_point = Cell[c].face_quadrature_point(f,q); 

	            U(0) = evaluate_weno_polynomial(coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point, h);
	            U(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point, h);
	            U(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point, h);
	            U(3) = evaluate_weno_polynomial(coeffs_RHO_W[c], WENO_poly_consts[c], face_quadrature_point, h);
	            U(4) = evaluate_weno_polynomial(coeffs_E[c], WENO_poly_consts[c], face_quadrature_point, h);

				claw.conserved_to_primitive(U,W);


	            if (W(0) < 0.0 || W(4) < 0.0) {
	
    	            for (unsigned int i = 1; i < 20; i++) {
    	              coeffs_RHO[c](i) = 0.0; 
    	              coeffs_RHO_U[c](i) = 0.0;
    	              coeffs_RHO_V[c](i) = 0.0;
    	              coeffs_RHO_W[c](i) = 0.0;
    	              coeffs_E[c](i) = 0.0;
    	            }
					is_1st_order[c] = true;
					negative = true;
    	        }

				if(negative) break;
    	    }

			if(negative) break;
		}
*/
//		////pcout<<"1000:\n";
//	if(c == 764) ////pcout<<"recon: "<<coeffs_RHO[c]<<std::endl;
    } // End of cell loop 
/*
    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
		std::ofstream fout_convergence ; 
		fout_convergence.flags( std::ios::dec | std::ios::scientific ) ; 
		fout_convergence.precision(8) ;

    	const std::string filename = "recon_quad.dat";
		fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);
	
	   	fout_convergence << weight_ls
   					<< std::endl;
		fout_convergence.close();
	}
*/
//	pcout<<"global_gradRho_max: "<<global_gradRho_max<<"\nglobal_rho_slip_max: "<<global_rho_slip_max<<"\nglobal_u_slip_max: "<<global_u_slip_max<<std::endl;
} // End of function 

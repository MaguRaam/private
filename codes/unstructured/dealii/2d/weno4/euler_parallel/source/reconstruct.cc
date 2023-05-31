#include "../include/Weno432.h"


// Perform the actual reconstruction 

void Weno4_2D::reconstruct() {
	
    Vector<double> W1(4), U1(4); 
    Vector<double> W2(4), U2(4);

    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;
    Point<2> face_quadrature_point_1;
    Point<2> face_quadrature_point_2;
/*
	// Gauss Quadrature 
    unsigned int N_gp = 2;  // No. of quadrature points
    QGauss<2> quadrature_formula(N_gp);
    FEValues<2> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);
    Point<2> q_point;
  
    QGauss<2-1> face_quadrature_formula(2);
    FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);
*/
//    Tensor<1,2> face_normal_vector; // Face normal vector

//    double nx, ny;   // Face normal vectors
     
    unsigned int no_stencils = 10;
    unsigned int p = 4; 
    double epsilon = 1.0e-12; 
    double h; // Measure of cell size 
    
    double fourth_order_wt = 0.8;
	double third_order_wt = 0.03; 
	double second_order_wt = 0.0125; 
	double sum_gamma; 
	
	Vector<double> gamma(no_stencils);
	
	gamma(0) = fourth_order_wt; 
	gamma(1) = third_order_wt; 
	gamma(2) = third_order_wt;
	gamma(3) = third_order_wt;
	gamma(4) = third_order_wt;
	gamma(5) = third_order_wt;
	gamma(6) = second_order_wt;
	gamma(7) = second_order_wt;
	gamma(8) = second_order_wt;
	gamma(9) = second_order_wt;
    
    // Variables for reconstruction of RHO
    
	double rho0; 
    
	/* Fourth order stencil */ 
	Vector<double> d_rho_4(4); 
	Vector<double> b_rho_4; 
	Vector<double> rho_coeff_4(9); 
    
	/* Third order centered stencil */ 
	Vector<double> d_rho_3(4); 
	Vector<double> b_rho_3; 
	Vector<double> rho_coeff_3(5); 
    
	/* First third order one-sided stencil */ 
	Vector<double> b_rho_31;         
	Vector<double> d_rho_31;
	Vector<double> rho_coeff_31(5);
    
	/* Second third order one-sided stencil */
	Vector<double> b_rho_32;       
	Vector<double> d_rho_32;
	Vector<double> rho_coeff_32(5);
    
    /* Third third order one-sided stencil */
	Vector<double> b_rho_33;          
	Vector<double> d_rho_33;
	Vector<double> rho_coeff_33(5);
    
    /* Fourth third order one-sided stencil */
    Vector<double> b_rho_34; 
    Vector<double> d_rho_34;
    Vector<double> rho_coeff_34(5);  
	
	Vector<double> rho_coeff_21(2); 
	Vector<double> rho_coeff_22(2); 
	Vector<double> rho_coeff_23(2); 
	Vector<double> rho_coeff_24(2);
	Vector<double> b_rho2(2);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_RHO(no_stencils); Vector<double> w_RHO(no_stencils); double sum_RHO;
    
    // Variables for reconstruction of RHO_U
    double rho_u0; 

	/* Fourth order stencil */ 
	Vector<double> d_rho_u_4(4); 
	Vector<double> b_rho_u_4; 
	Vector<double> rho_u_coeff_4(9);
    
	/* Third order centered stencil */ 
	Vector<double> d_rho_u_3(4); 
	Vector<double> b_rho_u_3; 
	Vector<double> rho_u_coeff_3(5); 

	/* First third order one-sided stencil */ 
	Vector<double> b_rho_u_31;         
	Vector<double> d_rho_u_31;
	Vector<double> rho_u_coeff_31(5);
    
	/* Second third order one-sided stencil */
	Vector<double> b_rho_u_32;       
	Vector<double> d_rho_u_32;
	Vector<double> rho_u_coeff_32(5);
    
	/* Third third order one-sided stencil */
	Vector<double> b_rho_u_33;          
	Vector<double> d_rho_u_33;
	Vector<double> rho_u_coeff_33(5);
    
	/* Fourth third order one-sided stencil */
    Vector<double> b_rho_u_34; 
    Vector<double> d_rho_u_34;
    Vector<double> rho_u_coeff_34(5);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_RHO_U(no_stencils); Vector<double> w_RHO_U(no_stencils);  double sum_RHO_U;
	
	Vector<double> rho_u_coeff_21(2); 
    Vector<double> rho_u_coeff_22(2); 
    Vector<double> rho_u_coeff_23(2); 
    Vector<double> rho_u_coeff_24(2); 
    Vector<double> b_rho_u_2(2);
    
    // Variables for reconstruction of RHO_V
    double rho_v0;
    
    /* Fourth order stencil */ 
    Vector<double> d_rho_v_4(4); 
    Vector<double> b_rho_v_4; 
    Vector<double> rho_v_coeff_4(9); 
    
    /* Third order centered stencil */ 
    Vector<double> d_rho_v_3(4); 
    Vector<double> b_rho_v_3; 
    Vector<double> rho_v_coeff_3(5); 
	
	/* First third order one-sided stencil */ 
	Vector<double> b_rho_v_31;         
	Vector<double> d_rho_v_31;
	Vector<double> rho_v_coeff_31(5);
    
	/* Second third order one-sided stencil */
    Vector<double> b_rho_v_32;       
    Vector<double> d_rho_v_32;
    Vector<double> rho_v_coeff_32(5);
    
    /* Third third order one-sided stencil */
    Vector<double> b_rho_v_33;          
    Vector<double> d_rho_v_33;
    Vector<double> rho_v_coeff_33(5);
    
    /* Fourth third order one-sided stencil */
    Vector<double> b_rho_v_34; 
    Vector<double> d_rho_v_34;
    Vector<double> rho_v_coeff_34(5);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_RHO_V(no_stencils); Vector<double> w_RHO_V(no_stencils); double sum_RHO_V;
	
	Vector<double> rho_v_coeff_21(2); 
    Vector<double> rho_v_coeff_22(2); 
    Vector<double> rho_v_coeff_23(2); 
    Vector<double> rho_v_coeff_24(2); 
    Vector<double> b_rho_v_2(2);
    
    // Variables for reconstruction of E
    double e0; 

    /* Fourth order stencil */ 
    Vector<double> d_e_4(4); 
    Vector<double> b_e_4; 
    Vector<double> e_coeff_4(9); 
    
    /* Third order centered stencil */ 
    Vector<double> d_e_3(4); 
    Vector<double> b_e_3; 
    Vector<double> e_coeff_3(5); 
    
    /* First third order one-sided stencil */ 
    Vector<double> b_e_31;         
    Vector<double> d_e_31;
    Vector<double> e_coeff_31(5); 
    
    /* Second third order one-sided stencil */
    Vector<double> b_e_32;       
    Vector<double> d_e_32;
    Vector<double> e_coeff_32(5);
    
    /* Third third order one-sided stencil */
    Vector<double> b_e_33;          
    Vector<double> d_e_33;
    Vector<double> e_coeff_33(5);
    
    /* Fourth third order one-sided stencil */
    Vector<double> b_e_34; 
    Vector<double> d_e_34;
    Vector<double> e_coeff_34(5);
    
	/* Smoothness Indicators */ 
	Vector<double> IS_E(no_stencils); Vector<double> w_E(no_stencils); double sum_E;  
	
	Vector<double> e_coeff_21(2); 
	Vector<double> e_coeff_22(2); 
	Vector<double> e_coeff_23(2); 
	Vector<double> e_coeff_24(2);
	Vector<double> b_e_2(2);
    
    // Iterate over all the cells 
    
    DoFHandler<2>::active_cell_iterator cell, neighbor;
    
    unsigned int index, ROWS, g_i; 

	unsigned int WW_index, NN_index, EE_index, SS_index;

	unsigned int neighbor_p1, neighbor_p2, neighbor_p3, neighbor_p4;
	
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		g_i = local_to_global_index_map[c];
		cell = local_index_to_iterator[c];

	    bool WW = false, EE = false, NN = false, SS = false;

        rho0   =   RHO(g_i);
        rho_u0 =   RHO_U(g_i);
        rho_v0 =   RHO_V(g_i);
        e0     =   E(g_i);
        h = std::sqrt(Cell[c].measure()); 

		if(!cell->face(0)->at_boundary()) {
	        neighbor = cell->neighbor(0);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p1 = local_neighbor_dof_indices[0];
		}

		if(!cell->face(3)->at_boundary()) {
	        neighbor = cell->neighbor(3);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p2 = local_neighbor_dof_indices[0];
		}

		if(!cell->face(1)->at_boundary()) {
	        neighbor = cell->neighbor(1);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p3 = local_neighbor_dof_indices[0];
		}

		if(!cell->face(2)->at_boundary()) {
	        neighbor = cell->neighbor(2);
	        neighbor->get_dof_indices(local_neighbor_dof_indices);
			neighbor_p4 = local_neighbor_dof_indices[0];
		}

		if (cell_neighbor_neighbor_index[c][0].size() > 0) {
			WW = true;
			WW_index = cell_neighbor_neighbor_index[c][0][0] ;
		}

		if (cell_neighbor_neighbor_index[c][1].size() > 0) {
			EE = true;
			EE_index = cell_neighbor_neighbor_index[c][1][0] ;
		}

		if (cell_neighbor_neighbor_index[c][2].size() > 0) {
			SS = true;
			SS_index = cell_neighbor_neighbor_index[c][2][0] ;
		}

		if (cell_neighbor_neighbor_index[c][3].size() > 0) {
			NN = true;
			NN_index = cell_neighbor_neighbor_index[c][3][0] ;
		}
        
        if ( !(cell->at_boundary()) ) {
//			pcout<<"interior"<<std::endl;
			
            d_rho_31.reinit(4);   d_rho_32.reinit(4);   d_rho_33.reinit(4);   d_rho_34.reinit(4);
            d_rho_u_31.reinit(4); d_rho_u_32.reinit(4); d_rho_u_33.reinit(4); d_rho_u_34.reinit(4);
            d_rho_v_31.reinit(4); d_rho_v_32.reinit(4); d_rho_v_33.reinit(4); d_rho_v_34.reinit(4);
            d_e_31.reinit(4);     d_e_32.reinit(4);     d_e_33.reinit(4);     d_e_34.reinit(4);
            
            // =====================================================================
            // r = 4 stencil 
            // =====================================================================
            
			d_rho_4(0)  = (RHO(neighbor_p1) - rho0);   // W neighbor 
			d_rho_4(1)  = (RHO(neighbor_p2) - rho0);   // N neighbor
			d_rho_4(2)  = (RHO(neighbor_p3) - rho0);   // E neighbor
			d_rho_4(3)  = (RHO(neighbor_p4) - rho0);   // S neighbor
			
			d_rho_u_4(0)  = (RHO_U(neighbor_p1) - rho_u0);   // W neighbor 
			d_rho_u_4(1)  = (RHO_U(neighbor_p2) - rho_u0);    // N neighbor
			d_rho_u_4(2)  = (RHO_U(neighbor_p3) - rho_u0);    // E neighbor
			d_rho_u_4(3)  = (RHO_U(neighbor_p4) - rho_u0);    // S neighbor
			
			d_rho_v_4(0)  = (RHO_V(neighbor_p1) - rho_v0);   // W neighbor 
			d_rho_v_4(1)  = (RHO_V(neighbor_p2) - rho_v0);   // N neighbor
			d_rho_v_4(2)  = (RHO_V(neighbor_p3) - rho_v0);   // S neighbor
			d_rho_v_4(3)  = (RHO_V(neighbor_p4) - rho_v0);   // E neighbor
			
			d_e_4(0)  = (E(neighbor_p1) - e0);   // W neighbor 
			d_e_4(1)  = (E(neighbor_p2) - e0);   // N neighbor
			d_e_4(2)  = (E(neighbor_p3) - e0);   // S neighbor
			d_e_4(3)  = (E(neighbor_p4) - e0);   // E neighbor
            
            index = 0; 
            
            // Least Squares Part  
		
			ROWS = cell_diagonal_neighbor_index[c].size();
        			
			// neighbor of neighbors 

			if (WW) {
				ROWS++; 
			}

			if (EE) {
				ROWS++; 
			}

			if (SS) {
				ROWS++; 
			}

			if (NN) {
				ROWS++; 
			}
            
            b_rho_4.reinit(ROWS); b_rho_u_4.reinit(ROWS); b_rho_v_4.reinit(ROWS); b_e_4.reinit(ROWS);
			
			// vertex neighbors
			
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];

				b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
				b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
				b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
				b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
				index++; 
			}
			
			// neighbors of neighbors 
			
			if (WW) {

				local_neighbor_dof_indices[0] = WW_index;

				b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
				b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
				b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
				b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
				index++; 
			}
			
			if (NN) {

				local_neighbor_dof_indices[0] = NN_index;

				b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
				b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
				b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
				b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
				index++; 
			}
			
			if (EE) {

				local_neighbor_dof_indices[0] = EE_index;

				b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
				b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
				b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
				b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
				index++; 
			}
			
			if (SS) {

				local_neighbor_dof_indices[0] = SS_index;

				b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
				b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
				b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
				b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
				index++; 
			}

//			pcout<<"cls_R4 interior"<<std::endl;
			
			CLS_R4[c].solve(b_rho_4, d_rho_4, rho_coeff_4);  
			CLS_R4[c].solve(b_rho_u_4, d_rho_u_4, rho_u_coeff_4);
			CLS_R4[c].solve(b_rho_v_4, d_rho_v_4, rho_v_coeff_4);
			CLS_R4[c].solve(b_e_4, d_e_4, e_coeff_4);

//			if(c == 1) std::cout<<"b_e_4: "<<b_e_4<<std::endl<<"d_e_4: "<<d_e_4<<std::endl<<"e_coeff_4: "<<e_coeff_4<<std::endl;
//			if(c == 1) std::cout<<"b_e_4: "<<b_e_4(0)<<"\t"<<b_e_4(1)<<"\t"<<b_e_4(2)<<"\t"<<b_e_4(3)<<std::endl;
//			pcout<<"cls_R4 interior end"<<std::endl;
			
			// =====================================================================
            // r = 3 stencil (Centered stencil)
            // =====================================================================
			
			// constraint part (consists of face neighbours)
			
			d_rho_3(0)  = (RHO(neighbor_p1) - rho0);                // W neighbor 
			d_rho_3(1)  = (RHO(neighbor_p2) - rho0);                // N neighbor
			d_rho_3(2)  = (RHO(neighbor_p3) - rho0);                // E neighbor
			d_rho_3(3)  = (RHO(neighbor_p4) - rho0);                // S neighbor
			
			d_rho_u_3(0)  = (RHO_U(neighbor_p1) - rho_u0);          // W neighbor 
			d_rho_u_3(1)  = (RHO_U(neighbor_p2) - rho_u0);          // N neighbor
			d_rho_u_3(2)  = (RHO_U(neighbor_p3) - rho_u0);          // E neighbor
			d_rho_u_3(3)  = (RHO_U(neighbor_p4) - rho_u0);          // S neighbor
			
			d_rho_v_3(0)  = (RHO_V(neighbor_p1) - rho_v0);          // W neighbor 
			d_rho_v_3(1)  = (RHO_V(neighbor_p2) - rho_v0);          // N neighbor
			d_rho_v_3(2)  = (RHO_V(neighbor_p3) - rho_v0);          // E neighbor
			d_rho_v_3(3)  = (RHO_V(neighbor_p4) - rho_v0);          // S neighbor
			
			d_e_3(0)  = (E(neighbor_p1) - e0);                      // W neighbor 
			d_e_3(1)  = (E(neighbor_p2) - e0);                      // N neighbor
			d_e_3(2)  = (E(neighbor_p3) - e0);                      // E neighbor
			d_e_3(3)  = (E(neighbor_p4) - e0);                      // S neighbor
					
			ROWS = cell_diagonal_neighbor_index[c].size();
			index = 0; 
            
            b_rho_3.reinit(ROWS);
            b_rho_u_3.reinit(ROWS);
            b_rho_v_3.reinit(ROWS);
            b_e_3.reinit(ROWS);

			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {

				local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];

				b_rho_3(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
				b_rho_u_3(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
				b_rho_v_3(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
				b_e_3(index)     = E(local_neighbor_dof_indices[0] ) - e0;
				index++; 
			}

//			std::cout<<"cls_R3 interior"<<std::endl;			 
          
			CLS_R3[c].solve(b_rho_3, d_rho_3, rho_coeff_3);
			CLS_R3[c].solve(b_rho_u_3, d_rho_u_3, rho_u_coeff_3); 
			CLS_R3[c].solve(b_rho_v_3, d_rho_v_3, rho_v_coeff_3);
			CLS_R3[c].solve(b_e_3, d_e_3, e_coeff_3);

//			std::cout<<"cls_R3 interior end"<<std::endl;			 
			
			// =====================================================================
			// r = 3 stencil 1
			// =====================================================================
            
			if (is_admissible_R31[c]) {
            
				d_rho_31(0)  = RHO(neighbor_p1) - rho0;                // W neighbor 
				d_rho_31(1)  = RHO(neighbor_p2) - rho0;                // N neighbor
				d_rho_31(2)  = RHO(neighbor_p4) - rho0;                // S neighbor
				d_rho_31(3)  = RHO(WW_index) - rho0;               // WW neighbor
                
				d_rho_u_31(0)  = RHO_U(neighbor_p1) - rho_u0;          // W neighbor 
				d_rho_u_31(1)  = RHO_U(neighbor_p2) - rho_u0;          // N neighbor
				d_rho_u_31(2)  = RHO_U(neighbor_p4) - rho_u0;          // S neighbor
				d_rho_u_31(3)  = RHO_U(WW_index) - rho_u0;         // WW neighbor
                
				d_rho_v_31(0)  = RHO_V(neighbor_p1) - rho_v0;          // W neighbor 
				d_rho_v_31(1)  = RHO_V(neighbor_p2) - rho_v0;          // N neighbor
				d_rho_v_31(2)  = RHO_V(neighbor_p4) - rho_v0;          // S neighbor
				d_rho_v_31(3)  = RHO_V(WW_index) - rho_v0;         // WW neighbor
                
				d_e_31(0)  = E(neighbor_p1) - e0;                      // W neighbor 
				d_e_31(1)  = E(neighbor_p2) - e0;                      // N neighbor
				d_e_31(2)  = E(neighbor_p4) - e0;                      // S neighbor
				d_e_31(3)  = E(WW_index) - e0;                     // WW neighbor
				
                ROWS = cell_neighbor_index[c][0].size(); 
				
				b_rho_31.reinit(ROWS); b_rho_u_31.reinit(ROWS); b_rho_v_31.reinit(ROWS); b_e_31.reinit(ROWS);

                index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][0].size(); ++d) {
	
					local_neighbor_dof_indices[0] = cell_neighbor_index[c][0][d];

					b_rho_31(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_31(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_31(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_31(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
//			std::cout<<"cls_R31 interior"<<std::endl;			 				
                CLS_R31[c].solve(b_rho_31, d_rho_31, rho_coeff_31);
                CLS_R31[c].solve(b_rho_u_31, d_rho_u_31, rho_u_coeff_31);
                CLS_R31[c].solve(b_rho_v_31, d_rho_v_31, rho_v_coeff_31);
                CLS_R31[c].solve(b_e_31, d_e_31, e_coeff_31); 
//			pcout<<"cls_R31 interior end"<<std::endl;			 				
			}
            
            // If stencil is not available fallback to third order centered stencil 
            
            else {
				rho_coeff_31 = rho_coeff_3; 
				rho_u_coeff_31 = rho_u_coeff_3;
				rho_v_coeff_31 = rho_v_coeff_3;
				e_coeff_31 = e_coeff_3;
			}
            
            // =====================================================================
            // r = 3 stencil 2
            // =====================================================================
            
			if (is_admissible_R32[c]) {
            
				d_rho_32(0)  = RHO(neighbor_p1) - rho0;                // W neighbor 
				d_rho_32(1)  = RHO(neighbor_p2) - rho0;                // N neighbor
				d_rho_32(2)  = RHO(neighbor_p3) - rho0;                // E neighbor
				d_rho_32(3)  = RHO(NN_index) - rho0;               // NN neighbor
                
				d_rho_u_32(0)  = RHO_U(neighbor_p1) - rho_u0;          // W neighbor 
				d_rho_u_32(1)  = RHO_U(neighbor_p2) - rho_u0;          // N neighbor
				d_rho_u_32(2)  = RHO_U(neighbor_p3) - rho_u0;          // E neighbor
				d_rho_u_32(3)  = RHO_U(NN_index) - rho_u0;           // NN neighbor
                
				d_rho_v_32(0)  = RHO_V(neighbor_p1) - rho_v0;          // W neighbor 
				d_rho_v_32(1)  = RHO_V(neighbor_p2) - rho_v0;          // N neighbor
				d_rho_v_32(2)  = RHO_V(neighbor_p3) - rho_v0;          // E neighbor
				d_rho_v_32(3)  = RHO_V(NN_index) - rho_v0;           // NN neighbor
                
				d_e_32(0)  = E(neighbor_p1) - e0;                      // W neighbor 
				d_e_32(1)  = E(neighbor_p2) - e0;                      // N neighbor
				d_e_32(2)  = E(neighbor_p3) - e0;                      // E neighbor
				d_e_32(3)  = E(NN_index) - e0;                     // NN neighbor

                ROWS = cell_neighbor_index[c][3].size(); 

				b_rho_32.reinit(ROWS); b_rho_u_32.reinit(ROWS); b_rho_v_32.reinit(ROWS); b_e_32.reinit(ROWS);

				index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][3].size(); ++d) {
	
					local_neighbor_dof_indices[0] = cell_neighbor_index[c][3][d];

					b_rho_32(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_32(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_32(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_32(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
//			std::cout<<"cls_R32 interior"<<std::endl;			 				
                CLS_R32[c].solve(b_rho_32, d_rho_32, rho_coeff_32);
                CLS_R32[c].solve(b_rho_u_32, d_rho_u_32, rho_u_coeff_32);
                CLS_R32[c].solve(b_rho_v_32, d_rho_v_32, rho_v_coeff_32);
                CLS_R32[c].solve(b_e_32, d_e_32, e_coeff_32); 
//			pcout<<"cls_R32 interior end"<<std::endl;			 				
			}
            // If stencil is not available fallback to third order centered stencil
            
            else {
				rho_coeff_32 = rho_coeff_3; 
				rho_u_coeff_32 = rho_u_coeff_3;
				rho_v_coeff_32 = rho_v_coeff_3;
				e_coeff_32 = e_coeff_3;
            }  
            
            // =====================================================================
            // r = 3 stencil 3
            // =====================================================================
            
			if (is_admissible_R33[c]) {
            
				d_rho_33(0)  = RHO(neighbor_p1) - rho0;                // W neighbor 
				d_rho_33(1)  = RHO(neighbor_p3) - rho0;                // E neighbor
				d_rho_33(2)  = RHO(neighbor_p4) - rho0;                // S neighbor
				d_rho_33(3)  = RHO(SS_index) - rho0;               // SS neighbor
                
				d_rho_u_33(0)  = RHO_U(neighbor_p1) - rho_u0;          // W neighbor 
				d_rho_u_33(1)  = RHO_U(neighbor_p3) - rho_u0;          // E neighbor
				d_rho_u_33(2)  = RHO_U(neighbor_p4) - rho_u0;          // S neighbor
				d_rho_u_33(3)  = RHO_U(SS_index) - rho_u0;         // SS neighbor
                
				d_rho_v_33(0)  = RHO_V(neighbor_p1) - rho_v0;          // W neighbor 
				d_rho_v_33(1)  = RHO_V(neighbor_p3) - rho_v0;          // E neighbor
				d_rho_v_33(2)  = RHO_V(neighbor_p4) - rho_v0;          // S neighbor
				d_rho_v_33(3)  = RHO_V(SS_index) - rho_v0;         // SS neighbor
                
				d_e_33(0)  = E(neighbor_p1) - e0;                      // W neighbor 
				d_e_33(1)  = E(neighbor_p3) - e0;                      // E neighbor
				d_e_33(2)  = E(neighbor_p4) - e0;                      // S neighbor
				d_e_33(3)  = E(SS_index) - e0;                     // SS neighbor

                ROWS = cell_neighbor_index[c][2].size(); 

				b_rho_33.reinit(ROWS); b_rho_u_33.reinit(ROWS); b_rho_v_33.reinit(ROWS); b_e_33.reinit(ROWS);

                index = 0; 			

				for (unsigned int d = 0; d < cell_neighbor_index[c][2].size(); ++d) {
	
					local_neighbor_dof_indices[0] = cell_neighbor_index[c][2][d];

					b_rho_33(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_33(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_33(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_33(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
//			std::cout<<"cls_R33 interior"<<std::endl;			 				
                CLS_R33[c].solve(b_rho_33,   d_rho_33,   rho_coeff_33);
                CLS_R33[c].solve(b_rho_u_33, d_rho_u_33, rho_u_coeff_33);
                CLS_R33[c].solve(b_rho_v_33, d_rho_v_33, rho_v_coeff_33);
                CLS_R33[c].solve(b_e_33,     d_e_33,     e_coeff_33); 
//			pcout<<"cls_R33 interior end"<<std::endl;			 				
			}
            // If stencil is not available fallback to third order centered stencil
            
            else {
				rho_coeff_33 = rho_coeff_3; 
				rho_u_coeff_33 = rho_u_coeff_3;
				rho_v_coeff_33 = rho_v_coeff_3;
				e_coeff_33 = e_coeff_3;
			} 
            
			// =====================================================================
            // r = 3 stencil 4
            // =====================================================================
            
			if (is_admissible_R34[c]) {
            
				d_rho_34(0)  = RHO(neighbor_p2) - rho0;                // N neighbor 
				d_rho_34(1)  = RHO(neighbor_p3) - rho0;                // E neighbor
				d_rho_34(2)  = RHO(neighbor_p4) - rho0;                // S neighbor
				d_rho_34(3)  = RHO(EE_index) - rho0;                // S neighbor
                
				d_rho_u_34(0)  = RHO_U(neighbor_p2) - rho_u0;          // N neighbor 
				d_rho_u_34(1)  = RHO_U(neighbor_p3) - rho_u0;          // E neighbor
				d_rho_u_34(2)  = RHO_U(neighbor_p4) - rho_u0;          // S neighbor
				d_rho_u_34(3)  = RHO_U(EE_index) - rho_u0;          // S neighbor
                
				d_rho_v_34(0)  = RHO_V(neighbor_p2) - rho_v0;          // N neighbor 
				d_rho_v_34(1)  = RHO_V(neighbor_p3) - rho_v0;          // E neighbor
				d_rho_v_34(2)  = RHO_V(neighbor_p4) - rho_v0;          // S neighbor
				d_rho_v_34(3)  = RHO_V(EE_index) - rho_v0;          // S neighbor
                
				d_e_34(0)  = E(neighbor_p2) - e0;                      // N neighbor 
				d_e_34(1)  = E(neighbor_p3) - e0;                      // E neighbor
				d_e_34(2)  = E(neighbor_p4) - e0;                      // S neighbor
				d_e_34(3)  = E(EE_index) - e0;                      // S neighbor

                ROWS = cell_neighbor_index[c][1].size(); 

                b_rho_34.reinit(ROWS); b_rho_u_34.reinit(ROWS); b_rho_v_34.reinit(ROWS); b_e_34.reinit(ROWS);

                index = 0; 

				for (unsigned int d = 0; d < cell_neighbor_index[c][1].size(); ++d) {
	
					local_neighbor_dof_indices[0] = cell_neighbor_index[c][1][d];

					b_rho_34(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_34(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_34(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_34(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
//			if (c == 241) pcout<<"cls_R34 interior rho"<<"\tc: "<<c<<std::endl<<b_rho_34<<std::endl<<d_rho_34<<std::endl;			 				
				CLS_R34[c].solve(b_rho_34, d_rho_34, rho_coeff_34);
	//		if (c == 241) pcout<<"cls_R34 interior rho u"<<"\tc: "<<c<<std::endl;			 				
				CLS_R34[c].solve(b_rho_u_34, d_rho_u_34, rho_u_coeff_34);
		//	if (c == 241) pcout<<"cls_R34 interior rho v"<<"\tc: "<<c<<std::endl;			 				
				CLS_R34[c].solve(b_rho_v_34, d_rho_v_34, rho_v_coeff_34);
	//		if (c == 241) pcout<<"cls_R34 interior e"<<"\tc: "<<c<<std::endl;			 				
				CLS_R34[c].solve(b_e_34, d_e_34, e_coeff_34); 
		//	if (c == 241) pcout<<"cls_R34 interior end"<<std::endl;			 				
            }
            
            // If stencil is not available fallback to third order centered stencil
            
            else {
				rho_coeff_34 = rho_coeff_3; 
				rho_u_coeff_34 = rho_u_coeff_3;
				rho_v_coeff_34 = rho_v_coeff_3;
				e_coeff_34 = e_coeff_3;
            } 
            
            
			// =====================================================================
			// r = 2 stencil 1
			// ===================================================================== 

			b_rho2(0)  = (RHO(neighbor_p1) - rho0);          // W neighbor 
			b_rho2(1)  = (RHO(neighbor_p2) - rho0);          // N neighbor
			LU_R21[c].solve(b_rho2, rho_coeff_21);

			b_rho_u_2(0)  = (RHO_U(neighbor_p1) - rho_u0);   // W neighbor 
			b_rho_u_2(1)  = (RHO_U(neighbor_p2) - rho_u0);   // N neighbor
			LU_R21[c].solve(b_rho_u_2, rho_u_coeff_21);

			b_rho_v_2(0)  = (RHO_V(neighbor_p1) - rho_v0);   // W neighbor 
			b_rho_v_2(1)  = (RHO_V(neighbor_p2) - rho_v0);   // N neighbor
			LU_R21[c].solve(b_rho_v_2, rho_v_coeff_21);

			b_e_2(0)  = (E(neighbor_p1) - e0);               // W neighbor 
			b_e_2(1)  = (E(neighbor_p2) - e0);               // N neighbor
			LU_R21[c].solve(b_e_2, e_coeff_21);
            
            // =====================================================================
            // r = 2 stencil 2 
            // ===================================================================== 

			b_rho2(0) = (RHO(neighbor_p2) - rho0);          // N neighbor
			b_rho2(1) = (RHO(neighbor_p3) - rho0);          // E neighbor
			LU_R22[c].solve(b_rho2, rho_coeff_22);

			b_rho_u_2(0) = (RHO_U(neighbor_p2) - rho_u0);   // N neighbor
			b_rho_u_2(1) = (RHO_U(neighbor_p3) - rho_u0);   // E neighbor
			LU_R22[c].solve(b_rho_u_2, rho_u_coeff_22);

			b_rho_v_2(0) = (RHO_V(neighbor_p2) - rho_v0);   // N neighbor
			b_rho_v_2(1) = (RHO_V(neighbor_p3) - rho_v0);   // E neighbor
			LU_R22[c].solve(b_rho_v_2, rho_v_coeff_22);

			b_e_2(0) = (E(neighbor_p2) - e0);               // N neighbor
			b_e_2(1) = (E(neighbor_p3) - e0);               // E neighbor
			LU_R22[c].solve(b_e_2, e_coeff_22);
		
			// =====================================================================
			// r = 2 stencil 3
			// =====================================================================

			b_rho2(0) = (RHO(neighbor_p3) - rho0);          // E neighbor
			b_rho2(1) = (RHO(neighbor_p4) - rho0);          // S neighbor                                          
			LU_R23[c].solve(b_rho2, rho_coeff_23); 

			b_rho_u_2(0) = (RHO_U(neighbor_p3) - rho_u0);   // E neighbor
			b_rho_u_2(1) = (RHO_U(neighbor_p4) - rho_u0);   // S neighbor                                   
			LU_R23[c].solve(b_rho_u_2, rho_u_coeff_23); 

			b_rho_v_2(0) = (RHO_V(neighbor_p3) - rho_v0);   // E neighbor
			b_rho_v_2(1) = (RHO_V(neighbor_p4) - rho_v0);   // S neighbor                                   
			LU_R23[c].solve(b_rho_v_2, rho_v_coeff_23);

			b_e_2(0) = (E(neighbor_p3) - e0);               // E neighbor
			b_e_2(1) = (E(neighbor_p4) - e0);               // S neighbor                                              
			LU_R23[c].solve(b_e_2, e_coeff_23);  
			
			// =====================================================================
			// r = 2 stencil 4
			// =====================================================================

			b_rho2(0) = (RHO(neighbor_p4) - rho0);          // S neighbor
			b_rho2(1) = (RHO(neighbor_p1) - rho0);          // W neighbor
			LU_R24[c].solve(b_rho2, rho_coeff_24);

			b_rho_u_2(0) = (RHO_U(neighbor_p4) - rho_u0);   // S neighbor
			b_rho_u_2(1) = (RHO_U(neighbor_p1) - rho_u0);   // W neighbor
			LU_R24[c].solve(b_rho_u_2, rho_u_coeff_24);

			b_rho_v_2(0) = (RHO_V(neighbor_p4) - rho_v0);   // S neighbor
			b_rho_v_2(1) = (RHO_V(neighbor_p1) - rho_v0);   // W neighbor
			LU_R24[c].solve(b_rho_v_2, rho_v_coeff_24);

			b_e_2(0) = (E(neighbor_p4) - e0);               // S neighbor
			b_e_2(1) = (E(neighbor_p1) - e0);               // W neighbor
			LU_R24[c].solve(b_e_2, e_coeff_24);
 //  			pcout<<"interior end"<<std::endl;    
        } // End of interior cell loop 
        
        
        else {

            bool W_face = false, E_face = false, N_face = false, S_face = false;
            
            if (!(is_corner_cell[c])) {

//	   			pcout<<"boundary "<<std::endl;    
                
				if(cell->face(0)->at_boundary()) { W_face = true;}
				if(cell->face(1)->at_boundary()) { E_face = true;}
				if(cell->face(2)->at_boundary()) { S_face = true;}
				if(cell->face(3)->at_boundary()) { N_face = true;}
                
                // =====================================================================
                // r = 4 stencil (boundary)
                // =====================================================================
            
                Vector<double> d_rho; Vector<double> d_rho_u; Vector<double> d_rho_v; Vector<double> d_e;
                
                d_rho.reinit(7); d_rho_u.reinit(7); d_rho_v.reinit(7); d_e.reinit(7);
                
                index = 0;
                
                if (!W_face) {
                    
                        d_rho(index)   = RHO(neighbor_p1) - rho0;
                        d_rho_u(index) = RHO_U(neighbor_p1) - rho_u0;
                        d_rho_v(index) = RHO_V(neighbor_p1) - rho_v0;
                        d_e(index)     = E(neighbor_p1) - e0;
                        
                        index++; 
                }
                
                if (!N_face) {
                    
                    d_rho(index)   = RHO(neighbor_p2) - rho0;
                    d_rho_u(index) = RHO_U(neighbor_p2) - rho_u0;
                    d_rho_v(index) = RHO_V(neighbor_p2) - rho_v0;
                    d_e(index)     = E(neighbor_p2) - e0;
                    
                    index++; 
                    
                }
                
                if (!E_face) {
                    
                    d_rho(index)   = RHO(neighbor_p3) - rho0;
                    d_rho_u(index) = RHO_U(neighbor_p3) - rho_u0;
                    d_rho_v(index) = RHO_V(neighbor_p3) - rho_v0;
                    d_e(index)     = E(neighbor_p3) - e0;
                    
                    index++; 
                    
                }
                
                if (!S_face) {
                    
                    d_rho(index)   = RHO(neighbor_p4) - rho0;
                    d_rho_u(index) = RHO_U(neighbor_p4) - rho_u0;
                    d_rho_v(index) = RHO_V(neighbor_p4) - rho_v0;
                    d_e(index)     = E(neighbor_p4) - e0;
                    
                    index++; 
                    
                }
                
                d_rho(3)   = 0.0; d_rho(4)   = 0.0; d_rho(5)   = 0.0; d_rho(6)   = 0.0;
                d_rho_u(3) = 0.0; d_rho_u(4) = 0.0; d_rho_u(5) = 0.0; d_rho_u(6) = 0.0;
                d_rho_v(3) = 0.0; d_rho_v(4) = 0.0; d_rho_v(5) = 0.0; d_rho_v(6) = 0.0;
                d_e(3)     = 0.0; d_e(4)     = 0.0; d_e(5)     = 0.0; d_e(6)     = 0.0;
                
                index = 0;
                
				// Least Squares Part  
			
				ROWS = cell_diagonal_neighbor_index[c].size(); 
				
				// neighbor of neighbors 
				
				if (WW) {
					ROWS++; 
				}
	
				if (EE) {
					ROWS++; 
				}
		
				if (SS) {
					ROWS++; 
				}

				if (NN) {
					ROWS++; 
				}
				
				b_rho_4.reinit(ROWS); b_rho_u_4.reinit(ROWS); b_rho_v_4.reinit(ROWS); b_e_4.reinit(ROWS);
				
				// vertex neighbors
				
				for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			

					local_neighbor_dof_indices[0] = cell_diagonal_neighbor_index[c][d];	

					b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
								
				// neighbors of neighbors 
				
				if (WW) {

					local_neighbor_dof_indices[0] = WW_index;

					b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
				
				if (NN) {

					local_neighbor_dof_indices[0] = NN_index;

					b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
				
				if (EE) {

					local_neighbor_dof_indices[0] = EE_index;

					b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
				
				if (SS) {

					local_neighbor_dof_indices[0] = SS_index;

					b_rho_4(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
					b_rho_u_4(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
					b_rho_v_4(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
					b_e_4(index)     = E(local_neighbor_dof_indices[0] ) - e0;
					index++; 
				}
                
                
                CLS_R4[c].solve(b_rho_4, d_rho, rho_coeff_4);
                CLS_R4[c].solve(b_rho_u_4, d_rho_u, rho_u_coeff_4);
                CLS_R4[c].solve(b_rho_v_4, d_rho_v, rho_v_coeff_4);
                CLS_R4[c].solve(b_e_4, d_e, e_coeff_4);

//				if(c == 1) std::cout<<"b_e_4: "<<b_e_4<<std::endl<<"d_e_4: "<<d_e_4<<std::endl<<"e_coeff_4: "<<e_coeff_4<<std::endl;
//				if(c == 1) std::cout<<"b_e_4: "<<b_e_4(0)<<"\t"<<b_e_4(1)<<"\t"<<b_e_4(2)<<"\t"<<b_e_4(3)<<std::endl;                
//				if(c == 1) std::cout<<"d_e: "<<d_e(0)<<"\t"<<d_e(1)<<"\t"<<d_e(2)<<"\t"<<d_e(3)<<std::endl;
//				if(c == 1) std::cout<<"e_coeff_4: "<<e_coeff_4(0)<<"\t"<<e_coeff_4(1)<<"\t"<<e_coeff_4(2)<<"\t"<<e_coeff_4(3)<<std::endl;                
				// =====================================================================
                // r = 3 center stencil (boundary)
                // =====================================================================
                
                b_rho_3.reinit(5); b_rho_u_3.reinit(5); b_rho_v_3.reinit(5); b_e_3.reinit(5);  
                
                index = 0; 
                
                if (!W_face) {
                    
                    // P1 neighbor
                    
                    b_rho_3(index)    = (RHO(neighbor_p1) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(neighbor_p1) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(neighbor_p1) - rho_v0);                
                    b_e_3(index)      = (E(neighbor_p1) - e0);                
                    
                    index++; 
                }
                
                if (!N_face) {

                    // P2 neighbor
                    
                    b_rho_3(index)    = (RHO(neighbor_p2) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(neighbor_p2) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(neighbor_p2) - rho_v0);                
                    b_e_3(index)      = (E(neighbor_p2) - e0);                
                    
                    index++; 
                
                }
                
                if (!E_face) {

                    // P3 neighbor
                    
                    b_rho_3(index)    = (RHO(neighbor_p3) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(neighbor_p3) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(neighbor_p3) - rho_v0);                
                    b_e_3(index)      = (E(neighbor_p3) - e0);                
                    
                    index++; 
                }
                
                if (!S_face) {

                    // P4 neighbor
                    
                    b_rho_3(index)    = (RHO(neighbor_p4) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(neighbor_p4) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(neighbor_p4) - rho_v0);                
                    b_e_3(index)      = (E(neighbor_p4) - e0);                
                    
                    index++; 
                }
                
                
                // Transmissive boundary conditions 
                
				b_rho_3(3) = 0.0; b_rho_u_3(3) = 0.0; b_rho_v_3(3) = 0.0; b_e_3(3) = 0.0;
				b_rho_3(4) = 0.0; b_rho_u_3(4) = 0.0; b_rho_v_3(4) = 0.0; b_e_3(4) = 0.0;
                
            
				CLS_R3[c].solve(b_rho_3, rho_coeff_3); 
				CLS_R3[c].solve(b_rho_u_3, rho_u_coeff_3);
				CLS_R3[c].solve(b_rho_v_3, rho_v_coeff_3);
				CLS_R3[c].solve(b_e_3, e_coeff_3);

//				if(c == 1) std::cout<<"b_e_3: "<<b_e_3(0)<<"\t"<<b_e_3(1)<<"\t"<<b_e_3(2)<<"\t"<<b_e_3(3)<<std::endl;						
//				if(c == 1) std::cout<<"e_coeff_3: "<<e_coeff_3(0)<<"\t"<<e_coeff_3(1)<<"\t"<<e_coeff_3(2)<<"\t"<<e_coeff_3(3)<<std::endl;						                                
				// =====================================================================
                // r = 3 stencil 1 (boundary)
                // =====================================================================
                
                if (is_admissible_R31[c]) {
                    
                    index = 0; 
                    if ( N_face || S_face ) {                    
                        d_rho_31.reinit(5); d_rho_u_31.reinit(5); d_rho_v_31.reinit(5); d_e_31.reinit(5);  
                    }
            
                    else {
                        d_rho_31.reinit(4); d_rho_u_31.reinit(4); d_rho_v_31.reinit(4); d_e_31.reinit(4);
                    }
                    
                    if (!W_face) {
                    
                        d_rho_31(index)  = RHO(neighbor_p1) - rho0;                // P1 neighbor 
                        d_rho_u_31(index)  = RHO_U(neighbor_p1) - rho_u0;          // P1 neighbor 
                        d_rho_v_31(index)  = RHO_V(neighbor_p1) - rho_v0;          // P1 neighbor 
                        d_e_31(index)  = E(neighbor_p1) - e0;                      // P1 neighbor 
                        
                        index++;
                    }
                    
                    if (!N_face) {
                    
                        d_rho_31(index)  = RHO(neighbor_p2) - rho0;                // P2 neighbor 
                        d_rho_u_31(index)  = RHO_U(neighbor_p2) - rho_u0;          // P2 neighbor 
                        d_rho_v_31(index)  = RHO_V(neighbor_p2) - rho_v0;          // P2 neighbor 
                        d_e_31(index)  = E(neighbor_p2) - e0;                      // P2 neighbor 
                        
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_rho_31(index)  = RHO(neighbor_p4) - rho0;                // P4 neighbor 
                        d_rho_u_31(index)  = RHO_U(neighbor_p4) - rho_u0;          // P4 neighbor 
                        d_rho_v_31(index)  = RHO_V(neighbor_p4) - rho_v0;          // P4 neighbor 
                        d_e_31(index)  = E(neighbor_p4) - e0;                      // P4 neighbor 
                        
                        index++;
                    }
                    
                    if (WW) {

						local_neighbor_dof_indices[0] = WW_index;
                    
                        d_rho_31(index)    = RHO(local_neighbor_dof_indices[0] ) - rho0;              // S1 neighbor 
                        d_rho_u_31(index)  = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;          // S1 neighbor 
                        d_rho_v_31(index)  = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;          // S1 neighbor 
                        d_e_31(index)      = E(local_neighbor_dof_indices[0] ) - e0;                  // S1 neighbor 
                        
                        index++;
                    }
                    if ( S_face || N_face ) {                    
                        
                        // Transmissive boundary conditions 
                        
                        d_rho_31(index)    = 0.0;   d_rho_31(index+1)    = 0.0;            
                        d_rho_u_31(index)  = 0.0;   d_rho_u_31(index+1)  = 0.0;      
                        d_rho_v_31(index)  = 0.0;   d_rho_v_31(index+1)  = 0.0;    
                        d_e_31(index)      = 0.0;   d_e_31(index+1)  = 0.0;
                    }
/*
					if(c == 2){

						std::cout<<"d: "<<std::endl;
						for(unsigned int i = 0; i < index+2; ++i )
							std::cout<<d_rho_31(i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < index+2; ++i )
							std::cout<<d_e_31(i)<<"\t";
						std::cout<<std::endl;

					}
*/                    
	                ROWS = cell_neighbor_index[c][0].size(); 
				
					b_rho_31.reinit(ROWS); b_rho_u_31.reinit(ROWS); b_rho_v_31.reinit(ROWS); b_e_31.reinit(ROWS);

					index = 0; 
					
					// vertex neighbor of cell at face 0 
					for (unsigned int d = 0; d < cell_neighbor_index[c][0].size(); ++d) {						

						local_neighbor_dof_indices[0] = cell_neighbor_index[c][0][d];

						b_rho_31(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
						b_rho_u_31(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
						b_rho_v_31(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
						b_e_31(index)     = E(local_neighbor_dof_indices[0] ) - e0;
						index++; 
					}

					CLS_R31[c].solve(b_rho_31, d_rho_31, rho_coeff_31);
					CLS_R31[c].solve(b_rho_u_31, d_rho_u_31, rho_u_coeff_31);
					CLS_R31[c].solve(b_rho_v_31, d_rho_v_31, rho_v_coeff_31);
					CLS_R31[c].solve(b_e_31, d_e_31, e_coeff_31);
/*
					if(c == 2){

						std::cout<<"b: "<<std::endl;
						for(unsigned int i = 0; i < index; ++i )
							std::cout<<b_rho_31(i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < index; ++i )
							std::cout<<b_e_31(i)<<"\t";
						std::cout<<std::endl;

						std::cout<<"coeff: "<<std::endl;
						for(unsigned int i = 0; i < 5; ++i )
							std::cout<<rho_coeff_31(i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 5; ++i )
							std::cout<<rho_u_coeff_31(i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 5; ++i )
							std::cout<<rho_v_coeff_31(i)<<"\t";
						std::cout<<std::endl;
						for(unsigned int i = 0; i < 5; ++i )
							std::cout<<e_coeff_31(i)<<"\t";
						std::cout<<std::endl;
					}
*/					
				}
                
                else {
					rho_coeff_31 = rho_coeff_3; 
					rho_u_coeff_31 = rho_u_coeff_3;
					rho_v_coeff_31 = rho_v_coeff_3;
					e_coeff_31 = e_coeff_3;
                }
                
				// =====================================================================
                // r = 3 stencil 2 (boundary)
                // =====================================================================
                
                if (is_admissible_R32[c]) {
                    
                    index = 0; 
                    if ( W_face || E_face ) {                    
                        d_rho_32.reinit(5); d_rho_u_32.reinit(5); d_rho_v_32.reinit(5); d_e_32.reinit(5);  
                    }
            
                    else {
                        d_rho_32.reinit(4); d_rho_u_32.reinit(4); d_rho_v_32.reinit(4); d_e_32.reinit(4);
                    }
                    
                    if (!W_face) {
                    
                        d_rho_32(index)  = RHO(neighbor_p1) - rho0;                // W neighbor 
                        d_rho_u_32(index)  = RHO_U(neighbor_p1) - rho_u0;          // W neighbor 
                        d_rho_v_32(index)  = RHO_V(neighbor_p1) - rho_v0;          // W neighbor 
                        d_e_32(index)  = E(neighbor_p1) - e0;                      // W neighbor 
                        
                        index++;
                    }
                    
                    if (!N_face) {
                    
                        d_rho_32(index)  = RHO(neighbor_p2) - rho0;                // N neighbor 
                        d_rho_u_32(index)  = RHO_U(neighbor_p2) - rho_u0;          // N neighbor 
                        d_rho_v_32(index)  = RHO_V(neighbor_p2) - rho_v0;          // N neighbor 
                        d_e_32(index)  = E(neighbor_p2) - e0;                      // N neighbor 
                        
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_rho_32(index)  = RHO(neighbor_p3) - rho0;                // P3 neighbor 
                        d_rho_u_32(index)  = RHO_U(neighbor_p3) - rho_u0;          // P3 neighbor 
                        d_rho_v_32(index)  = RHO_V(neighbor_p3) - rho_v0;          // P3 neighbor 
                        d_e_32(index)  = E(neighbor_p3) - e0;                      // P3 neighbor 
                        
                        index++;
                    }
                    
                    if (NN) {

						local_neighbor_dof_indices[0] = NN_index;
                     
                        d_rho_32(index)    = RHO(local_neighbor_dof_indices[0] ) - rho0;              // S3 neighbor 
                        d_rho_u_32(index)  = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;          // S3 neighbor 
                        d_rho_v_32(index)  = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;          // S3 neighbor 
                        d_e_32(index)      = E(local_neighbor_dof_indices[0] ) - e0;                  // S3 neighbor 
                        
                        index++;
                    }
                    
                    if ( W_face || E_face ) {
                        
                        // Transmissive boundary conditions 
                        
                        d_rho_32(index)    = 0.0;   d_rho_32(index+1)    = 0.0;            
                        d_rho_u_32(index)  = 0.0;   d_rho_u_32(index+1)  = 0.0;      
                        d_rho_v_32(index)  = 0.0;   d_rho_v_32(index+1)  = 0.0;    
                        d_e_32(index)      = 0.0;   d_e_32(index+1)  = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][3].size(); 

					b_rho_32.reinit(ROWS); b_rho_u_32.reinit(ROWS); b_rho_v_32.reinit(ROWS); b_e_32.reinit(ROWS);

					index = 0; 		

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][3].size(); ++d) {						

						local_neighbor_dof_indices[0] = cell_neighbor_index[c][3][d];

						b_rho_32(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
						b_rho_u_32(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
						b_rho_v_32(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
						b_e_32(index)     = E(local_neighbor_dof_indices[0] ) - e0;
						index++; 
					}
                    
					CLS_R32[c].solve(b_rho_32, d_rho_32, rho_coeff_32);
					CLS_R32[c].solve(b_rho_u_32, d_rho_u_32, rho_u_coeff_32);
					CLS_R32[c].solve(b_rho_v_32, d_rho_v_32, rho_v_coeff_32);
					CLS_R32[c].solve(b_e_32, d_e_32, e_coeff_32); 
					
				}
                
                else {
					rho_coeff_32 = rho_coeff_3; 
					rho_u_coeff_32 = rho_u_coeff_3;
					rho_v_coeff_32 = rho_v_coeff_3;
					e_coeff_32 = e_coeff_3;
				}
                
                // =====================================================================
                // r = 3 stencil 3 (boundary)
                // =====================================================================
                
                if (is_admissible_R33[c]) {
                    
                    index = 0; 
                    if ( W_face || E_face ) {                    
                        d_rho_33.reinit(5); d_rho_u_33.reinit(5); d_rho_v_33.reinit(5); d_e_33.reinit(5);  
                    }
            
                    else {
                        d_rho_33.reinit(4); d_rho_u_33.reinit(4); d_rho_v_33.reinit(4); d_e_33.reinit(4);
                    }
                    
                    if (!W_face) {
                    
                        d_rho_33(index)  = RHO(neighbor_p1) - rho0;                // W neighbor 
                        d_rho_u_33(index)  = RHO_U(neighbor_p1) - rho_u0;          // W neighbor 
                        d_rho_v_33(index)  = RHO_V(neighbor_p1) - rho_v0;          // W neighbor 
                        d_e_33(index)  = E(neighbor_p1) - e0;                      // W neighbor 
                        
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_rho_33(index)  = RHO(neighbor_p3) - rho0;                // E neighbor 
                        d_rho_u_33(index)  = RHO_U(neighbor_p3) - rho_u0;          // E neighbor 
                        d_rho_v_33(index)  = RHO_V(neighbor_p3) - rho_v0;          // E neighbor 
                        d_e_33(index)  = E(neighbor_p3) - e0;                      // E neighbor 
                        
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_rho_33(index)  = RHO(neighbor_p4) - rho0;                // S neighbor 
                        d_rho_u_33(index)  = RHO_U(neighbor_p4) - rho_u0;          // S neighbor 
                        d_rho_v_33(index)  = RHO_V(neighbor_p4) - rho_v0;          // S neighbor 
                        d_e_33(index)  = E(neighbor_p4) - e0;                      // S neighbor 
                        
                        index++;
                    }
                    
                    if (SS) {

						local_neighbor_dof_indices[0] = SS_index;
                    
                        d_rho_33(index)    = RHO(local_neighbor_dof_indices[0] ) - rho0;              // SS neighbor 
                        d_rho_u_33(index)  = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;          // SS neighbor 
                        d_rho_v_33(index)  = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;          // SS neighbor 
                        d_e_33(index)      = E(local_neighbor_dof_indices[0] ) - e0;                  // SS neighbor 
                        
                        index++;
                    }
                    
                    if ( W_face || E_face ) {                    
                        
                        // Transmissive boundary conditions 
                        
                        d_rho_33(index)    = 0.0;   d_rho_33(index+1)    = 0.0;            
                        d_rho_u_33(index)  = 0.0;   d_rho_u_33(index+1)  = 0.0;      
                        d_rho_v_33(index)  = 0.0;   d_rho_v_33(index+1)  = 0.0;    
                        d_e_33(index)      = 0.0;   d_e_33(index+1)  = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][2].size(); 

					b_rho_33.reinit(ROWS); b_rho_u_33.reinit(ROWS); b_rho_v_33.reinit(ROWS); b_e_33.reinit(ROWS);

					index = 0; 							

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][2].size(); ++d) {						

						local_neighbor_dof_indices[0] = cell_neighbor_index[c][2][d];

						b_rho_33(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
						b_rho_u_33(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
						b_rho_v_33(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
						b_e_33(index)     = E(local_neighbor_dof_indices[0] ) - e0;
						index++; 
					}

					CLS_R33[c].solve(b_rho_33, d_rho_33, rho_coeff_33);
					CLS_R33[c].solve(b_rho_u_33, d_rho_u_33, rho_u_coeff_33);
					CLS_R33[c].solve(b_rho_v_33, d_rho_v_33, rho_v_coeff_33);
					CLS_R33[c].solve(b_e_33, d_e_33, e_coeff_33); 
				}
                
                else {
					rho_coeff_33 = rho_coeff_3; 
					rho_u_coeff_33 = rho_u_coeff_3;
					rho_v_coeff_33 = rho_v_coeff_3;
					e_coeff_33 = e_coeff_3;
                }
                
                // =====================================================================
                // r = 3 stencil 4 (boundary)
                // =====================================================================
                
                if (is_admissible_R34[c]) {
                    
                    index = 0; 
                    if ( N_face || S_face ) {                                        
                        d_rho_34.reinit(5); d_rho_u_34.reinit(5); d_rho_v_34.reinit(5); d_e_34.reinit(5);  
                    }
            
                    else {
                        d_rho_34.reinit(4); d_rho_u_34.reinit(4); d_rho_v_34.reinit(4); d_e_34.reinit(4);
                    }
                    
                    if (!N_face) {
                    
                        d_rho_34(index)    = RHO(neighbor_p2) - rho0;              // N neighbor 
                        d_rho_u_34(index)  = RHO_U(neighbor_p2) - rho_u0;          // N neighbor 
                        d_rho_v_34(index)  = RHO_V(neighbor_p2) - rho_v0;          // N neighbor 
                        d_e_34(index)      = E(neighbor_p2) - e0;                  // N neighbor 
                        
                        index++;
                    }
                    
                    if (!E_face) {
                    
                        d_rho_34(index)    = RHO(neighbor_p3) - rho0;              // E neighbor 
                        d_rho_u_34(index)  = RHO_U(neighbor_p3) - rho_u0;          // E neighbor 
                        d_rho_v_34(index)  = RHO_V(neighbor_p3) - rho_v0;          // E neighbor 
                        d_e_34(index)      = E(neighbor_p3) - e0;                  // E neighbor 
                        
                        index++;
                    }
                    
                    if (!S_face) {
                    
                        d_rho_34(index)    = RHO(neighbor_p4) - rho0;              // S neighbor 
                        d_rho_u_34(index)  = RHO_U(neighbor_p4) - rho_u0;          // S neighbor 
                        d_rho_v_34(index)  = RHO_V(neighbor_p4) - rho_v0;          // S neighbor 
                        d_e_34(index)      = E(neighbor_p4) - e0;                  // S neighbor 
                        
                        index++;
                    }
                    
                    if (EE) {

						local_neighbor_dof_indices[0] = EE_index;
                    
                        d_rho_34(index)    = RHO(local_neighbor_dof_indices[0] ) - rho0;              // S6 neighbor 
                        d_rho_u_34(index)  = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;          // S6 neighbor 
                        d_rho_v_34(index)  = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;          // S6 neighbor 
                        d_e_34(index)      = E(local_neighbor_dof_indices[0] ) - e0;                  // S6 neighbor 
                        
                        index++;
                    }
                    
                    if ( N_face || S_face ) {                                        
                        
                        // Transmissive boundary conditions 
                        
                        d_rho_34(index)    = 0.0;   d_rho_34(index+1)    = 0.0;            
                        d_rho_u_34(index)  = 0.0;   d_rho_u_34(index+1)  = 0.0;      
                        d_rho_v_34(index)  = 0.0;   d_rho_v_34(index+1)  = 0.0;    
                        d_e_34(index)      = 0.0;   d_e_34(index+1)  = 0.0;
                    }
                    
	                ROWS = cell_neighbor_index[c][1].size(); 

					b_rho_34.reinit(ROWS); b_rho_u_34.reinit(ROWS); b_rho_v_34.reinit(ROWS); b_e_34.reinit(ROWS);

					index = 0; 		

					// vertex neighbor of cell at face 3	
			
					for (unsigned int d = 0; d < cell_neighbor_index[c][1].size(); ++d) {						

						local_neighbor_dof_indices[0] = cell_neighbor_index[c][1][d];

						b_rho_34(index)   = RHO(local_neighbor_dof_indices[0] ) - rho0; 
						b_rho_u_34(index) = RHO_U(local_neighbor_dof_indices[0] ) - rho_u0;
						b_rho_v_34(index) = RHO_V(local_neighbor_dof_indices[0] ) - rho_v0;
						b_e_34(index)     = E(local_neighbor_dof_indices[0] ) - e0;
						index++; 
					}                   
					
                    CLS_R34[c].solve(b_rho_34, d_rho_34, rho_coeff_34);
                    CLS_R34[c].solve(b_rho_u_34, d_rho_u_34, rho_u_coeff_34);
                    CLS_R34[c].solve(b_rho_v_34, d_rho_v_34, rho_v_coeff_34);
                    CLS_R34[c].solve(b_e_34, d_e_34, e_coeff_34); 
                }
                
                else {
					
					rho_coeff_34 = rho_coeff_3; 
					rho_u_coeff_34 = rho_u_coeff_3;
					rho_v_coeff_34 = rho_v_coeff_3;
					e_coeff_34 = e_coeff_3;
                }
              
                // =====================================================================
                // r = 2 stencil 1 (boundary)
                // =====================================================================
                
                if (W_face) { // WEST boundary 
                
                    b_rho2(0) = (RHO(neighbor_p2) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p2) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p2) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p2) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 
                }
                
                else if (N_face) { // NORTH boundary 
                
                    b_rho2(0) = (RHO(neighbor_p1) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p1) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p1) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p1) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 
                }
                
                else {
                    
                    b_rho2(0)  = (RHO(neighbor_p1) - rho0);   // P1 neighbor 
                    b_rho2(1)  = (RHO(neighbor_p2) - rho0);   // P2 neighbor

                    b_rho_u_2(0)  = (RHO_U(neighbor_p1) - rho_u0);   // P1 neighbor 
                    b_rho_u_2(1)  = (RHO_U(neighbor_p2) - rho_u0);   // P2 neighbor

                    b_rho_v_2(0)  = (RHO_V(neighbor_p1) - rho_v0);   // P1 neighbor 
                    b_rho_v_2(1)  = (RHO_V(neighbor_p2) - rho_v0);   // P2 neighbor

                    b_e_2(0)  = (E(neighbor_p1) - e0);   // P1 neighbor 
                    b_e_2(1)  = (E(neighbor_p2) - e0);   // P2 neighbor

                }
                
				LU_R21[c].solve(b_rho2, rho_coeff_21);
                LU_R21[c].solve(b_rho_u_2, rho_u_coeff_21);
                LU_R21[c].solve(b_rho_v_2, rho_v_coeff_21);
                LU_R21[c].solve(b_e_2, e_coeff_21);

                // =====================================================================
                // r = 2 stencil 2 (boundary)
                // =====================================================================
                
                if (E_face) { // EAST boundary 
                
                    b_rho2(0) = (RHO(neighbor_p2) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p2) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p2) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p2) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 
                }
                
                else if (N_face) { // NORTH boundary 
                
                    b_rho2(0) = (RHO(neighbor_p3) - rho0);   // P3 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p3) - rho_u0);   // P3 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p3) - rho_v0);   // P3 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p3) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 

                }
                
                else {
                    
                    b_rho2(0) = (RHO(neighbor_p2) - rho0);   // P2 neighbor
                    b_rho2(1) = (RHO(neighbor_p3) - rho0);   // P3 neighbor
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p2) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = (RHO_U(neighbor_p3) - rho_u0);   // P3 neighbor

                    b_rho_v_2(0) = (RHO_V(neighbor_p2) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = (RHO_V(neighbor_p3) - rho_v0);   // P3 neighbor

                    b_e_2(0) = (E(neighbor_p2) - e0);   // P2 neighbor
                    b_e_2(1) = (E(neighbor_p3) - e0);   // P3 neighbor

                }
                
                LU_R22[c].solve(b_rho2, rho_coeff_22);
                LU_R22[c].solve(b_rho_u_2, rho_u_coeff_22);
                LU_R22[c].solve(b_rho_v_2, rho_v_coeff_22);
                LU_R22[c].solve(b_e_2, e_coeff_22);
                

                // =====================================================================
                // r = 2 stencil 3 (boundary)
                // =====================================================================
                
                if (E_face) { // EAST boundary 
                
                    b_rho2(0) = (RHO(neighbor_p4) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p4) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p4) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p4) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 

                }
                
                else if (S_face) { // SOUTH boundary 
                
                    b_rho2(0) = (RHO(neighbor_p3) - rho0);   // P3 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p3) - rho_u0);   // P3 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p3) - rho_v0);   // P3 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p3) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 

                }
                
                else {
                
                    b_rho2(0) = (RHO(neighbor_p3) - rho0);   // P3 neighbor
                    b_rho2(1) = (RHO(neighbor_p4) - rho0);   // P4 neighbor                                            // P4 neighbor
                     
                    b_rho_u_2(0) = (RHO_U(neighbor_p3) - rho_u0);   // P3 neighbor
                    b_rho_u_2(1) = (RHO_U(neighbor_p4) - rho_u0);   // P4 neighbor                                            // P4 neighbor
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p3) - rho_v0);   // P3 neighbor
                    b_rho_v_2(1) = (RHO_V(neighbor_p4) - rho_v0);   // P4 neighbor                                            // P4 neighbor
                    
                    b_e_2(0) = (E(neighbor_p3) - e0);   // P3 neighbor
                    b_e_2(1) = (E(neighbor_p4) - e0);   // P4 neighbor                                            // P4 neighbor

                }
                
                LU_R23[c].solve(b_rho2, rho_coeff_23);
                LU_R23[c].solve(b_rho_u_2, rho_u_coeff_23);
                LU_R23[c].solve(b_rho_v_2, rho_v_coeff_23);
                LU_R23[c].solve(b_e_2, e_coeff_23);
                
                // =====================================================================
                // r = 2 stencil 4 (boundary)
                // =====================================================================
                
                if (W_face) { // WEST boundary 
                
                    b_rho2(0) = (RHO(neighbor_p4) - rho0);   // P4 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p4) - rho_u0);   // P4 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p4) - rho_v0);   // P4 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p4) - e0);   // P4 neighbor
                    b_e_2(1) = 0.0; 

                    
                }
                
                else if (S_face) { // SOUTH boundary 
                
                    b_rho2(0) = (RHO(neighbor_p1) - rho0);   // P1 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p1) - rho_u0);   // P1 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p1) - rho_v0);   // P1 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p1) - e0);   // P1 1eighbor
                    b_e_2(1) = 0.0; 


                }
                
                else {
                    b_rho2(0) = (RHO(neighbor_p4) - rho0);   // P4 neighbor
                    b_rho2(1) = (RHO(neighbor_p1) - rho0);   // P1 neighbor

                    b_rho_u_2(0) = (RHO_U(neighbor_p4) - rho_u0);   // P4 neighbor
                    b_rho_u_2(1) = (RHO_U(neighbor_p1) - rho_u0);   // P1 neighbor
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p4) - rho_v0);   // P4 neighbor
                    b_rho_v_2(1) = (RHO_V(neighbor_p1) - rho_v0);   // P1 neighbor
                    
                    b_e_2(0) = (E(neighbor_p4) - e0);   // P4 neighbor
                    b_e_2(1) = (E(neighbor_p1) - e0);   // P1 neighbor
                    
                }
                
                LU_R24[c].solve(b_rho2, rho_coeff_24);
                LU_R24[c].solve(b_rho_u_2, rho_u_coeff_24);
                LU_R24[c].solve(b_rho_v_2, rho_v_coeff_24);
                LU_R24[c].solve(b_e_2, e_coeff_24);

				/*
                
                // Overwrite the x and y momentum components on the slip boundaries 
                
                for (unsigned int f=0; f < GeometryInfo<2>::faces_per_cell; ++f) {
                
                    if (cell->face(f)->at_boundary()) {
                    
                        if (cell->face(f)->boundary_id() == 2) { // boundary_id 2 corresponds to the reflecting wall 
                            
                            fv_face_values.reinit(cell, f);
                            face_normal_vector = fv_face_values.normal_vector(0);
                            nx = face_normal_vector[0]; ny = face_normal_vector[1];
                        
                            // =====================================================================
                            // r = 4 stencil (wall boundaries)
                            // =====================================================================
                            
                            Vector <double> d_4(7); Vector<double> uv_4(18);
                            
                            // Zero wall normal velocity
                            d_4(0) = -nx*rho_u0 - ny*rho_v0;  
                            d_4(1) = -nx*rho_u0 - ny*rho_v0; 
                            d_4(2) = -nx*rho_u0 - ny*rho_v0;
                            // Zero derivative tangential velocity
                            d_4(3) = 0.0; d_4(4) = 0.0; d_4(5) = 0.0;
                            // Zero third derivative tangential velocity
                            d_4(6) = 0.0; 
                            
                            ROWS = 0; index = 0; 

                            if (neighbor_p1 != -1) {
                                ROWS++; // P1 cell
                            } 

                            if (neighbor_p2 != -1) {
                                ROWS++; // P2 cell 
                            }

                            if (neighbor_p3 != -1) {
                                ROWS++; // P3 cell
                            }

                            if (neighbor_p4 != -1) {
                                ROWS++; // P4 cell 
                            }

                            if (neighbor_p1 != -1) {
                                if (cell->neighbor(0)->neighbor_index(0) != -1) {
                                    ROWS++;   // S1 cell 
                                }
                            }

                            if (neighbor_p1 != -1) {
                                if (cell->neighbor(0)->neighbor_index(3) != -1) {
                                    ROWS++; // S2 cell 
                                }
                            }

                            if (neighbor_p2 != -1) {
                                if (cell->neighbor(3)->neighbor_index(3) != -1) {
                                    ROWS++; // S3 cell
                                }
                            }

                            if (neighbor_p2 != -1) {
                                if (cell->neighbor(3)->neighbor_index(1) != -1) {
                                    ROWS++; // S4 cell 
                                }
                            }

                            if (neighbor_p3 != -1) {
                                if (cell->neighbor(1)->neighbor_index(1) != -1) {
                                    ROWS++; // S5 cell 
                                }
                            }

                            if (neighbor_p3 != -1) {
                                if (cell->neighbor(1)->neighbor_index(2) != -1) {
                                    ROWS++; // S6 cell 
                                }
                            }

                            if (neighbor_p4 != -1) {
                                if (cell->neighbor(2)->neighbor_index(2) != -1) {
                                    ROWS++; // S7 cell
                                }
                            }

                            if (neighbor_p4 != -1) {
                                if (cell->neighbor(2)->neighbor_index(0) != -1) {
                                    ROWS++; // S8 cell
                                }
                            }
                            
                            Vector<double> b_4(2*ROWS); 
                            
                            // Fill the vector 
                            
                            if (neighbor_p1 != -1) {
                            
                                // (P1 cell)
                                
                                b_4(index)      =  (RHO_U(neighbor_p1) - rho_u0); // P1 (u)
                                b_4(index+ROWS) =  (RHO_V(neighbor_p1) - rho_v0); // P1 (u)
                                
                                index++;

                            } 

                            if (neighbor_p2 != -1) {
                            
                                // (P2 cell)
                                
                                b_4(index)      =  (RHO_U(neighbor_p2) - rho_u0); // P1 (u)
                                b_4(index+ROWS) =  (RHO_V(neighbor_p2) - rho_v0); // P1 (u)
                                
                                index++;

                            }

                            if (neighbor_p3 != -1) {
                                // (P3 cell)

                                b_4(index)      =  (RHO_U(neighbor_p3) - rho_u0); // P1 (u)
                                b_4(index+ROWS) =  (RHO_V(neighbor_p3) - rho_v0); // P1 (u)
                                
                                index++;
                            }

                            if (neighbor_p4 != -1) {
                                // (P4 cell)
                                
                                b_4(index)      =  (RHO_U(neighbor_p4) - rho_u0); // P1 (u)
                                b_4(index+ROWS) =  (RHO_V(neighbor_p4) - rho_v0); // P1 (u)
                                
                                index++; 
                            }


                            // S1 cell 

                            if (neighbor_p1 != -1) {
                                if (cell->neighbor(0)->neighbor_index(0) != -1) {

                                    b_4(index)      =  (RHO_U(cell->neighbor(0)->neighbor_index(0)) - rho_u0); // P1 (u)
                                    b_4(index+ROWS) =  (RHO_V(cell->neighbor(0)->neighbor_index(0)) - rho_v0); // P1 (u)
                                    
                                    index++;
                                }
                            }

                            // S2 cell 

                            if (neighbor_p1 != -1) {
                                if (cell->neighbor(0)->neighbor_index(3) != -1) {
                                    
                                    b_4(index)      =  (RHO_U(cell->neighbor(0)->neighbor_index(3)) - rho_u0); // P1 (u)
                                    b_4(index+ROWS) =  (RHO_V(cell->neighbor(0)->neighbor_index(3)) - rho_v0); // P1 (u)
                                    index++; 
                                }
                            }

                            // S3 cell

                            if (neighbor_p2 != -1) {
                                if (cell->neighbor(3)->neighbor_index(3) != -1) {

                                    b_4(index)      =  (RHO_U(cell->neighbor(3)->neighbor_index(3)) - rho_u0); // P1 (u)
                                    b_4(index+ROWS) =  (RHO_V(cell->neighbor(3)->neighbor_index(3)) - rho_v0); // P1 (u)

                                    index++; 
                                }
                            }

                            // S4 cell 

                            if (neighbor_p2 != -1) {
                                if (cell->neighbor(3)->neighbor_index(1) != -1) {

                                    b_4(index)      =  (RHO_U(cell->neighbor(3)->neighbor_index(1)) - rho_u0); // P1 (u)
                                    b_4(index+ROWS) =  (RHO_V(cell->neighbor(3)->neighbor_index(1)) - rho_v0); // P1 (u)

                                    index++; 
                                }
                            }

                            // S5 cell

                            if (neighbor_p3 != -1) {
                                if (cell->neighbor(1)->neighbor_index(1) != -1) {

                                    b_4(index)      =  (RHO_U(cell->neighbor(1)->neighbor_index(1)) - rho_u0); // P1 (u)
                                    b_4(index+ROWS) =  (RHO_V(cell->neighbor(1)->neighbor_index(1)) - rho_v0); // P1 (u)

                                    index++; 
                                }
                            }

                            // S6 cell

                            if (neighbor_p3 != -1) {
                                if (cell->neighbor(1)->neighbor_index(2) != -1) {

                                    b_4(index)      =  (RHO_U(cell->neighbor(1)->neighbor_index(2)) - rho_u0); // P1 (u)
                                    b_4(index+ROWS) =  (RHO_V(cell->neighbor(1)->neighbor_index(2)) - rho_v0); // P1 (u)

                                    index++; 
                                }
                            }

                            // S7 cell

                            if (neighbor_p4 != -1) {
                                if (cell->neighbor(2)->neighbor_index(2) != -1) {

                                    b_4(index)      =  (RHO_U(cell->neighbor(2)->neighbor_index(2)) - rho_u0); // P1 (u)
                                    b_4(index+ROWS) =  (RHO_V(cell->neighbor(2)->neighbor_index(2)) - rho_v0); // P1 (u)

                                    index++; 
                                }
                            }

                            // S8 cell

                            if (neighbor_p4 != -1) {
                                if (cell->neighbor(2)->neighbor_index(0) != -1) {

                                    b_4(index)      =  (RHO_U(cell->neighbor(2)->neighbor_index(0)) - rho_u0); // P1 (u)
                                    b_4(index+ROWS) =  (RHO_V(cell->neighbor(2)->neighbor_index(0)) - rho_v0); // P1 (u)

                                    index++; 
                                }
                    
                            }
                            
                            CLS_R4_slip[wall_boundary_global_index_map[c]].solve(b_4, d_4, uv_4);
                            
                            for (unsigned int i = 0; i < 9; i++) {
                                rho_u_coeff_4[i] = uv_4[i]; 
                                rho_v_coeff_4[i] = uv_4[i+9];
                            }
                            
                            // =====================================================================
                            // r = 3 
                            // =====================================================================
                            
                            Vector <double> b_3(6); Vector <double> d_3(5); Vector<double> uv_3(10); 
                            
                            index = 0; 
                            
                            
                            if (neighbor_p1 != -1 ) {
                                b_3(index)   = (RHO_U(neighbor_p1) - rho_u0); // P1 (rho_u)
                                b_3(index+3) = (RHO_V(neighbor_p1) - rho_v0); // P1 (rho_v)
                                index++; 
                            }
                            
                            if (neighbor_p2 != -1 ) {
                                b_3(index)   = (RHO_U(neighbor_p2) - rho_u0);  // P2 (u)
                                b_3(index+3) = (RHO_V(neighbor_p2) - rho_v0);  // P2 (v)
                                index++; 
                            }
                            
                            if (neighbor_p3 != -1 ) {
                                b_3(index)   = (RHO_U(neighbor_p3) - rho_u0);  // P3 (u)
                                b_3(index+3) = (RHO_V(neighbor_p3) - rho_v0);  // P3 (v)
                                index++; 
                            }
                            
                            if (neighbor_p4 != -1 ) {
                                b_3(index)   = (RHO_U(neighbor_p4) - rho_u0);  // P4 (u)
                                b_3(index+3) = (RHO_V(neighbor_p4) - rho_v0);  // P4 (v)
                                index++; 
                            }
                            
                            fv_face_values.reinit(cell, f);
                            face_normal_vector = fv_face_values.normal_vector(0);
                            nx = face_normal_vector[0]; ny = face_normal_vector[1];

                            //face_center = cell->face(f)->center(); 
                
                            // Zero wall normal velocity
                            d_3(0) = -nx*rho_u0 - ny*rho_v0;  
                            d_3(1) = -nx*rho_u0 - ny*rho_v0; 
                            d_3(2) = -nx*rho_u0 - ny*rho_v0;
                            // Zero derivative tangential velocity
                            d_3(3) = 0.0; d_3(4) = 0.0; 
                            
                            CLS_R3_slip[wall_boundary_global_index_map[c]].solve(b_3, d_3, uv_3); 
                            
                            for (unsigned int i = 0; i < 5; i++) {
                                rho_u_coeff_3[i] = uv_3[i]; 
                                rho_v_coeff_3[i] = uv_3[i+5]; 
                            } 
                            
                            // Set other coefficients to centered third order coefficients 
                            
                            rho_u_coeff_31 = rho_u_coeff_3; rho_v_coeff_31 = rho_v_coeff_3;
                            rho_u_coeff_31 = rho_u_coeff_3; rho_v_coeff_31 = rho_v_coeff_3;
                            rho_u_coeff_33 = rho_u_coeff_3; rho_v_coeff_33 = rho_v_coeff_3;
                            rho_u_coeff_34 = rho_u_coeff_3; rho_v_coeff_34 = rho_v_coeff_3;
                            rho_u_coeff_35 = rho_u_coeff_3; rho_v_coeff_35 = rho_v_coeff_3;
                            rho_u_coeff_36 = rho_u_coeff_3; rho_v_coeff_36 = rho_v_coeff_3;
                            rho_u_coeff_34 = rho_u_coeff_3; rho_v_coeff_34 = rho_v_coeff_3;
                            rho_u_coeff_38 = rho_u_coeff_3; rho_v_coeff_38 = rho_v_coeff_3;
                        }
                    }
                } // End of wall boundary loop  */ 
//	   			pcout<<"boundary end"<<std::endl;                    
            } // Non-corner boundary cell loop 
            
            else {
                
                // Corner boundary cells - reduce to first order 
                
                rho_coeff_4 = 0.0;
				rho_coeff_3 = 0.0;
				rho_coeff_31 = 0.0;
				rho_coeff_32 = 0.0;
                rho_coeff_33 = 0.0;  
                rho_coeff_34 = 0.0;  
				rho_coeff_21 = 0.0;
				rho_coeff_22 = 0.0;
				rho_coeff_23 = 0.0;
				rho_coeff_24 = 0.0;
				
				rho_u_coeff_4 = 0.0;
				rho_u_coeff_3 = 0.0;
				rho_u_coeff_31 = 0.0;
				rho_u_coeff_32 = 0.0;
                rho_u_coeff_32 = 0.0;  
                rho_u_coeff_34 = 0.0;  
				rho_u_coeff_21 = 0.0;
				rho_u_coeff_22 = 0.0;
				rho_u_coeff_23 = 0.0;
				rho_u_coeff_24 = 0.0;
				
				rho_v_coeff_4 = 0.0;
				rho_v_coeff_3 = 0.0;
				rho_v_coeff_31 = 0.0;
				rho_v_coeff_32 = 0.0;
                rho_v_coeff_33 = 0.0;  
                rho_v_coeff_34 = 0.0;  
				rho_v_coeff_21 = 0.0;
				rho_v_coeff_22 = 0.0;
				rho_v_coeff_23 = 0.0;
				rho_v_coeff_24 = 0.0;
				
				e_coeff_4 = 0.0;
				e_coeff_3 = 0.0;
				e_coeff_31 = 0.0;
				e_coeff_32 = 0.0;
                e_coeff_33 = 0.0;  
                e_coeff_34 = 0.0;  
				e_coeff_21 = 0.0;
				e_coeff_22 = 0.0;
				e_coeff_23 = 0.0;
				e_coeff_24 = 0.0;
                
            } // Corner boundary cell loop 
            
        } // End of boundary cell loop 
        
        // Find the smoothness indicators 
        
		// =====================================================================
		// r = 4 Stencil 
		// =====================================================================
        
		IS_RHO(0)   = compute_fourth_order_smoothness_indicator(rho_coeff_4, IS_constants[c], h);
		IS_RHO_U(0) = compute_fourth_order_smoothness_indicator(rho_u_coeff_4, IS_constants[c], h);
		IS_RHO_V(0) = compute_fourth_order_smoothness_indicator(rho_v_coeff_4, IS_constants[c], h);
		IS_E(0)     = compute_fourth_order_smoothness_indicator(e_coeff_4, IS_constants[c], h);

		// =====================================================================
		// r = 3
		// =====================================================================
        
		IS_RHO(1)   = compute_third_order_smoothness_indicator(rho_coeff_3, IS_constants[c], h); 
		IS_RHO_U(1) = compute_third_order_smoothness_indicator(rho_u_coeff_3, IS_constants[c], h);
		IS_RHO_V(1) = compute_third_order_smoothness_indicator(rho_v_coeff_3, IS_constants[c], h);
		IS_E(1)     = compute_third_order_smoothness_indicator(e_coeff_3, IS_constants[c], h);
        
		// =====================================================================
		// r = 3  stencil 1
		// =====================================================================
        
		IS_RHO(2)   = compute_third_order_smoothness_indicator(rho_coeff_31, IS_constants[c], h); 
		IS_RHO_U(2) = compute_third_order_smoothness_indicator(rho_u_coeff_31, IS_constants[c], h);
		IS_RHO_V(2) = compute_third_order_smoothness_indicator(rho_v_coeff_31, IS_constants[c], h);
		IS_E(2)     = compute_third_order_smoothness_indicator(e_coeff_31, IS_constants[c], h);
        
        
		// =====================================================================
		// r = 3  stencil 2
		// =====================================================================
        
		IS_RHO(3)   = compute_third_order_smoothness_indicator(rho_coeff_32, IS_constants[c], h); 
		IS_RHO_U(3) = compute_third_order_smoothness_indicator(rho_u_coeff_32, IS_constants[c], h);
		IS_RHO_V(3) = compute_third_order_smoothness_indicator(rho_v_coeff_32, IS_constants[c], h);
		IS_E(3)     = compute_third_order_smoothness_indicator(e_coeff_32, IS_constants[c], h);
        
		// =====================================================================
		// r = 3  stencil 3
		// =====================================================================
        
		IS_RHO(4)   = compute_third_order_smoothness_indicator(rho_coeff_33, IS_constants[c], h); 
		IS_RHO_U(4) = compute_third_order_smoothness_indicator(rho_u_coeff_33, IS_constants[c], h);
		IS_RHO_V(4) = compute_third_order_smoothness_indicator(rho_v_coeff_33, IS_constants[c], h);
		IS_E(4)     = compute_third_order_smoothness_indicator(e_coeff_33, IS_constants[c], h);
        
		// =====================================================================
		// r = 3  stencil 4
		// =====================================================================
        
		IS_RHO(5)   = compute_third_order_smoothness_indicator(rho_coeff_34, IS_constants[c], h); 
		IS_RHO_U(5) = compute_third_order_smoothness_indicator(rho_u_coeff_34, IS_constants[c], h);
		IS_RHO_V(5) = compute_third_order_smoothness_indicator(rho_v_coeff_34, IS_constants[c], h);
		IS_E(5)     = compute_third_order_smoothness_indicator(e_coeff_34, IS_constants[c], h);
        
		// =====================================================================
		// r = 2  stencil 1
		// =====================================================================
		
		IS_RHO(6)   = compute_second_order_smoothness_indicator(rho_coeff_21, IS_constants[c], h); 
		IS_RHO_U(6) = compute_second_order_smoothness_indicator(rho_u_coeff_21, IS_constants[c], h);
		IS_RHO_V(6) = compute_second_order_smoothness_indicator(rho_v_coeff_21, IS_constants[c], h);
		IS_E(6)     = compute_second_order_smoothness_indicator(e_coeff_21, IS_constants[c], h);
		
		// =====================================================================
		// r = 2  stencil 2
		// =====================================================================
		
		IS_RHO(7)   = compute_second_order_smoothness_indicator(rho_coeff_22, IS_constants[c], h); 
		IS_RHO_U(7) = compute_second_order_smoothness_indicator(rho_u_coeff_22, IS_constants[c], h);
		IS_RHO_V(7) = compute_second_order_smoothness_indicator(rho_v_coeff_22, IS_constants[c], h);
		IS_E(7)     = compute_second_order_smoothness_indicator(e_coeff_22, IS_constants[c], h);
		
		// =====================================================================
		// r = 2  stencil 3
		// =====================================================================
		
		IS_RHO(8)   = compute_second_order_smoothness_indicator(rho_coeff_23, IS_constants[c], h); 
		IS_RHO_U(8) = compute_second_order_smoothness_indicator(rho_u_coeff_23, IS_constants[c], h);
		IS_RHO_V(8) = compute_second_order_smoothness_indicator(rho_v_coeff_23, IS_constants[c], h);
		IS_E(8)     = compute_second_order_smoothness_indicator(e_coeff_23, IS_constants[c], h);
		
		// =====================================================================
		// r = 2  stencil 4
		// =====================================================================
		
		IS_RHO(9)   = compute_second_order_smoothness_indicator(rho_coeff_24, IS_constants[c], h); 
		IS_RHO_U(9) = compute_second_order_smoothness_indicator(rho_u_coeff_24, IS_constants[c], h);
		IS_RHO_V(9) = compute_second_order_smoothness_indicator(rho_v_coeff_24, IS_constants[c], h);
		IS_E(9)     = compute_second_order_smoothness_indicator(e_coeff_24, IS_constants[c], h);
		
		// Combine the polynomials 
		
		sum_RHO = 0.0; sum_RHO_U = 0.0; sum_RHO_V = 0.0; sum_E = 0.0; sum_gamma = 0.0; 
		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			
			w_RHO(j)   = gamma(j)/(std::pow((IS_RHO(j) + epsilon), p));
			w_RHO_U(j) = gamma(j)/(std::pow((IS_RHO_U(j) + epsilon), p));
			w_RHO_V(j) = gamma(j)/(std::pow((IS_RHO_V(j) + epsilon), p));
			w_E(j)     = gamma(j)/(std::pow((IS_E(j) + epsilon), p));
			
			sum_RHO += w_RHO(j); sum_RHO_U += w_RHO_U(j); sum_RHO_V += w_RHO_V(j); sum_E += w_E(j);
			
			sum_gamma += gamma(j); 
		}
		
		// Normalize the weights 
		
		for (unsigned int j = 0; j < no_stencils; j++ ) {
			w_RHO(j) = w_RHO(j)/sum_RHO; 
			w_RHO_U(j) = w_RHO_U(j)/sum_RHO_U;
			w_RHO_V(j) = w_RHO_V(j)/sum_RHO_V;
			w_E(j) = w_E(j)/sum_E;
			gamma(j) = gamma(j)/sum_gamma; 
		}
		
		// Density 
		
		coeffs_RHO[c](0) = rho0; 
		
		coeffs_RHO[c](1) = (w_RHO(0)/gamma(0))*(rho_coeff_4(0) - 
												
												gamma(1)*rho_coeff_3(0)  - 
												gamma(2)*rho_coeff_31(0) - 
												gamma(3)*rho_coeff_32(0) - 
												gamma(4)*rho_coeff_33(0) - 
												gamma(5)*rho_coeff_34(0) - 
												gamma(6)*rho_coeff_21(0) -
												gamma(7)*rho_coeff_22(0) -
												gamma(8)*rho_coeff_23(0) -
												gamma(9)*rho_coeff_24(0) ) + 
												
												w_RHO(1)*rho_coeff_3(0) +
												w_RHO(2)*rho_coeff_31(0) + 
												w_RHO(3)*rho_coeff_32(0) +
												w_RHO(4)*rho_coeff_33(0) +
												w_RHO(5)*rho_coeff_34(0) +
												w_RHO(6)*rho_coeff_21(0) +
												w_RHO(7)*rho_coeff_22(0) +
												w_RHO(8)*rho_coeff_23(0) +
												w_RHO(9)*rho_coeff_24(0) ;// u_x
		
		coeffs_RHO[c](2) = (w_RHO(0)/gamma(0))*(rho_coeff_4(1) - 
												
												gamma(1)*rho_coeff_3(1) - 
												gamma(2)*rho_coeff_31(1) - 
												gamma(3)*rho_coeff_32(1) -
												gamma(4)*rho_coeff_33(1) -
												gamma(5)*rho_coeff_34(1) - 
												gamma(6)*rho_coeff_21(1) -
												gamma(7)*rho_coeff_22(1) -
												gamma(8)*rho_coeff_23(1) -
												gamma(9)*rho_coeff_24(1) ) + 
												
												w_RHO(1)*rho_coeff_3(1) +
												w_RHO(2)*rho_coeff_31(1) +
												w_RHO(3)*rho_coeff_32(1) +
												w_RHO(4)*rho_coeff_33(1) +
												w_RHO(5)*rho_coeff_34(1) +
												w_RHO(6)*rho_coeff_21(1) +
												w_RHO(7)*rho_coeff_22(1) +
												w_RHO(8)*rho_coeff_23(1) +
												w_RHO(9)*rho_coeff_24(1) ; // u_x

		coeffs_RHO[c](3) = (w_RHO(0)/gamma(0))*(rho_coeff_4(2) - 
												
												gamma(1)*rho_coeff_3(2) - 
												gamma(2)*rho_coeff_31(2) - 
												gamma(3)*rho_coeff_32(2) - 
												gamma(4)*rho_coeff_33(2) - 
												gamma(5)*rho_coeff_34(2) ) + 
												
												w_RHO(1)*rho_coeff_3(2) +
												w_RHO(2)*rho_coeff_31(2) +
												w_RHO(3)*rho_coeff_32(2) + 
												w_RHO(4)*rho_coeff_33(2) + 
												w_RHO(5)*rho_coeff_34(2) ; // u_xx
												
		coeffs_RHO[c](4) = (w_RHO(0)/gamma(0))*(rho_coeff_4(3) - 
												
												gamma(1)*rho_coeff_3(3) - 
												gamma(2)*rho_coeff_31(3) - 
												gamma(3)*rho_coeff_32(3) - 
												gamma(4)*rho_coeff_33(3) - 
												gamma(5)*rho_coeff_34(3) ) + 
												
												w_RHO(1)*rho_coeff_3(3) +
												w_RHO(2)*rho_coeff_31(3) + 
												w_RHO(3)*rho_coeff_32(3) + 
												w_RHO(4)*rho_coeff_33(3) + 
												w_RHO(5)*rho_coeff_34(3) ; // u_yy
		
		coeffs_RHO[c](5) = (w_RHO(0)/gamma(0))*(rho_coeff_4(4) - 
												
												gamma(1)*rho_coeff_3(4) - 
												gamma(2)*rho_coeff_31(4) - 
												gamma(3)*rho_coeff_32(4) - 
												gamma(4)*rho_coeff_33(4) - 
												gamma(5)*rho_coeff_34(4) ) + 
												
												w_RHO(1)*rho_coeff_3(4) +
												w_RHO(2)*rho_coeff_31(4) + 
												w_RHO(3)*rho_coeff_32(4) + 
												w_RHO(4)*rho_coeff_33(4) + 
												w_RHO(5)*rho_coeff_34(4) ; // u_xy
		
		coeffs_RHO[c](6) = (w_RHO(0)/gamma(0))*(rho_coeff_4(5)); // u_xxx

		coeffs_RHO[c](7) = (w_RHO(0)/gamma(0))*(rho_coeff_4(6)); // u_yyy

		coeffs_RHO[c](8) = (w_RHO(0)/gamma(0))*(rho_coeff_4(7)); // u_xxy

		coeffs_RHO[c](9) = (w_RHO(0)/gamma(0))*(rho_coeff_4(8)); // u_xyy
		
		// x-momentum
		
		coeffs_RHO_U[c](0) = rho_u0; 
		
		coeffs_RHO_U[c](1) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(0) - 
												
												gamma(1)*rho_u_coeff_3(0)  - 
												gamma(2)*rho_u_coeff_31(0) - 
												gamma(3)*rho_u_coeff_32(0) - 
												gamma(4)*rho_u_coeff_33(0) - 
												gamma(5)*rho_u_coeff_34(0) - 
												gamma(6)*rho_u_coeff_21(0) -
												gamma(7)*rho_u_coeff_22(0) -
												gamma(8)*rho_u_coeff_23(0) -
												gamma(9)*rho_u_coeff_24(0) ) + 
												
												w_RHO_U(1)*rho_u_coeff_3(0) +
												w_RHO_U(2)*rho_u_coeff_31(0) + 
												w_RHO_U(3)*rho_u_coeff_32(0) +
												w_RHO_U(4)*rho_u_coeff_33(0) +
												w_RHO_U(5)*rho_u_coeff_34(0) +
												w_RHO_U(6)*rho_u_coeff_21(0) +
												w_RHO_U(7)*rho_u_coeff_22(0) +
												w_RHO_U(8)*rho_u_coeff_23(0) +
												w_RHO_U(9)*rho_u_coeff_24(0) ;// u_x
		
		coeffs_RHO_U[c](2) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(1) - 
												
												gamma(1)*rho_u_coeff_3(1) - 
												gamma(2)*rho_u_coeff_31(1) - 
												gamma(3)*rho_u_coeff_32(1) -
												gamma(4)*rho_u_coeff_33(1) -
												gamma(5)*rho_u_coeff_34(1) - 
												gamma(6)*rho_u_coeff_21(1) -
												gamma(7)*rho_u_coeff_22(1) -
												gamma(8)*rho_u_coeff_23(1) -
												gamma(9)*rho_u_coeff_24(1) ) + 
												
												w_RHO_U(1)*rho_u_coeff_3(1) +
												w_RHO_U(2)*rho_u_coeff_31(1) +
												w_RHO_U(3)*rho_u_coeff_32(1) +
												w_RHO_U(4)*rho_u_coeff_33(1) +
												w_RHO_U(5)*rho_u_coeff_34(1) +
												w_RHO_U(6)*rho_u_coeff_21(1) +
												w_RHO_U(7)*rho_u_coeff_22(1) +
												w_RHO_U(8)*rho_u_coeff_23(1) +
												w_RHO_U(9)*rho_u_coeff_24(1) ; // u_x

		coeffs_RHO_U[c](3) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(2) - 
												
												gamma(1)*rho_u_coeff_3(2) - 
												gamma(2)*rho_u_coeff_31(2) - 
												gamma(3)*rho_u_coeff_32(2) - 
												gamma(4)*rho_u_coeff_33(2) - 
												gamma(5)*rho_u_coeff_34(2) ) + 
												
												w_RHO_U(1)*rho_u_coeff_3(2) +
												w_RHO_U(2)*rho_u_coeff_31(2) +
												w_RHO_U(3)*rho_u_coeff_32(2) + 
												w_RHO_U(4)*rho_u_coeff_33(2) + 
												w_RHO_U(5)*rho_u_coeff_34(2) ; // u_xx
												
		coeffs_RHO_U[c](4) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(3) - 
												
												gamma(1)*rho_u_coeff_3(3) - 
												gamma(2)*rho_u_coeff_31(3) - 
												gamma(3)*rho_u_coeff_32(3) - 
												gamma(4)*rho_u_coeff_33(3) - 
												gamma(5)*rho_u_coeff_34(3) ) + 
												
												w_RHO_U(1)*rho_u_coeff_3(3) +
												w_RHO_U(2)*rho_u_coeff_31(3) + 
												w_RHO_U(3)*rho_u_coeff_32(3) + 
												w_RHO_U(4)*rho_u_coeff_33(3) + 
												w_RHO_U(5)*rho_u_coeff_34(3) ; // u_yy
		
		coeffs_RHO_U[c](5) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(4) - 
												
												gamma(1)*rho_u_coeff_3(4) - 
												gamma(2)*rho_u_coeff_31(4) - 
												gamma(3)*rho_u_coeff_32(4) - 
												gamma(4)*rho_u_coeff_33(4) - 
												gamma(5)*rho_u_coeff_34(4) ) + 
												
												w_RHO_U(1)*rho_u_coeff_3(4) +
												w_RHO_U(2)*rho_u_coeff_31(4) + 
												w_RHO_U(3)*rho_u_coeff_32(4) + 
												w_RHO_U(4)*rho_u_coeff_33(4) + 
												w_RHO_U(5)*rho_u_coeff_34(4) ; // u_xy
		
		coeffs_RHO_U[c](6) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(5)); // u_xxx

		coeffs_RHO_U[c](7) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(6)); // u_yyy

		coeffs_RHO_U[c](8) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(7)); // u_xxy

		coeffs_RHO_U[c](9) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_4(8)); // u_xyy

		// y-momentum
		
		coeffs_RHO_V[c](0) = rho_v0; 
		
		coeffs_RHO_V[c](1) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(0) - 
												
												gamma(1)*rho_v_coeff_3(0)  - 
												gamma(2)*rho_v_coeff_31(0) - 
												gamma(3)*rho_v_coeff_32(0) - 
												gamma(4)*rho_v_coeff_33(0) - 
												gamma(5)*rho_v_coeff_34(0) - 
												gamma(6)*rho_v_coeff_21(0) -
												gamma(7)*rho_v_coeff_22(0) -
												gamma(8)*rho_v_coeff_23(0) -
												gamma(9)*rho_v_coeff_24(0) ) + 
												
												w_RHO_V(1)*rho_v_coeff_3(0) +
												w_RHO_V(2)*rho_v_coeff_31(0) + 
												w_RHO_V(3)*rho_v_coeff_32(0) +
												w_RHO_V(4)*rho_v_coeff_33(0) +
												w_RHO_V(5)*rho_v_coeff_34(0) +
												w_RHO_V(6)*rho_v_coeff_21(0) +
												w_RHO_V(7)*rho_v_coeff_22(0) +
												w_RHO_V(8)*rho_v_coeff_23(0) +
												w_RHO_V(9)*rho_v_coeff_24(0) ;// v_x
		
		coeffs_RHO_V[c](2) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(1) - 
												
												gamma(1)*rho_v_coeff_3(1) - 
												gamma(2)*rho_v_coeff_31(1) - 
												gamma(3)*rho_v_coeff_32(1) -
												gamma(4)*rho_v_coeff_33(1) -
												gamma(5)*rho_v_coeff_34(1) - 
												gamma(6)*rho_v_coeff_21(1) -
												gamma(7)*rho_v_coeff_22(1) -
												gamma(8)*rho_v_coeff_23(1) -
												gamma(9)*rho_v_coeff_24(1) ) + 
												
												w_RHO_V(1)*rho_v_coeff_3(1) +
												w_RHO_V(2)*rho_v_coeff_31(1) +
												w_RHO_V(3)*rho_v_coeff_32(1) +
												w_RHO_V(4)*rho_v_coeff_33(1) +
												w_RHO_V(5)*rho_v_coeff_34(1) +
												w_RHO_V(6)*rho_v_coeff_21(1) +
												w_RHO_V(7)*rho_v_coeff_22(1) +
												w_RHO_V(8)*rho_v_coeff_23(1) +
												w_RHO_V(9)*rho_v_coeff_24(1) ; // v_x

		coeffs_RHO_V[c](3) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(2) - 
												
												gamma(1)*rho_v_coeff_3(2) - 
												gamma(2)*rho_v_coeff_31(2) - 
												gamma(3)*rho_v_coeff_32(2) - 
												gamma(4)*rho_v_coeff_33(2) - 
												gamma(5)*rho_v_coeff_34(2) ) + 
												
												w_RHO_V(1)*rho_v_coeff_3(2) +
												w_RHO_V(2)*rho_v_coeff_31(2) +
												w_RHO_V(3)*rho_v_coeff_32(2) + 
												w_RHO_V(4)*rho_v_coeff_33(2) + 
												w_RHO_V(5)*rho_v_coeff_34(2) ; // v_xx
												
		coeffs_RHO_V[c](4) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(3) - 
												
												gamma(1)*rho_v_coeff_3(3) - 
												gamma(2)*rho_v_coeff_31(3) - 
												gamma(3)*rho_v_coeff_32(3) - 
												gamma(4)*rho_v_coeff_33(3) - 
												gamma(5)*rho_v_coeff_34(3) ) + 
												
												w_RHO_V(1)*rho_v_coeff_3(3) +
												w_RHO_V(2)*rho_v_coeff_31(3) + 
												w_RHO_V(3)*rho_v_coeff_32(3) + 
												w_RHO_V(4)*rho_v_coeff_33(3) + 
												w_RHO_V(5)*rho_v_coeff_34(3) ; // v_yy
		
		coeffs_RHO_V[c](5) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(4) - 
												
												gamma(1)*rho_v_coeff_3(4) - 
												gamma(2)*rho_v_coeff_31(4) - 
												gamma(3)*rho_v_coeff_32(4) - 
												gamma(4)*rho_v_coeff_33(4) - 
												gamma(5)*rho_v_coeff_34(4) ) + 
												
												w_RHO_V(1)*rho_v_coeff_3(4) +
												w_RHO_V(2)*rho_v_coeff_31(4) + 
												w_RHO_V(3)*rho_v_coeff_32(4) + 
												w_RHO_V(4)*rho_v_coeff_33(4) + 
												w_RHO_V(5)*rho_v_coeff_34(4) ; // v_xy
		
		coeffs_RHO_V[c](6) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(5)); // v_xxx

		coeffs_RHO_V[c](7) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(6)); // v_yyy

		coeffs_RHO_V[c](8) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(7)); // v_xxy

		coeffs_RHO_V[c](9) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_4(8)); // v_xyy
		
		// Total energy  
		
		coeffs_E[c](0) = e0; 
		
		coeffs_E[c](1) = (w_E(0)/gamma(0))*(e_coeff_4(0) - 
												
												gamma(1)*e_coeff_3(0)  - 
												gamma(2)*e_coeff_31(0) - 
												gamma(3)*e_coeff_32(0) - 
												gamma(4)*e_coeff_33(0) - 
												gamma(5)*e_coeff_34(0) - 
												gamma(6)*e_coeff_21(0) -
												gamma(7)*e_coeff_22(0) -
												gamma(8)*e_coeff_23(0) -
												gamma(9)*e_coeff_24(0) ) + 
												
												w_E(1)*e_coeff_3(0) +
												w_E(2)*e_coeff_31(0) + 
												w_E(3)*e_coeff_32(0) +
												w_E(4)*e_coeff_33(0) +
												w_E(5)*e_coeff_34(0) +
												w_E(6)*e_coeff_21(0) +
												w_E(7)*e_coeff_22(0) +
												w_E(8)*e_coeff_23(0) +
												w_E(9)*e_coeff_24(0) ;// v_x
		
		coeffs_E[c](2) = (w_E(0)/gamma(0))*(e_coeff_4(1) - 
												
												gamma(1)*e_coeff_3(1) - 
												gamma(2)*e_coeff_31(1) - 
												gamma(3)*e_coeff_32(1) -
												gamma(4)*e_coeff_33(1) -
												gamma(5)*e_coeff_34(1) - 
												gamma(6)*e_coeff_21(1) -
												gamma(7)*e_coeff_22(1) -
												gamma(8)*e_coeff_23(1) -
												gamma(9)*e_coeff_24(1) ) + 
												
												w_E(1)*e_coeff_3(1) +
												w_E(2)*e_coeff_31(1) +
												w_E(3)*e_coeff_32(1) +
												w_E(4)*e_coeff_33(1) +
												w_E(5)*e_coeff_34(1) +
												w_E(6)*e_coeff_21(1) +
												w_E(7)*e_coeff_22(1) +
												w_E(8)*e_coeff_23(1) +
												w_E(9)*e_coeff_24(1) ; // v_x

		coeffs_E[c](3) = (w_E(0)/gamma(0))*(e_coeff_4(2) - 
												
												gamma(1)*e_coeff_3(2) - 
												gamma(2)*e_coeff_31(2) - 
												gamma(3)*e_coeff_32(2) - 
												gamma(4)*e_coeff_33(2) - 
												gamma(5)*e_coeff_34(2) ) + 
												
												w_E(1)*e_coeff_3(2) +
												w_E(2)*e_coeff_31(2) +
												w_E(3)*e_coeff_32(2) + 
												w_E(4)*e_coeff_33(2) + 
												w_E(5)*e_coeff_34(2) ; // v_xx
												
		coeffs_E[c](4) = (w_E(0)/gamma(0))*(e_coeff_4(3) - 
												
												gamma(1)*e_coeff_3(3) - 
												gamma(2)*e_coeff_31(3) - 
												gamma(3)*e_coeff_32(3) - 
												gamma(4)*e_coeff_33(3) - 
												gamma(5)*e_coeff_34(3) ) + 
												
												w_E(1)*e_coeff_3(3) +
												w_E(2)*e_coeff_31(3) + 
												w_E(3)*e_coeff_32(3) + 
												w_E(4)*e_coeff_33(3) + 
												w_E(5)*e_coeff_34(3) ; // v_yy
		
		coeffs_E[c](5) = (w_E(0)/gamma(0))*(e_coeff_4(4) - 
												
												gamma(1)*e_coeff_3(4) - 
												gamma(2)*e_coeff_31(4) - 
												gamma(3)*e_coeff_32(4) - 
												gamma(4)*e_coeff_33(4) - 
												gamma(5)*e_coeff_34(4) ) + 
												
												w_E(1)*e_coeff_3(4) +
												w_E(2)*e_coeff_31(4) + 
												w_E(3)*e_coeff_32(4) + 
												w_E(4)*e_coeff_33(4) + 
												w_E(5)*e_coeff_34(4) ; // v_xy
		
		coeffs_E[c](6) = (w_E(0)/gamma(0))*(e_coeff_4(5)); // u_xxx

		coeffs_E[c](7) = (w_E(0)/gamma(0))*(e_coeff_4(6)); // u_yyy

		coeffs_E[c](8) = (w_E(0)/gamma(0))*(e_coeff_4(7)); // u_xxy

		coeffs_E[c](9) = (w_E(0)/gamma(0))*(e_coeff_4(8)); // u_xyy

		// check bounds 
        if(!Use_ader)
        for (unsigned int f = 0; f < faces_per_cell; ++f) {

            face_quadrature_point_1 = Cell[c].face_quadrature_point1(f);
            face_quadrature_point_2 = Cell[c].face_quadrature_point2(f);

            U1(0) = evaluate_weno_polynomial(coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point_1, h);
            U1(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point_1, h);
            U1(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point_1, h);
            U1(3) = evaluate_weno_polynomial(coeffs_E[c], WENO_poly_consts[c], face_quadrature_point_1, h);

            U2(0) = evaluate_weno_polynomial(coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point_2, h);
            U2(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point_2, h);
            U2(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point_2, h);
            U2(3) = evaluate_weno_polynomial(coeffs_E[c], WENO_poly_consts[c], face_quadrature_point_2, h);

            W1 = conserved_to_primitive(U1);
            W2 = conserved_to_primitive(U2);

            if (W1(0) < 0.0 || W1(3) < 0.0 || W2(0) < 0.0 || W2(3) < 0.0) {

                coeffs_RHO[c](0) = rho0;  coeffs_RHO_U[c](0) = rho_u0;    coeffs_RHO_V[c](0) = rho_v0;    coeffs_E[c](0) = e0; 

                for (unsigned int i = 1; i < 10; i++) {
                  coeffs_RHO[c](i) = 0.0; 
                  coeffs_RHO_U[c](i) = 0.0;
                  coeffs_RHO_V[c](i) = 0.0;
                  coeffs_E[c](i) = 0.0;
                }

            }
        }
        
    } // End of cell loop 
    
} // End of function 

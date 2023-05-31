#include "../include/Weno32.h"

// Perform the actual reconstruction 

void Weno3_2D::reconstruct() {
    
    unsigned int no_stencils = 5;
    unsigned int p = 4; 
    //double tau_RHO, tau_RHO_U, tau_RHO_V, tau_E;
	double epsilon = 1.0e-12; 
    double h; // Measure of cell size 
    
    double gammaHi = 0.85;
    double gammaLo = 1.0 - gammaHi;  
    
    Vector<double> gamma(no_stencils); 
    
    gamma(0) = gammaHi; 
    
    for (unsigned int i = 1; i < no_stencils; i++) {
        gamma(i) = 0.25*gammaLo;  
    }

    Vector<double> UL1(4), UL2(4);  Vector<double> WL1(4), WL2(4);

	Vector<double> local_RHO_coeffs(6);
	Vector<double> local_RHO_U_coeffs(6);
	Vector<double> local_RHO_V_coeffs(6);
	Vector<double> local_E_coeffs(6); 

	Vector<double> local_WENO_poly_consts(5);

    // Variables for reconstruction of RHO
    
    double rho0; 
    
    /* Third order centered stencil */ 
    Vector<double> d_rho_3(4); 
    Vector<double> b_rho_3; 
    Vector<double> rho_coeff_3(5); 
    
    Vector<double> rho_coeff_21(2); 
    Vector<double> rho_coeff_22(2); 
    Vector<double> rho_coeff_23(2); 
    Vector<double> rho_coeff_24(2);
    Vector<double> b_rho2(2);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_RHO(no_stencils); Vector<double> w_RHO(no_stencils); double sum_RHO;
    
    // Variables for reconstruction of RHO_U
    double rho_u0; 
    
    /* Third order centered stencil */ 
    Vector<double> d_rho_u_3(4); 
    Vector<double> b_rho_u_3; 
    Vector<double> rho_u_coeff_3(5); 
    
    Vector<double> rho_u_coeff_21(2); 
    Vector<double> rho_u_coeff_22(2); 
    Vector<double> rho_u_coeff_23(2); 
    Vector<double> rho_u_coeff_24(2); 
    Vector<double> b_rho_u_2(2);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_RHO_U(no_stencils); Vector<double> w_RHO_U(no_stencils);  double sum_RHO_U;
    
    // Variables for reconstruction of RHO_V
    double rho_v0;
    
    /* Third order centered stencil */ 
    Vector<double> d_rho_v_3(4); 
    Vector<double> b_rho_v_3; 
    Vector<double> rho_v_coeff_3(5); 
    
    Vector<double> rho_v_coeff_21(2); 
    Vector<double> rho_v_coeff_22(2); 
    Vector<double> rho_v_coeff_23(2); 
    Vector<double> rho_v_coeff_24(2); 
    Vector<double> b_rho_v_2(2);
     
    
    /* Smoothness Indicators */ 
    Vector<double> IS_RHO_V(no_stencils); Vector<double> w_RHO_V(no_stencils); double sum_RHO_V;
    
    // Variables for reconstruction of E
    double e0; 
    
    /* Third order centered stencil */ 
    Vector<double> d_e_3(4); 
    Vector<double> b_e_3; 
    Vector<double> e_coeff_3(5);
    
    Vector<double> e_coeff_21(2); 
    Vector<double> e_coeff_22(2); 
    Vector<double> e_coeff_23(2); 
    Vector<double> e_coeff_24(2);
    Vector<double> b_e_2(2);
    
    /* Smoothness Indicators */ 
    Vector<double> IS_E(no_stencils); Vector<double> w_E(no_stencils); double sum_E; 
    
    // Iterate over all the cells 
    
    DoFHandler<2>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
    
    unsigned int index, ROWS, c, g_i; 
    
    for (; cell != endc; ++cell) {

	if (cell->is_locally_owned()){
        
        cell->get_dof_indices(local_dof_indices);
		c = global_to_local_index_map[local_dof_indices[0] ];
		g_i = local_dof_indices[0];

        rho0   =   RHO(g_i);
        rho_u0 =   RHO_U(g_i);
        rho_v0 =   RHO_V(g_i);
        e0     =   E(g_i);
        h = std::sqrt(Cell[c].measure()); 

		unsigned int neighbor_p1, neighbor_p2, neighbor_p3, neighbor_p4;
        typename DoFHandler<2>::cell_iterator neighbor = cell;

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

        
        if ( !(cell->at_boundary()) ) {
            
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
			
			ROWS = cell_neighbor_iterator[c].size();
			index = 0; 
            
            b_rho_3.reinit(ROWS);
            b_rho_u_3.reinit(ROWS);
            b_rho_v_3.reinit(ROWS);
            b_e_3.reinit(ROWS);
			
				typedef typename std::set<DoFHandler<2>::active_cell_iterator>::iterator ver_cell_iter;

		        typename DoFHandler<2>::cell_iterator neighbor = cell;

				ver_cell_iter diagonal_cell = cell_neighbor_iterator[c].begin();
				ver_cell_iter end_cell = cell_neighbor_iterator[c].end();
				
				// SW vertex

				for (; diagonal_cell != end_cell; ++diagonal_cell) {

					neighbor = *diagonal_cell;
					neighbor->get_dof_indices(local_neighbor_dof_indices);
//					neighbor->get_dof_indices(local_neighbor_dof_indices);

					b_rho_3(index)   = RHO(local_neighbor_dof_indices[0]) - rho0; 
					b_rho_u_3(index) = RHO_U(local_neighbor_dof_indices[0]) - rho_u0;
					b_rho_v_3(index) = RHO_V(local_neighbor_dof_indices[0]) - rho_v0;
					b_e_3(index)     = E(local_neighbor_dof_indices[0]) - e0;
					index++; 
				}
            
			CLS_R3[c].solve(b_rho_3, d_rho_3, rho_coeff_3);
			CLS_R3[c].solve(b_rho_u_3, d_rho_u_3, rho_u_coeff_3); 
			CLS_R3[c].solve(b_rho_v_3, d_rho_v_3, rho_v_coeff_3);
			CLS_R3[c].solve(b_e_3, d_e_3, e_coeff_3);
            
            // =====================================================================
            // r = 2 stencil 1
            // ===================================================================== 

            b_rho2(0)  = (RHO(neighbor_p1) - rho0);   // P1 neighbor 
            b_rho2(1)  = (RHO(neighbor_p2) - rho0);   // P2 neighbor
            LU_R21[c].solve(b_rho2, rho_coeff_21);

            b_rho_u_2(0)  = (RHO_U(neighbor_p1) - rho_u0);   // P1 neighbor 
            b_rho_u_2(1)  = (RHO_U(neighbor_p2) - rho_u0);   // P2 neighbor
            LU_R21[c].solve(b_rho_u_2, rho_u_coeff_21);

            b_rho_v_2(0)  = (RHO_V(neighbor_p1) - rho_v0);   // P1 neighbor 
            b_rho_v_2(1)  = (RHO_V(neighbor_p2) - rho_v0);   // P2 neighbor
            LU_R21[c].solve(b_rho_v_2, rho_v_coeff_21);

            b_e_2(0)  = (E(neighbor_p1) - e0);   // P1 neighbor 
            b_e_2(1)  = (E(neighbor_p2) - e0);   // P2 neighbor
            LU_R21[c].solve(b_e_2, e_coeff_21);
            
            // =====================================================================
            // r = 2 stencil 2 
            // ===================================================================== 

            b_rho2(0) = (RHO(neighbor_p2) - rho0);   // P2 neighbor
            b_rho2(1) = (RHO(neighbor_p3) - rho0);   // P3 neighbor
            LU_R22[c].solve(b_rho2, rho_coeff_22);

            b_rho_u_2(0) = (RHO_U(neighbor_p2) - rho_u0);   // P2 neighbor
            b_rho_u_2(1) = (RHO_U(neighbor_p3) - rho_u0);   // P3 neighbor
            LU_R22[c].solve(b_rho_u_2, rho_u_coeff_22);

            b_rho_v_2(0) = (RHO_V(neighbor_p2) - rho_v0);   // P2 neighbor
            b_rho_v_2(1) = (RHO_V(neighbor_p3) - rho_v0);   // P3 neighbor
            LU_R22[c].solve(b_rho_v_2, rho_v_coeff_22);

            b_e_2(0) = (E(neighbor_p2) - e0);   // P2 neighbor
            b_e_2(1) = (E(neighbor_p3) - e0);   // P3 neighbor
            LU_R22[c].solve(b_e_2, e_coeff_22);
        
  
            // =====================================================================
            // r = 2 stencil 3
            // =====================================================================

            b_rho2(0) = (RHO(neighbor_p3) - rho0);   // P3 neighbor
            b_rho2(1) = (RHO(neighbor_p4) - rho0);   // P4 neighbor                                           // P4 neighbor
            LU_R23[c].solve(b_rho2, rho_coeff_23); 

            b_rho_u_2(0) = (RHO_U(neighbor_p3) - rho_u0);   // P3 neighbor
            b_rho_u_2(1) = (RHO_U(neighbor_p4) - rho_u0);   // P4 neighbor                                   // P4 neighbor
            LU_R23[c].solve(b_rho_u_2, rho_u_coeff_23); 

            b_rho_v_2(0) = (RHO_V(neighbor_p3) - rho_v0);   // P3 neighbor
            b_rho_v_2(1) = (RHO_V(neighbor_p4) - rho_v0);   // P4 neighbor                                   // P4 neighbor
            LU_R23[c].solve(b_rho_v_2, rho_v_coeff_23);

            b_e_2(0) = (E(neighbor_p3) - e0);   // P3 neighbor
            b_e_2(1) = (E(neighbor_p4) - e0);   // P4 neighbor                                               // P4 neighbor
            LU_R23[c].solve(b_e_2, e_coeff_23);  
            

            // =====================================================================
            // r = 2 stencil 4
            // =====================================================================

            b_rho2(0) = (RHO(neighbor_p4) - rho0);   // P4 neighbor
            b_rho2(1) = (RHO(neighbor_p1) - rho0);   // P1 neighbor
            LU_R24[c].solve(b_rho2, rho_coeff_24);

            b_rho_u_2(0) = (RHO_U(neighbor_p4) - rho_u0);   // P4 neighbor
            b_rho_u_2(1) = (RHO_U(neighbor_p1) - rho_u0);   // P1 neighbor
            LU_R24[c].solve(b_rho_u_2, rho_u_coeff_24);

            b_rho_v_2(0) = (RHO_V(neighbor_p4) - rho_v0);   // P4 neighbor
            b_rho_v_2(1) = (RHO_V(neighbor_p1) - rho_v0);   // P1 neighbor
            LU_R24[c].solve(b_rho_v_2, rho_v_coeff_24);

            b_e_2(0) = (E(neighbor_p4) - e0);   // P4 neighbor
            b_e_2(1) = (E(neighbor_p1) - e0);   // P1 neighbor
            LU_R24[c].solve(b_e_2, e_coeff_24);
            
        } // End of interior cell loop 
        
        else {
            
            if (!(is_corner_cell[c])) {
                
                // =====================================================================
                // r = 3 center stencil (boundary)
                // =====================================================================
                
                b_rho_3.reinit(5); b_rho_u_3.reinit(5); b_rho_v_3.reinit(5); b_e_3.reinit(5);  
                
                index = 0; 
                
                if (!cell->face(0)->at_boundary()) {
                    
                    // P1 neighbor
                    
                    b_rho_3(index)    = (RHO(neighbor_p1) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(neighbor_p1) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(neighbor_p1) - rho_v0);                
                    b_e_3(index)      = (E(neighbor_p1) - e0);                
                    
                    index++; 
                }
                
                if (!cell->face(3)->at_boundary()) {

                    // P2 neighbor
                    
                    b_rho_3(index)    = (RHO(neighbor_p2) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(neighbor_p2) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(neighbor_p2) - rho_v0);                
                    b_e_3(index)      = (E(neighbor_p2) - e0);                
                    
                    index++; 
                
                }
                
                if (!cell->face(1)->at_boundary()) {

                    // P3 neighbor
                    
                    b_rho_3(index)    = (RHO(neighbor_p3) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(neighbor_p3) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(neighbor_p3) - rho_v0);                
                    b_e_3(index)      = (E(neighbor_p3) - e0);                
                    
                    index++; 
                }
                
                if (!cell->face(2)->at_boundary()) {

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
                
                // =====================================================================
                // r = 2 stencil 1 (boundary)
                // =====================================================================
                
                if (cell->face(0)->at_boundary()) { // WEST boundary 
                
                    b_rho2(0) = (RHO(neighbor_p2) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p2) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p2) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p2) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 
                }
                
                else if (cell->face(3)->at_boundary()) { // NORTH boundary 
                
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
                
                if (cell->face(1)->at_boundary()) { // EAST boundary 
                
                    b_rho2(0) = (RHO(neighbor_p2) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p2) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p2) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p2) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 
                }
                
                else if (cell->face(3)->at_boundary()) { // NORTH boundary 
                
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
                
                if (cell->face(1)->at_boundary()) { // EAST boundary 
                
                    b_rho2(0) = (RHO(neighbor_p4) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p4) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p4) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p4) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 

                }
                
                else if (cell->face(2)->at_boundary()) { // SOUTH boundary 
                
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
                
                if (cell->face(0)->at_boundary()) { // WEST boundary 
                
                    b_rho2(0) = (RHO(neighbor_p4) - rho0);   // P4 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(neighbor_p4) - rho_u0);   // P4 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(neighbor_p4) - rho_v0);   // P4 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(neighbor_p4) - e0);   // P4 neighbor
                    b_e_2(1) = 0.0; 

                    
                }
                
                else if (cell->face(2)->at_boundary()) { // SOUTH boundary 
                
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
                            // r = 3 
                            // =====================================================================
                            
                            Vector <double> b_3(6); Vector <double> d_3(5); Vector<double> uv_3(10); 
                            
                            index = 0; 
                            
                            
                            if (cell->neighbor_index(0) != -1 ) {
                                b_3(index)   = (RHO_U(neighbor_p1) - rho_u0); // P1 (rho_u)
                                b_3(index+3) = (RHO_V(neighbor_p1) - rho_v0); // P1 (rho_v)
                                index++; 
                            }
                            
                            if (cell->neighbor_index(3) != -1 ) {
                                b_3(index)   = (RHO_U(neighbor_p2) - rho_u0);  // P2 (u)
                                b_3(index+3) = (RHO_V(neighbor_p2) - rho_v0);  // P2 (v)
                                index++; 
                            }
                            
                            if (cell->neighbor_index(1) != -1 ) {
                                b_3(index)   = (RHO_U(neighbor_p3) - rho_u0);  // P3 (u)
                                b_3(index+3) = (RHO_V(neighbor_p3) - rho_v0);  // P3 (v)
                                index++; 
                            }
                            
                            if (cell->neighbor_index(2) != -1 ) {
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
                            
                        }
                    }
                } // End of wall boundary loop 
                
                */ 
                
            } // Non-corner boundary cell loop 
            
            else {
                
                // Corner boundary cells - reduce to first order 
                
                rho_coeff_3 = 0.0;
                rho_coeff_21 = 0.0;  rho_coeff_22 = 0.0;
                rho_coeff_23 = 0.0;  rho_coeff_24 = 0.0;
                
                rho_u_coeff_3 = 0.0;
                rho_u_coeff_21 = 0.0;  rho_u_coeff_22 = 0.0;
                rho_u_coeff_23 = 0.0;  rho_u_coeff_24 = 0.0;
                
                rho_v_coeff_3 = 0.0;
                rho_v_coeff_21 = 0.0;  rho_v_coeff_22 = 0.0;
                rho_v_coeff_23 = 0.0;  rho_v_coeff_24 = 0.0;
                
                e_coeff_3 = 0.0;
                e_coeff_21 = 0.0;  e_coeff_22 = 0.0;
                e_coeff_23 = 0.0;  e_coeff_24 = 0.0;
                
                
            } // Corner boundary cell loop 
            
        } // End of boundary cell loop 
        
        // Find smoothness indicators 
        
        // =====================================================================
        // r = 3
        // =====================================================================
        
        IS_RHO(0)   = compute_third_order_smoothness_indicator(rho_coeff_3, IS_constants[c], h); 
        IS_RHO_U(0) = compute_third_order_smoothness_indicator(rho_u_coeff_3, IS_constants[c], h);
        IS_RHO_V(0) = compute_third_order_smoothness_indicator(rho_v_coeff_3, IS_constants[c], h);
        IS_E(0)     = compute_third_order_smoothness_indicator(e_coeff_3, IS_constants[c], h);
        
        // =====================================================================
        // r = 2  stencil 1
        // =====================================================================
        
        IS_RHO(1)   = compute_second_order_smoothness_indicator(rho_coeff_21, IS_constants[c], h); 
        IS_RHO_U(1) = compute_second_order_smoothness_indicator(rho_u_coeff_21, IS_constants[c], h);
        IS_RHO_V(1) = compute_second_order_smoothness_indicator(rho_v_coeff_21, IS_constants[c], h);
        IS_E(1)     = compute_second_order_smoothness_indicator(e_coeff_21, IS_constants[c], h);
        
        // =====================================================================
        // r = 2  stencil 2
        // =====================================================================
        
        IS_RHO(2)   = compute_second_order_smoothness_indicator(rho_coeff_22, IS_constants[c], h); 
        IS_RHO_U(2) = compute_second_order_smoothness_indicator(rho_u_coeff_22, IS_constants[c], h);
        IS_RHO_V(2) = compute_second_order_smoothness_indicator(rho_v_coeff_22, IS_constants[c], h);
        IS_E(2)     = compute_second_order_smoothness_indicator(e_coeff_22, IS_constants[c], h);
        
        // =====================================================================
        // r = 2  stencil 3
        // =====================================================================
        
        IS_RHO(3)   = compute_second_order_smoothness_indicator(rho_coeff_23, IS_constants[c], h); 
        IS_RHO_U(3) = compute_second_order_smoothness_indicator(rho_u_coeff_23, IS_constants[c], h);
        IS_RHO_V(3) = compute_second_order_smoothness_indicator(rho_v_coeff_23, IS_constants[c], h);
        IS_E(3)     = compute_second_order_smoothness_indicator(e_coeff_23, IS_constants[c], h);
        
        // =====================================================================
        // r = 2  stencil 4
        // =====================================================================
        
        IS_RHO(4)   = compute_second_order_smoothness_indicator(rho_coeff_24, IS_constants[c], h); 
        IS_RHO_U(4) = compute_second_order_smoothness_indicator(rho_u_coeff_24, IS_constants[c], h);
        IS_RHO_V(4) = compute_second_order_smoothness_indicator(rho_v_coeff_24, IS_constants[c], h);
        IS_E(4)     = compute_second_order_smoothness_indicator(e_coeff_24, IS_constants[c], h);
    
        
        // Combine the polynomials 
        
        sum_RHO = 0.0; sum_RHO_U = 0.0; sum_RHO_V = 0.0; sum_E = 0.0;
		/*
		tau_RHO = (0.25)*(std::abs(IS_RHO[0] - IS_RHO[1]) + std::abs(IS_RHO[0] - IS_RHO[2]) + 
							std::abs(IS_RHO[0] - IS_RHO[3]) + std::abs(IS_RHO[0] - IS_RHO[4]) );
		
		tau_RHO = tau_RHO*tau_RHO; 
		
		tau_RHO_U = (0.25)*(std::abs(IS_RHO_U[0] - IS_RHO_U[1]) + std::abs(IS_RHO_U[0] - IS_RHO_U[2]) + 
							std::abs(IS_RHO_U[0] - IS_RHO_U[3]) + std::abs(IS_RHO_U[0] - IS_RHO_U[4]) );
		
		tau_RHO_U = tau_RHO_U*tau_RHO_U;
		
		tau_RHO_V = (0.25)*(std::abs(IS_RHO_V[0] - IS_RHO_V[1]) + std::abs(IS_RHO_V[0] - IS_RHO_V[2]) + 
							std::abs(IS_RHO_V[0] - IS_RHO_V[3]) + std::abs(IS_RHO_V[0] - IS_RHO_V[4]) );
		
		tau_RHO_V = tau_RHO_V*tau_RHO_V;
		
		tau_E = (0.25)*(std::abs(IS_E[0] - IS_E[1]) + std::abs(IS_E[0] - IS_E[2]) + 
						std::abs(IS_E[0] - IS_E[3]) + std::abs(IS_E[0] - IS_E[4]) );
		
		tau_E = tau_E*tau_E;
        */ 
        for (unsigned int j = 0; j < no_stencils; j++ ) {
			
			//w_RHO(j)   = gamma(j)*(1.0 + (tau_RHO/(epsilon+IS_RHO[j])) );
			//w_RHO_U(j) = gamma(j)*(1.0 + (tau_RHO_U/(epsilon+IS_RHO_U[j])) );
			//w_RHO_V(j) = gamma(j)*(1.0 + (tau_RHO_V/(epsilon+IS_RHO_V[j])) );
			//w_E(j)     = gamma(j)*(1.0 + (tau_E/(epsilon+IS_E[j])) );
            
            w_RHO(j)   = gamma(j)/(std::pow((IS_RHO(j) + epsilon), p));
            w_RHO_U(j) = gamma(j)/(std::pow((IS_RHO_U(j) + epsilon), p));
            w_RHO_V(j) = gamma(j)/(std::pow((IS_RHO_V(j) + epsilon), p));
            w_E(j)     = gamma(j)/(std::pow((IS_E(j) + epsilon), p));
            
            sum_RHO += w_RHO(j); sum_RHO_U += w_RHO_U(j); sum_RHO_V += w_RHO_V(j); sum_E += w_E(j);
        }
        
        for (unsigned int j = 0; j < no_stencils; j++ ) {
            w_RHO(j) = w_RHO(j)/sum_RHO; 
            w_RHO_U(j) = w_RHO_U(j)/sum_RHO_U;
            w_RHO_V(j) = w_RHO_V(j)/sum_RHO_V;
            w_E(j) = w_E(j)/sum_E;
        }
    
        // Density 
        
        local_coeffs_x_RHO(g_i) = (w_RHO(0)/gamma(0))*(rho_coeff_3(0) - 
                                                gamma(1)*rho_coeff_21(0) -
                                                gamma(2)*rho_coeff_22(0) -
                                                gamma(3)*rho_coeff_23(0) -
                                                gamma(4)*rho_coeff_24(0) ) + 
                                                
                                                w_RHO(1)*rho_coeff_21(0) +
                                                w_RHO(2)*rho_coeff_22(0) +
                                                w_RHO(3)*rho_coeff_23(0) +
                                                w_RHO(4)*rho_coeff_24(0) ; // u_x
        
        local_coeffs_y_RHO(g_i) = (w_RHO(0)/gamma(0))*(rho_coeff_3(1) - 
                                                gamma(1)*rho_coeff_21(1) -
                                                gamma(2)*rho_coeff_22(1) -
                                                gamma(3)*rho_coeff_23(1) -
                                                gamma(4)*rho_coeff_24(1) ) + 
                                                
                                                w_RHO(1)*rho_coeff_21(1) +
                                                w_RHO(2)*rho_coeff_22(1) +
                                                w_RHO(3)*rho_coeff_23(1) +
                                                w_RHO(4)*rho_coeff_24(1) ; // u_y

        local_coeffs_xx_RHO(g_i) = (w_RHO(0)/gamma(0))*(rho_coeff_3(2)); // u_xx
                                        
                                                
        local_coeffs_yy_RHO(g_i) = (w_RHO(0)/gamma(0))*(rho_coeff_3(3)); // u_yy
        
        local_coeffs_xy_RHO(g_i) = (w_RHO(0)/gamma(0))*(rho_coeff_3(4)); // u_xy
        
        
        // x-momentum
        
        local_coeffs_x_RHO_U(g_i) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(0) - 
                                                gamma(1)*rho_u_coeff_21(0) -
                                                gamma(2)*rho_u_coeff_22(0) -
                                                gamma(3)*rho_u_coeff_23(0) -
                                                gamma(4)*rho_u_coeff_24(0) ) + 
                                                
                                                w_RHO_U(1)*rho_u_coeff_21(0) +
                                                w_RHO_U(2)*rho_u_coeff_22(0) +
                                                w_RHO_U(3)*rho_u_coeff_23(0) +
                                                w_RHO_U(4)*rho_u_coeff_24(0) ; // u_x
        
        local_coeffs_y_RHO_U(g_i) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(1) - 
                                                gamma(1)*rho_u_coeff_21(1) -
                                                gamma(2)*rho_u_coeff_22(1) -
                                                gamma(3)*rho_u_coeff_23(1) -
                                                gamma(4)*rho_u_coeff_24(1) ) + 
                                                
                                                w_RHO_U(1)*rho_u_coeff_21(1) +
                                                w_RHO_U(2)*rho_u_coeff_22(1) +
                                                w_RHO_U(3)*rho_u_coeff_23(1) +
                                                w_RHO_U(4)*rho_u_coeff_24(1) ; // u_y

        local_coeffs_xx_RHO_U(g_i) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(2)); // u_xx
                                        
                                                
        local_coeffs_yy_RHO_U(g_i) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(3)); // u_yy
        
        local_coeffs_xy_RHO_U(g_i) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(4)); // u_xy
        
        
        // y-momentum       
        
        local_coeffs_x_RHO_V(g_i) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(0) - 
                                                gamma(1)*rho_v_coeff_21(0) -
                                                gamma(2)*rho_v_coeff_22(0) -
                                                gamma(3)*rho_v_coeff_23(0) -
                                                gamma(4)*rho_v_coeff_24(0) ) + 
                                                
                                                w_RHO_V(1)*rho_v_coeff_21(0) +
                                                w_RHO_V(2)*rho_v_coeff_22(0) +
                                                w_RHO_V(3)*rho_v_coeff_23(0) +
                                                w_RHO_V(4)*rho_v_coeff_24(0) ; // v_x
        
        local_coeffs_y_RHO_V(g_i) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(1) - 
                                                gamma(1)*rho_v_coeff_21(1) -
                                                gamma(2)*rho_v_coeff_22(1) -
                                                gamma(3)*rho_v_coeff_23(1) -
                                                gamma(4)*rho_v_coeff_24(1) ) + 
                                                
                                                w_RHO_V(1)*rho_v_coeff_21(1) +
                                                w_RHO_V(2)*rho_v_coeff_22(1) +
                                                w_RHO_V(3)*rho_v_coeff_23(1) +
                                                w_RHO_V(4)*rho_v_coeff_24(1) ; // v_y

        local_coeffs_xx_RHO_V(g_i) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(2)); // v_xx
                                        
                                                
        local_coeffs_yy_RHO_V(g_i) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(3)); // v_yy
        
        local_coeffs_xy_RHO_V(g_i) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(4)); // v_xy
        
        // Total energy  
               
        local_coeffs_x_E(g_i) = (w_E(0)/gamma(0))*(e_coeff_3(0) - 
                                                gamma(1)*e_coeff_21(0) -
                                                gamma(2)*e_coeff_22(0) -
                                                gamma(3)*e_coeff_23(0) -
                                                gamma(4)*e_coeff_24(0) ) + 
                                                
                                                w_E(1)*e_coeff_21(0) +
                                                w_E(2)*e_coeff_22(0) +
                                                w_E(3)*e_coeff_23(0) +
                                                w_E(4)*e_coeff_24(0) ; // v_x
        
        local_coeffs_y_E(g_i) = (w_E(0)/gamma(0))*(e_coeff_3(1) - 
                                                gamma(1)*e_coeff_21(1) -
                                                gamma(2)*e_coeff_22(1) -
                                                gamma(3)*e_coeff_23(1) -
                                                gamma(4)*e_coeff_24(1) ) + 
                                                
                                                w_E(1)*e_coeff_21(1) +
                                                w_E(2)*e_coeff_22(1) +
                                                w_E(3)*e_coeff_23(1) +
                                                w_E(4)*e_coeff_24(1) ; // v_y

        local_coeffs_xx_E(g_i) = (w_E(0)/gamma(0))*(e_coeff_3(2)); // v_xx
                                        
                                                
        local_coeffs_yy_E(g_i) = (w_E(0)/gamma(0))*(e_coeff_3(3)); // v_yy
        
        local_coeffs_xy_E(g_i) = (w_E(0)/gamma(0))*(e_coeff_3(4)); // v_xy
        /*
        
    	    for (unsigned int f = 0; f < 4; ++f) {  
	
				local_RHO_coeffs(0) = RHO(g_i);
				local_RHO_U_coeffs(0) = RHO_U(g_i);
				local_RHO_V_coeffs(0) = RHO_V(g_i);
				local_E_coeffs(0) = E(g_i);

				local_RHO_coeffs(1) = local_coeffs_x_RHO(g_i);
				local_RHO_U_coeffs(1) = local_coeffs_x_RHO_U(g_i);
				local_RHO_V_coeffs(1) = local_coeffs_x_RHO_V(g_i);
				local_E_coeffs(1) = local_coeffs_x_E(g_i);
	
				local_RHO_coeffs(2) = local_coeffs_y_RHO(g_i);
				local_RHO_U_coeffs(2) = local_coeffs_y_RHO_U(g_i);
				local_RHO_V_coeffs(2) = local_coeffs_y_RHO_V(g_i);
				local_E_coeffs(2) = local_coeffs_y_E(g_i); 
	
				local_RHO_coeffs(3) = local_coeffs_xx_RHO(g_i);
				local_RHO_U_coeffs(3) = local_coeffs_xx_RHO_U(g_i);
				local_RHO_V_coeffs(3) = local_coeffs_xx_RHO_V(g_i);
				local_E_coeffs(3) = local_coeffs_xx_E(g_i);

				local_RHO_coeffs(4) = local_coeffs_yy_RHO(g_i);
				local_RHO_U_coeffs(4) = local_coeffs_yy_RHO_U(g_i);
				local_RHO_V_coeffs(4) = local_coeffs_yy_RHO_V(g_i);
				local_E_coeffs(4) = local_coeffs_yy_E(g_i);
	
				local_RHO_coeffs(5) = local_coeffs_xy_RHO(g_i);
				local_RHO_U_coeffs(5) = local_coeffs_xy_RHO_U(g_i);
				local_RHO_V_coeffs(5) = local_coeffs_xy_RHO_V(g_i);
				local_E_coeffs(5) = local_coeffs_xy_E(g_i); 
	
				local_WENO_poly_consts(0) = WENO_poly_consts_x(g_i);
				local_WENO_poly_consts(1) = WENO_poly_consts_y(g_i);
				local_WENO_poly_consts(2) = WENO_poly_consts_xx(g_i);
				local_WENO_poly_consts(3) = WENO_poly_consts_yy(g_i);
				local_WENO_poly_consts(4) = WENO_poly_consts_xy(g_i);

                Point<2> face_quadrature_point_1 = Cell[c].face_quadrature_point1(f);
                Point<2> face_quadrature_point_2 = Cell[c].face_quadrature_point2(f);

                // Left face
                UL1(0) = evaluate_weno_polynomial(local_RHO_coeffs, local_WENO_poly_consts, face_quadrature_point_1);
                UL1(1) = evaluate_weno_polynomial(local_RHO_U_coeffs, local_WENO_poly_consts, face_quadrature_point_1);
                UL1(2) = evaluate_weno_polynomial(local_RHO_V_coeffs, local_WENO_poly_consts, face_quadrature_point_1);
                UL1(3) = evaluate_weno_polynomial(local_E_coeffs, local_WENO_poly_consts, face_quadrature_point_1);
                
                UL2(0) = evaluate_weno_polynomial(local_RHO_coeffs, local_WENO_poly_consts, face_quadrature_point_2);
                UL2(1) = evaluate_weno_polynomial(local_RHO_U_coeffs, local_WENO_poly_consts, face_quadrature_point_2);
                UL2(2) = evaluate_weno_polynomial(local_RHO_V_coeffs, local_WENO_poly_consts, face_quadrature_point_2);
                UL2(3) = evaluate_weno_polynomial(local_E_coeffs, local_WENO_poly_consts, face_quadrature_point_2);
               
                WL1 = conserved_to_primitive(UL1); WL2 = conserved_to_primitive(UL2);

    	        if (WL1(0) < 0.0 || WL1(3) < 0.0 || WL2(0) < 0.0 || WL2(3) < 0.0 ) {

					local_coeffs_x_RHO(g_i) = 0.0;	
					local_coeffs_x_RHO_U(g_i) = 0.0;	
					local_coeffs_x_RHO_V(g_i) = 0.0;		
					local_coeffs_x_E(g_i) = 0.0;	

					local_coeffs_y_RHO(g_i) = 0.0;	
					local_coeffs_y_RHO_U(g_i) = 0.0;	
					local_coeffs_y_RHO_V(g_i) = 0.0;		
					local_coeffs_y_E(g_i) = 0.0;

					local_coeffs_xx_RHO(g_i) = 0.0;	
					local_coeffs_xx_RHO_U(g_i) = 0.0;	
					local_coeffs_xx_RHO_V(g_i) = 0.0;		
					local_coeffs_xx_E(g_i) = 0.0;	

					local_coeffs_yy_RHO(g_i) = 0.0;	
					local_coeffs_yy_RHO_U(g_i) = 0.0;	
					local_coeffs_yy_RHO_V(g_i) = 0.0;		
					local_coeffs_yy_E(g_i) = 0.0;

					local_coeffs_xy_RHO(g_i) = 0.0;	
					local_coeffs_xy_RHO_U(g_i) = 0.0;	
					local_coeffs_xy_RHO_V(g_i) = 0.0;		
					local_coeffs_xy_E(g_i) = 0.0;

    	        }
			}
			*/ 

        }        
    } // End of cell loop 

	local_coeffs_x_RHO.compress(VectorOperation::insert);
	local_coeffs_x_RHO_U.compress(VectorOperation::insert);
	local_coeffs_x_RHO_V.compress(VectorOperation::insert);
	local_coeffs_x_E.compress(VectorOperation::insert);

	local_coeffs_y_RHO.compress(VectorOperation::insert);
	local_coeffs_y_RHO_U.compress(VectorOperation::insert);
	local_coeffs_y_RHO_V.compress(VectorOperation::insert);
	local_coeffs_y_E.compress(VectorOperation::insert);

	local_coeffs_xx_RHO.compress(VectorOperation::insert);
	local_coeffs_xx_RHO_U.compress(VectorOperation::insert);
	local_coeffs_xx_RHO_V.compress(VectorOperation::insert);
	local_coeffs_xx_E.compress(VectorOperation::insert);

	local_coeffs_yy_RHO.compress(VectorOperation::insert);
	local_coeffs_yy_RHO_U.compress(VectorOperation::insert);
	local_coeffs_yy_RHO_V.compress(VectorOperation::insert);
	local_coeffs_yy_E.compress(VectorOperation::insert);

	local_coeffs_xy_RHO.compress(VectorOperation::insert);
	local_coeffs_xy_RHO_U.compress(VectorOperation::insert);
	local_coeffs_xy_RHO_V.compress(VectorOperation::insert);
	local_coeffs_xy_E.compress(VectorOperation::insert);

	coeffs_x_RHO = local_coeffs_x_RHO;
	coeffs_x_RHO_U = local_coeffs_x_RHO_U;
	coeffs_x_RHO_V = local_coeffs_x_RHO_V;
	coeffs_x_E = local_coeffs_x_E;

	coeffs_y_RHO = local_coeffs_y_RHO;
	coeffs_y_RHO_U = local_coeffs_y_RHO_U;
	coeffs_y_RHO_V = local_coeffs_y_RHO_V;
	coeffs_y_E = local_coeffs_y_E;

	coeffs_xx_RHO = local_coeffs_xx_RHO;
	coeffs_xx_RHO_U = local_coeffs_xx_RHO_U;
	coeffs_xx_RHO_V = local_coeffs_xx_RHO_V;
	coeffs_xx_E = local_coeffs_xx_E;

	coeffs_yy_RHO = local_coeffs_yy_RHO;
	coeffs_yy_RHO_U = local_coeffs_yy_RHO_U;
	coeffs_yy_RHO_V = local_coeffs_yy_RHO_V;
	coeffs_yy_E = local_coeffs_yy_E;

	coeffs_xy_RHO = local_coeffs_xy_RHO;
	coeffs_xy_RHO_U = local_coeffs_xy_RHO_U;
	coeffs_xy_RHO_V = local_coeffs_xy_RHO_V;
	coeffs_xy_E = local_coeffs_xy_E; 
   
} // End of function 



/*             Test if fourth order stencil Transmissive boundary conditions are applied properly 
 * 
 * 
 * 
 * 
                unsigned int f; 

                // Find the boundary face 
                for (int i = 0; i < 4; i++) {
                    if (cell->face(i)->at_boundary()) {
                        f = i; 
                    }
                }
                
                QGauss<2-1> face_quadrature_formula(1);
                FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);
                
                fv_face_values.reinit(cell, f);
                
                double nx, ny; 
                Tensor<1,2> face_normal_vector1; // Face normal vector
                face_normal_vector1 = fv_face_values.normal_vector(0); 
                nx = face_normal_vector1[0]; ny = face_normal_vector1[1];
                
                Vector<double> coeffs(10); 
                
                coeffs(0) = rho0; 
                
                for (unsigned int i =1; i<10; i++) {
                    coeffs(i) = rho_coeff_4(i-1); 
                }
                
                
                Point<2> C; 
                
                C = cell->face(f)->center();
                double xg = C(0); double yg = C(1); 
                
                double value = nx*(rho_coeff_4(0) + 2*xg*rho_coeff_4(2) + yg*rho_coeff_4(4) + 
                                   3.0*xg*xg*rho_coeff_4(5) +      2*xg*yg*rho_coeff_4(7) + yg*yg*rho_coeff_4(8)) + 
                               ny*(rho_coeff_4(1) + 2*yg*rho_coeff_4(3) + xg*rho_coeff_4(4) + 
                                   3.0*xg*xg*rho_coeff_4(6) + xg*xg*rho_coeff_4(7) + 2*xg*yg*rho_coeff_4(8)); 
                
                double value2 = nx*rho_coeff_4(5) + ny*rho_coeff_4(6); 
                
                std::cout << "nx = " << nx << ", ny = " << ny << std::endl; 
                std::cout << "First Derivative = " << value << std::endl; 
                std::cout << "Third Derivative = " << value2 << std::endl;
                std::cout << "rho_coeff = " << rho_coeff_4 << std::endl;
                std::cout << "Poly_value = " << evaluate_weno_polynomial(coeffs, WENO_poly_consts[c], cell->face(f)->center()) << std::endl;
                
                std::cout << "==================================================================================" << std::endl;
*/ 



/* 
 *  Test if fourth order slip condition is working well 
 * 
 *                             //Point<2> face_center; 
                            
                            //face_center = cell->face(f)->center();
                            
                            //Vector<double> U_Coeffs(10); Vector<double> V_Coeffs(10); 
                            
                            //U_Coeffs(0) = rho_u0; V_Coeffs(0) = rho_v0; 

                            for (unsigned int i = 0; i < 9; i++) {
                                U_Coeffs[i+1] = uv_4[i];  
                                V_Coeffs[i+1] = uv_4[i+9];
                            }

                            //std::cout << U_Coeffs; 
                            
                            //double un = nx*evaluate_weno_polynomial(U_Coeffs, WENO_poly_consts[c], face_center) + 
                                        //ny*evaluate_weno_polynomial(V_Coeffs, WENO_poly_consts[c], face_center);
                                        
                            //std::cout << "Normal velocity = " << un << std::endl; 
*/

/* Test if third order slip condition is working well 
 *                             Vector<double> U_Coeffs(6); Vector<double> V_Coeffs(6); 
                            
                            U_Coeffs(0) = rho_u0; V_Coeffs(0) = rho_v0; 
                            
                            for (unsigned int i = 0; i < 5; i++) {
                                U_Coeffs[i+1] = rho_u_coeff_3[i]; 
                                V_Coeffs[i+1] = rho_v_coeff_3[i];
                            } 

                            //double un = nx*evaluate_weno_polynomial(U_Coeffs, WENO_poly_consts[c], face_center) + 
                            //         ny*evaluate_weno_polynomial(V_Coeffs, WENO_poly_consts[c], face_center);
                                        
                            //double dut_dn = -nx*ny*(u_coeff_3(0) + 2.0*face_center(0)*u_coeff_3(2) + face_center(1)*u_coeff_3(4)) +
                            //             nx*nx*(v_coeff_3(0) + 2.0*face_center(0)*v_coeff_3(2) + face_center(1)*v_coeff_3(4)) - 
                            //           ny*ny*(u_coeff_3(1) + 2.0*face_center(1)*u_coeff_3(3) + face_center(0)*u_coeff_3(4)) + 
                            //          ny*nx*(v_coeff_3(1) + 2.0*face_center(1)*v_coeff_3(3) + face_center(0)*v_coeff_3(4)) ; 
                            
                            
                            //std::cout << "un3 = " << un << std::endl; 
                            //std::cout << "dut_dn3 = " << dut_dn << std::endl; 
*/

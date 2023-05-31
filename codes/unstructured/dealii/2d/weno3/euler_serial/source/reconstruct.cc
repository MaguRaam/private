#include "../include/Weno32.h"

// Perform the actual reconstruction 

void Weno3_2D::reconstruct() {
    
    unsigned int no_stencils = 5;
    unsigned int p = 2; 
    double epsilon = 1.0e-12; 
    double h; // Measure of cell size 
    
    double gammaHi = 0.8;
    double gammaLo = 1.0 - gammaHi;  
    
    Vector<double> gamma(no_stencils); 
    
    gamma(0) = gammaHi; 
    
    for (unsigned int i = 1; i < no_stencils; i++) {
        gamma(i) = 0.25*gammaLo;  
    }

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
    
    unsigned int index, ROWS; 
    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        rho0   =   RHO(cell->active_cell_index());
        rho_u0 =   RHO_U(cell->active_cell_index());
        rho_v0 =   RHO_V(cell->active_cell_index());
        e0     =   E(cell->active_cell_index());
        h = std::sqrt(cell->measure()); 
        
        if ( !(cell->at_boundary()) ) {
            
            // =====================================================================
            // r = 3 stencil (Centered stencil)
            // =====================================================================
            
			// constraint part (consists of face neighbours)
			
			d_rho_3(0)  = (RHO(Stencil[c].W_index) - rho0);                // W neighbor 
			d_rho_3(1)  = (RHO(Stencil[c].N_index) - rho0);                // N neighbor
			d_rho_3(2)  = (RHO(Stencil[c].E_index) - rho0);                // E neighbor
			d_rho_3(3)  = (RHO(Stencil[c].S_index) - rho0);                // S neighbor
			
			d_rho_u_3(0)  = (RHO_U(Stencil[c].W_index) - rho_u0);          // W neighbor 
			d_rho_u_3(1)  = (RHO_U(Stencil[c].N_index) - rho_u0);          // N neighbor
			d_rho_u_3(2)  = (RHO_U(Stencil[c].E_index) - rho_u0);          // E neighbor
			d_rho_u_3(3)  = (RHO_U(Stencil[c].S_index) - rho_u0);          // S neighbor
			
			d_rho_v_3(0)  = (RHO_V(Stencil[c].W_index) - rho_v0);          // W neighbor 
			d_rho_v_3(1)  = (RHO_V(Stencil[c].N_index) - rho_v0);          // N neighbor
			d_rho_v_3(2)  = (RHO_V(Stencil[c].E_index) - rho_v0);          // E neighbor
			d_rho_v_3(3)  = (RHO_V(Stencil[c].S_index) - rho_v0);          // S neighbor
			
			d_e_3(0)  = (E(Stencil[c].W_index) - e0);                      // W neighbor 
			d_e_3(1)  = (E(Stencil[c].N_index) - e0);                      // N neighbor
			d_e_3(2)  = (E(Stencil[c].E_index) - e0);                      // E neighbor
			d_e_3(3)  = (E(Stencil[c].S_index) - e0);                      // S neighbor
			
            ROWS = Stencil[c].SW_index.size() + Stencil[c].SE_index.size() + Stencil[c].NW_index.size() + Stencil[c].NE_index.size(); 
			index = 0; 
            
            b_rho_3.reinit(ROWS);
            b_rho_u_3.reinit(ROWS);
            b_rho_v_3.reinit(ROWS);
            b_e_3.reinit(ROWS);
			
			// SW vertex

			for (unsigned int i = 0; i < Stencil[c].SW_index.size(); ++i) {
				b_rho_3(index)   = RHO(Stencil[c].SW_index[i]) - rho0; 
				b_rho_u_3(index) = RHO_U(Stencil[c].SW_index[i]) - rho_u0;
				b_rho_v_3(index) = RHO_V(Stencil[c].SW_index[i]) - rho_v0;
				b_e_3(index)     = E(Stencil[c].SW_index[i]) - e0;
				index++; 
			}
			
			// SE vertex 
			
			for (unsigned int i = 0; i < Stencil[c].SE_index.size(); ++i) {
				b_rho_3(index)   = RHO(Stencil[c].SE_index[i]) - rho0; 
				b_rho_u_3(index) = RHO_U(Stencil[c].SE_index[i]) - rho_u0;
				b_rho_v_3(index) = RHO_V(Stencil[c].SE_index[i]) - rho_v0;
				b_e_3(index)     = E(Stencil[c].SE_index[i]) - e0;
				index++; 
			}
			
			// NW vertex 
			
			for (unsigned int i = 0; i < Stencil[c].NW_index.size(); ++i) {
				b_rho_3(index)   = RHO(Stencil[c].NW_index[i]) - rho0; 
				b_rho_u_3(index) = RHO_U(Stencil[c].NW_index[i]) - rho_u0;
				b_rho_v_3(index) = RHO_V(Stencil[c].NW_index[i]) - rho_v0;
				b_e_3(index)     = E(Stencil[c].NW_index[i]) - e0;
				index++; 
			}
			
			// NE vertex 
			
			for (unsigned int i = 0; i < Stencil[c].NE_index.size(); ++i) {
				b_rho_3(index)   = RHO(Stencil[c].NE_index[i]) - rho0; 
				b_rho_u_3(index) = RHO_U(Stencil[c].NE_index[i]) - rho_u0;
				b_rho_v_3(index) = RHO_V(Stencil[c].NE_index[i]) - rho_v0;
				b_e_3(index)     = E(Stencil[c].NE_index[i]) - e0;
			}
            
			CLS_R3[c].solve(b_rho_3, d_rho_3, rho_coeff_3);
			CLS_R3[c].solve(b_rho_u_3, d_rho_u_3, rho_u_coeff_3); 
			CLS_R3[c].solve(b_rho_v_3, d_rho_v_3, rho_v_coeff_3);
			CLS_R3[c].solve(b_e_3, d_e_3, e_coeff_3);
            
            // =====================================================================
            // r = 2 stencil 1
            // ===================================================================== 

            b_rho2(0)  = (RHO(cell->neighbor_index(0)) - rho0);   // P1 neighbor 
            b_rho2(1)  = (RHO(cell->neighbor_index(3)) - rho0);   // P2 neighbor
            LU_R21[c].solve(b_rho2, rho_coeff_21);

            b_rho_u_2(0)  = (RHO_U(cell->neighbor_index(0)) - rho_u0);   // P1 neighbor 
            b_rho_u_2(1)  = (RHO_U(cell->neighbor_index(3)) - rho_u0);   // P2 neighbor
            LU_R21[c].solve(b_rho_u_2, rho_u_coeff_21);

            b_rho_v_2(0)  = (RHO_V(cell->neighbor_index(0)) - rho_v0);   // P1 neighbor 
            b_rho_v_2(1)  = (RHO_V(cell->neighbor_index(3)) - rho_v0);   // P2 neighbor
            LU_R21[c].solve(b_rho_v_2, rho_v_coeff_21);

            b_e_2(0)  = (E(cell->neighbor_index(0)) - e0);   // P1 neighbor 
            b_e_2(1)  = (E(cell->neighbor_index(3)) - e0);   // P2 neighbor
            LU_R21[c].solve(b_e_2, e_coeff_21);
            
            // =====================================================================
            // r = 2 stencil 2 
            // ===================================================================== 

            b_rho2(0) = (RHO(cell->neighbor_index(3)) - rho0);   // P2 neighbor
            b_rho2(1) = (RHO(cell->neighbor_index(1)) - rho0);   // P3 neighbor
            LU_R22[c].solve(b_rho2, rho_coeff_22);

            b_rho_u_2(0) = (RHO_U(cell->neighbor_index(3)) - rho_u0);   // P2 neighbor
            b_rho_u_2(1) = (RHO_U(cell->neighbor_index(1)) - rho_u0);   // P3 neighbor
            LU_R22[c].solve(b_rho_u_2, rho_u_coeff_22);

            b_rho_v_2(0) = (RHO_V(cell->neighbor_index(3)) - rho_v0);   // P2 neighbor
            b_rho_v_2(1) = (RHO_V(cell->neighbor_index(1)) - rho_v0);   // P3 neighbor
            LU_R22[c].solve(b_rho_v_2, rho_v_coeff_22);

            b_e_2(0) = (E(cell->neighbor_index(3)) - e0);   // P2 neighbor
            b_e_2(1) = (E(cell->neighbor_index(1)) - e0);   // P3 neighbor
            LU_R22[c].solve(b_e_2, e_coeff_22);
        
  
            // =====================================================================
            // r = 2 stencil 3
            // =====================================================================

            b_rho2(0) = (RHO(cell->neighbor_index(1)) - rho0);   // P3 neighbor
            b_rho2(1) = (RHO(cell->neighbor_index(2)) - rho0);   // P4 neighbor                                           // P4 neighbor
            LU_R23[c].solve(b_rho2, rho_coeff_23); 

            b_rho_u_2(0) = (RHO_U(cell->neighbor_index(1)) - rho_u0);   // P3 neighbor
            b_rho_u_2(1) = (RHO_U(cell->neighbor_index(2)) - rho_u0);   // P4 neighbor                                   // P4 neighbor
            LU_R23[c].solve(b_rho_u_2, rho_u_coeff_23); 

            b_rho_v_2(0) = (RHO_V(cell->neighbor_index(1)) - rho_v0);   // P3 neighbor
            b_rho_v_2(1) = (RHO_V(cell->neighbor_index(2)) - rho_v0);   // P4 neighbor                                   // P4 neighbor
            LU_R23[c].solve(b_rho_v_2, rho_v_coeff_23);

            b_e_2(0) = (E(cell->neighbor_index(1)) - e0);   // P3 neighbor
            b_e_2(1) = (E(cell->neighbor_index(2)) - e0);   // P4 neighbor                                               // P4 neighbor
            LU_R23[c].solve(b_e_2, e_coeff_23);  
            

            // =====================================================================
            // r = 2 stencil 4
            // =====================================================================

            b_rho2(0) = (RHO(cell->neighbor_index(2)) - rho0);   // P4 neighbor
            b_rho2(1) = (RHO(cell->neighbor_index(0)) - rho0);   // P1 neighbor
            LU_R24[c].solve(b_rho2, rho_coeff_24);

            b_rho_u_2(0) = (RHO_U(cell->neighbor_index(2)) - rho_u0);   // P4 neighbor
            b_rho_u_2(1) = (RHO_U(cell->neighbor_index(0)) - rho_u0);   // P1 neighbor
            LU_R24[c].solve(b_rho_u_2, rho_u_coeff_24);

            b_rho_v_2(0) = (RHO_V(cell->neighbor_index(2)) - rho_v0);   // P4 neighbor
            b_rho_v_2(1) = (RHO_V(cell->neighbor_index(0)) - rho_v0);   // P1 neighbor
            LU_R24[c].solve(b_rho_v_2, rho_v_coeff_24);

            b_e_2(0) = (E(cell->neighbor_index(2)) - e0);   // P4 neighbor
            b_e_2(1) = (E(cell->neighbor_index(0)) - e0);   // P1 neighbor
            LU_R24[c].solve(b_e_2, e_coeff_24);
            
        } // End of interior cell loop 
        
        else {
            
            if (!(is_corner_cell[c])) {
                
                // =====================================================================
                // r = 3 center stencil (boundary)
                // =====================================================================
                
                b_rho_3.reinit(5); b_rho_u_3.reinit(5); b_rho_v_3.reinit(5); b_e_3.reinit(5);  
                
                index = 0; 
                
                if (Stencil[c].is_W_present) {
                    
                    // P1 neighbor
                    
                    b_rho_3(index)    = (RHO(Stencil[c].W_index) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(Stencil[c].W_index) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(Stencil[c].W_index) - rho_v0);                
                    b_e_3(index)      = (E(Stencil[c].W_index) - e0);                
                    
                    index++; 
                }
                
                if (Stencil[c].is_N_present) {

                    // P2 neighbor
                    
                    b_rho_3(index)    = (RHO(Stencil[c].N_index) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(Stencil[c].N_index) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(Stencil[c].N_index) - rho_v0);                
                    b_e_3(index)      = (E(Stencil[c].N_index) - e0);                
                    
                    index++; 
                
                }
                
                if (Stencil[c].is_E_present) {

                    // P3 neighbor
                    
                    b_rho_3(index)    = (RHO(Stencil[c].E_index) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(Stencil[c].E_index) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(Stencil[c].E_index) - rho_v0);                
                    b_e_3(index)      = (E(Stencil[c].E_index) - e0);                
                    
                    index++; 
                }
                
                if (Stencil[c].is_S_present) {

                    // P4 neighbor
                    
                    b_rho_3(index)    = (RHO(Stencil[c].S_index) - rho0);                 
                    b_rho_u_3(index)  = (RHO_U(Stencil[c].S_index) - rho_u0);                
                    b_rho_v_3(index)  = (RHO_V(Stencil[c].S_index) - rho_v0);                
                    b_e_3(index)      = (E(Stencil[c].S_index) - e0);                
                    
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
                
                if (cell->neighbor_index(0) == -1) { // WEST boundary 
                
                    b_rho2(0) = (RHO(cell->neighbor_index(3)) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(3)) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(3)) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(cell->neighbor_index(3)) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 
                }
                
                else if (cell->neighbor_index(3) == -1) { // NORTH boundary 
                
                    b_rho2(0) = (RHO(cell->neighbor_index(0)) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(0)) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(0)) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(cell->neighbor_index(0)) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 
                }
                
                else {
                    
                    b_rho2(0)  = (RHO(cell->neighbor_index(0)) - rho0);   // P1 neighbor 
                    b_rho2(1)  = (RHO(cell->neighbor_index(3)) - rho0);   // P2 neighbor

                    b_rho_u_2(0)  = (RHO_U(cell->neighbor_index(0)) - rho_u0);   // P1 neighbor 
                    b_rho_u_2(1)  = (RHO_U(cell->neighbor_index(3)) - rho_u0);   // P2 neighbor

                    b_rho_v_2(0)  = (RHO_V(cell->neighbor_index(0)) - rho_v0);   // P1 neighbor 
                    b_rho_v_2(1)  = (RHO_V(cell->neighbor_index(3)) - rho_v0);   // P2 neighbor

                    b_e_2(0)  = (E(cell->neighbor_index(0)) - e0);   // P1 neighbor 
                    b_e_2(1)  = (E(cell->neighbor_index(3)) - e0);   // P2 neighbor

                }
                
                LU_R21[c].solve(b_rho2, rho_coeff_21);
                LU_R21[c].solve(b_rho_u_2, rho_u_coeff_21);
                LU_R21[c].solve(b_rho_v_2, rho_v_coeff_21);
                LU_R21[c].solve(b_e_2, e_coeff_21);

                
                // =====================================================================
                // r = 2 stencil 2 (boundary)
                // =====================================================================
                
                if (cell->neighbor_index(1) == -1) { // EAST boundary 
                
                    b_rho2(0) = (RHO(cell->neighbor_index(3)) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(3)) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(3)) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(cell->neighbor_index(3)) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 
                }
                
                else if (cell->neighbor_index(3) == -1) { // NORTH boundary 
                
                    b_rho2(0) = (RHO(cell->neighbor_index(1)) - rho0);   // P3 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(1)) - rho_u0);   // P3 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(1)) - rho_v0);   // P3 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(cell->neighbor_index(1)) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 

                }
                
                else {
                    
                    b_rho2(0) = (RHO(cell->neighbor_index(3)) - rho0);   // P2 neighbor
                    b_rho2(1) = (RHO(cell->neighbor_index(1)) - rho0);   // P3 neighbor
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(3)) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = (RHO_U(cell->neighbor_index(1)) - rho_u0);   // P3 neighbor

                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(3)) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = (RHO_V(cell->neighbor_index(1)) - rho_v0);   // P3 neighbor

                    b_e_2(0) = (E(cell->neighbor_index(3)) - e0);   // P2 neighbor
                    b_e_2(1) = (E(cell->neighbor_index(1)) - e0);   // P3 neighbor

                }
                
                LU_R22[c].solve(b_rho2, rho_coeff_22);
                LU_R22[c].solve(b_rho_u_2, rho_u_coeff_22);
                LU_R22[c].solve(b_rho_v_2, rho_v_coeff_22);
                LU_R22[c].solve(b_e_2, e_coeff_22);
                

                // =====================================================================
                // r = 2 stencil 3 (boundary)
                // =====================================================================
                
                if (cell->neighbor_index(1) == -1) { // EAST boundary 
                
                    b_rho2(0) = (RHO(cell->neighbor_index(2)) - rho0);   // P2 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(2)) - rho_u0);   // P2 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(2)) - rho_v0);   // P2 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(cell->neighbor_index(2)) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 

                }
                
                else if (cell->neighbor_index(2) == -1) { // SOUTH boundary 
                
                    b_rho2(0) = (RHO(cell->neighbor_index(1)) - rho0);   // P3 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(1)) - rho_u0);   // P3 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(1)) - rho_v0);   // P3 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(cell->neighbor_index(1)) - e0);   // P2 neighbor
                    b_e_2(1) = 0.0; 

                }
                
                else {
                
                    b_rho2(0) = (RHO(cell->neighbor_index(1)) - rho0);   // P3 neighbor
                    b_rho2(1) = (RHO(cell->neighbor_index(2)) - rho0);   // P4 neighbor                                            // P4 neighbor
                     
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(1)) - rho_u0);   // P3 neighbor
                    b_rho_u_2(1) = (RHO_U(cell->neighbor_index(2)) - rho_u0);   // P4 neighbor                                            // P4 neighbor
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(1)) - rho_v0);   // P3 neighbor
                    b_rho_v_2(1) = (RHO_V(cell->neighbor_index(2)) - rho_v0);   // P4 neighbor                                            // P4 neighbor
                    
                    b_e_2(0) = (E(cell->neighbor_index(1)) - e0);   // P3 neighbor
                    b_e_2(1) = (E(cell->neighbor_index(2)) - e0);   // P4 neighbor                                            // P4 neighbor

                }
                
                LU_R23[c].solve(b_rho2, rho_coeff_23);
                LU_R23[c].solve(b_rho_u_2, rho_u_coeff_23);
                LU_R23[c].solve(b_rho_v_2, rho_v_coeff_23);
                LU_R23[c].solve(b_e_2, e_coeff_23);
                
                // =====================================================================
                // r = 2 stencil 4 (boundary)
                // =====================================================================
                
                if (cell->neighbor_index(0) == -1) { // WEST boundary 
                
                    b_rho2(0) = (RHO(cell->neighbor_index(2)) - rho0);   // P4 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(2)) - rho_u0);   // P4 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(2)) - rho_v0);   // P4 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(cell->neighbor_index(2)) - e0);   // P4 neighbor
                    b_e_2(1) = 0.0; 

                    
                }
                
                else if (cell->neighbor_index(2) == -1) { // SOUTH boundary 
                
                    b_rho2(0) = (RHO(cell->neighbor_index(0)) - rho0);   // P1 neighbor
                    b_rho2(1) = 0.0; 
                    
                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(0)) - rho_u0);   // P1 neighbor
                    b_rho_u_2(1) = 0.0; 
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(0)) - rho_v0);   // P1 neighbor
                    b_rho_v_2(1) = 0.0; 
                    
                    b_e_2(0) = (E(cell->neighbor_index(0)) - e0);   // P1 1eighbor
                    b_e_2(1) = 0.0; 


                }
                
                else {
                    b_rho2(0) = (RHO(cell->neighbor_index(2)) - rho0);   // P4 neighbor
                    b_rho2(1) = (RHO(cell->neighbor_index(0)) - rho0);   // P1 neighbor

                    b_rho_u_2(0) = (RHO_U(cell->neighbor_index(2)) - rho_u0);   // P4 neighbor
                    b_rho_u_2(1) = (RHO_U(cell->neighbor_index(0)) - rho_u0);   // P1 neighbor
                    
                    b_rho_v_2(0) = (RHO_V(cell->neighbor_index(2)) - rho_v0);   // P4 neighbor
                    b_rho_v_2(1) = (RHO_V(cell->neighbor_index(0)) - rho_v0);   // P1 neighbor
                    
                    b_e_2(0) = (E(cell->neighbor_index(2)) - e0);   // P4 neighbor
                    b_e_2(1) = (E(cell->neighbor_index(0)) - e0);   // P1 neighbor
                    
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
                                b_3(index)   = (RHO_U(cell->neighbor_index(0)) - rho_u0); // P1 (rho_u)
                                b_3(index+3) = (RHO_V(cell->neighbor_index(0)) - rho_v0); // P1 (rho_v)
                                index++; 
                            }
                            
                            if (cell->neighbor_index(3) != -1 ) {
                                b_3(index)   = (RHO_U(cell->neighbor_index(3)) - rho_u0);  // P2 (u)
                                b_3(index+3) = (RHO_V(cell->neighbor_index(3)) - rho_v0);  // P2 (v)
                                index++; 
                            }
                            
                            if (cell->neighbor_index(1) != -1 ) {
                                b_3(index)   = (RHO_U(cell->neighbor_index(1)) - rho_u0);  // P3 (u)
                                b_3(index+3) = (RHO_V(cell->neighbor_index(1)) - rho_v0);  // P3 (v)
                                index++; 
                            }
                            
                            if (cell->neighbor_index(2) != -1 ) {
                                b_3(index)   = (RHO_U(cell->neighbor_index(2)) - rho_u0);  // P4 (u)
                                b_3(index+3) = (RHO_V(cell->neighbor_index(2)) - rho_v0);  // P4 (v)
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
        
        for (unsigned int j = 0; j < no_stencils; j++ ) {
            
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
        
        coeffs_RHO[c](0) = rho0; 
        
        coeffs_RHO[c](1) = (w_RHO(0)/gamma(0))*(rho_coeff_3(0) - 
                                                gamma(1)*rho_coeff_21(0) -
                                                gamma(2)*rho_coeff_22(0) -
                                                gamma(3)*rho_coeff_23(0) -
                                                gamma(4)*rho_coeff_24(0) ) + 
                                                
                                                w_RHO(1)*rho_coeff_21(0) +
                                                w_RHO(2)*rho_coeff_22(0) +
                                                w_RHO(3)*rho_coeff_23(0) +
                                                w_RHO(4)*rho_coeff_24(0) ; // u_x
        
        coeffs_RHO[c](2) = (w_RHO(0)/gamma(0))*(rho_coeff_3(1) - 
                                                gamma(1)*rho_coeff_21(1) -
                                                gamma(2)*rho_coeff_22(1) -
                                                gamma(3)*rho_coeff_23(1) -
                                                gamma(4)*rho_coeff_24(1) ) + 
                                                
                                                w_RHO(1)*rho_coeff_21(1) +
                                                w_RHO(2)*rho_coeff_22(1) +
                                                w_RHO(3)*rho_coeff_23(1) +
                                                w_RHO(4)*rho_coeff_24(1) ; // u_y

        coeffs_RHO[c](3) = (w_RHO(0)/gamma(0))*(rho_coeff_3(2)); // u_xx
                                        
                                                
        coeffs_RHO[c](4) = (w_RHO(0)/gamma(0))*(rho_coeff_3(3)); // u_yy
        
        coeffs_RHO[c](5) = (w_RHO(0)/gamma(0))*(rho_coeff_3(4)); // u_xy
        
        
        // x-momentum
        
        coeffs_RHO_U[c](0) = rho_u0; 
        
        coeffs_RHO_U[c](1) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(0) - 
                                                gamma(1)*rho_u_coeff_21(0) -
                                                gamma(2)*rho_u_coeff_22(0) -
                                                gamma(3)*rho_u_coeff_23(0) -
                                                gamma(4)*rho_u_coeff_24(0) ) + 
                                                
                                                w_RHO_U(1)*rho_u_coeff_21(0) +
                                                w_RHO_U(2)*rho_u_coeff_22(0) +
                                                w_RHO_U(3)*rho_u_coeff_23(0) +
                                                w_RHO_U(4)*rho_u_coeff_24(0) ; // u_x
        
        coeffs_RHO_U[c](2) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(1) - 
                                                gamma(1)*rho_u_coeff_21(1) -
                                                gamma(2)*rho_u_coeff_22(1) -
                                                gamma(3)*rho_u_coeff_23(1) -
                                                gamma(4)*rho_u_coeff_24(1) ) + 
                                                
                                                w_RHO_U(1)*rho_u_coeff_21(1) +
                                                w_RHO_U(2)*rho_u_coeff_22(1) +
                                                w_RHO_U(3)*rho_u_coeff_23(1) +
                                                w_RHO_U(4)*rho_u_coeff_24(1) ; // u_y

        coeffs_RHO_U[c](3) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(2)); // u_xx
                                        
                                                
        coeffs_RHO_U[c](4) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(3)); // u_yy
        
        coeffs_RHO_U[c](5) = (w_RHO_U(0)/gamma(0))*(rho_u_coeff_3(4)); // u_xy
        
        
        // y-momentum
        
        coeffs_RHO_V[c](0) = rho_v0; 
        
        coeffs_RHO_V[c](1) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(0) - 
                                                gamma(1)*rho_v_coeff_21(0) -
                                                gamma(2)*rho_v_coeff_22(0) -
                                                gamma(3)*rho_v_coeff_23(0) -
                                                gamma(4)*rho_v_coeff_24(0) ) + 
                                                
                                                w_RHO_V(1)*rho_v_coeff_21(0) +
                                                w_RHO_V(2)*rho_v_coeff_22(0) +
                                                w_RHO_V(3)*rho_v_coeff_23(0) +
                                                w_RHO_V(4)*rho_v_coeff_24(0) ; // v_x
        
        coeffs_RHO_V[c](2) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(1) - 
                                                gamma(1)*rho_v_coeff_21(1) -
                                                gamma(2)*rho_v_coeff_22(1) -
                                                gamma(3)*rho_v_coeff_23(1) -
                                                gamma(4)*rho_v_coeff_24(1) ) + 
                                                
                                                w_RHO_V(1)*rho_v_coeff_21(1) +
                                                w_RHO_V(2)*rho_v_coeff_22(1) +
                                                w_RHO_V(3)*rho_v_coeff_23(1) +
                                                w_RHO_V(4)*rho_v_coeff_24(1) ; // v_y

        coeffs_RHO_V[c](3) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(2)); // v_xx
                                        
                                                
        coeffs_RHO_V[c](4) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(3)); // v_yy
        
        coeffs_RHO_V[c](5) = (w_RHO_V(0)/gamma(0))*(rho_v_coeff_3(4)); // v_xy
        
        // Total energy  
        
        coeffs_E[c](0) = e0; 
        
        coeffs_E[c](1) = (w_E(0)/gamma(0))*(e_coeff_3(0) - 
                                                gamma(1)*e_coeff_21(0) -
                                                gamma(2)*e_coeff_22(0) -
                                                gamma(3)*e_coeff_23(0) -
                                                gamma(4)*e_coeff_24(0) ) + 
                                                
                                                w_E(1)*e_coeff_21(0) +
                                                w_E(2)*e_coeff_22(0) +
                                                w_E(3)*e_coeff_23(0) +
                                                w_E(4)*e_coeff_24(0) ; // v_x
        
        coeffs_E[c](2) = (w_E(0)/gamma(0))*(e_coeff_3(1) - 
                                                gamma(1)*e_coeff_21(1) -
                                                gamma(2)*e_coeff_22(1) -
                                                gamma(3)*e_coeff_23(1) -
                                                gamma(4)*e_coeff_24(1) ) + 
                                                
                                                w_E(1)*e_coeff_21(1) +
                                                w_E(2)*e_coeff_22(1) +
                                                w_E(3)*e_coeff_23(1) +
                                                w_E(4)*e_coeff_24(1) ; // v_y

        coeffs_E[c](3) = (w_E(0)/gamma(0))*(e_coeff_3(2)); // v_xx
                                        
                                                
        coeffs_E[c](4) = (w_E(0)/gamma(0))*(e_coeff_3(3)); // v_yy
        
        coeffs_E[c](5) = (w_E(0)/gamma(0))*(e_coeff_3(4)); // v_xy

        
    } // End of cell loop 
    
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

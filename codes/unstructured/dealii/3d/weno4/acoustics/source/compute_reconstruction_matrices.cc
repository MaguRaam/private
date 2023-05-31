#include "../include/Weno432.h"

// Precompute reconstruction matrices 

void Weno4_3D::precompute_matrices() {

    TimerOutput::Scope t(computing_timer, "Computing the stencil matrices");

//    auto start = std::chrono::system_clock::now();

    pcout << "Computing the stencil matrices" << std::endl;

    unsigned int N_gp = 2;  // No. of quadrature points
	unsigned int n_quad_points = N_gp*N_gp*N_gp;
    unsigned int n_faces = GeometryInfo<3>::faces_per_cell;  // No. of quadrature points
	
    Point<3> q_point;
    Point<3> C; 
    Point<3> P; 
	unsigned int  n_boundary_quad_points = 2, n_boundary_quad_points_edge = 1;    
	
    double V_neighbor, dx, dy, dz; 

    DoFHandler<3>::active_cell_iterator cell, neighbor;

    FullMatrix<double> A_R4; // Fourth Order Stencil (Least Squares Part)
    FullMatrix<double> C_R4; // Fourth Order Stencil (Constraint Part)
    
    FullMatrix<double> A_R3; // Least Squares Matrix for r=3 stencil 
    FullMatrix<double> C_R3; // Constraint Matrix for r=3 stencil 
    
    // One-sided stencils 

	std::vector < FullMatrix<double> > C_R3_d(6);		// Directional stencil
	std::vector < FullMatrix<double> > A_R3_d(6);

	std::vector < FullMatrix<double> > A_R2_d(8);

	std::vector<bool> n_n_face(6);
	std::vector<unsigned int> n_n_face_index(6);  

    unsigned int ROWS, index, local_index, fff, g_i;
	
	double x0, y0, z0, h, h2, h3, j_w; 
    
	for (unsigned int c = 0; c < n_relevant_cells; ++c) {

		cell = local_index_to_iterator[c];
		g_i = local_to_global_index_map[c];

//		pcout<<"g_i : "<<g_i<<"\tc: "<<c<<std::endl;
      	for (unsigned int f = 0; f < n_faces; f++) n_n_face[f] = false;
		
		x0 = WENO_poly_consts[c](0); 
        y0 = WENO_poly_consts[c](1);
        z0 = WENO_poly_consts[c](2);
		h = Cell[c].h();
		h2 = h*h; h3 = h2*h;

      	for (unsigned int f = 0; f < n_faces; f++){	
			if (cell_neighbor_neighbor_index[c][f].size() > 0) {
				n_n_face[f] = true;
				n_n_face_index[f] = global_to_local_index_map[cell_neighbor_neighbor_index[c][f][0] ];
			}
			else 
				n_n_face_index[f] = numbers::invalid_unsigned_int;
		}
       
        if ( !(cell->at_boundary()) ) {
			
			// =====================================================================
            // r = 4 stencil 
            // =====================================================================
			C_R4.reinit(6, 19); C_R4 = 0.0;
            	
           	x0 = WENO_poly_consts[c](0); 
           	y0 = WENO_poly_consts[c](1);
           	z0 = WENO_poly_consts[c](2);
            	
           	// Fill C_R4 (Constraint r = 4 matrix)	

           	for (unsigned int f = 0; f < n_faces; f++) {			
            	
				neighbor = cell->neighbor(f);
			    neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
	           	V_neighbor = Cell[local_index].measure();
            	
            	for (unsigned int i = 0; i < n_quad_points; i++) {
            	    q_point = Cell[local_index].cell_quadrature_point(i);
					dx = (q_point(0) - x0)/h;
					dy = (q_point(1) - y0)/h;
					dz = (q_point(2) - z0)/h;

					j_w = Cell[local_index].jxw(i);
            	    C_R4(f,0) += (1./V_neighbor)*j_w*dx;
            	    C_R4(f,1) += (1./V_neighbor)*j_w*dy;
            	    C_R4(f,2) += (1./V_neighbor)*j_w*dz;
            	    C_R4(f,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
            	    C_R4(f,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
            	    C_R4(f,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
            	    C_R4(f,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
            	    C_R4(f,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
            	    C_R4(f,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
            	    C_R4(f,9) += (1./V_neighbor)*j_w*(dx*dx*dx - WENO_poly_consts[c](9));
            	    C_R4(f,10) += (1./V_neighbor)*j_w*(dy*dy*dy - WENO_poly_consts[c](10));
            	    C_R4(f,11) += (1./V_neighbor)*j_w*(dz*dz*dz - WENO_poly_consts[c](11));
            	    C_R4(f,12) += (1./V_neighbor)*j_w*(dx*dx*dy - WENO_poly_consts[c](12));
            	    C_R4(f,13) += (1./V_neighbor)*j_w*(dx*dx*dz - WENO_poly_consts[c](13));
            	    C_R4(f,14) += (1./V_neighbor)*j_w*(dx*dy*dy - WENO_poly_consts[c](14));
            	    C_R4(f,15) += (1./V_neighbor)*j_w*(dz*dy*dy - WENO_poly_consts[c](15));
            	    C_R4(f,16) += (1./V_neighbor)*j_w*(dx*dz*dz - WENO_poly_consts[c](16));
            	    C_R4(f,17) += (1./V_neighbor)*j_w*(dy*dz*dz - WENO_poly_consts[c](17));
            	    C_R4(f,18) += (1./V_neighbor)*j_w*(dx*dy*dz - WENO_poly_consts[c](18));
            	}

			}            
            // Least Squares Part  
            
			// First get the number of rows in the least squares matrix part 
			
			// vertex neighbors 

			ROWS = cell_diagonal_neighbor_index[c].size(); 
			
			// neighbor of neighbors 

	      	for (unsigned int f = 0; f < n_faces; f++)
				if (n_n_face[f])
					ROWS++; 
				
			index = 0; 
            
            A_R4.reinit(ROWS, 19); A_R4 = 0.0;
//			if(c == 8570) pcout<<"A_R4 size: "<<ROWS<<std::endl; 						
			// First fill the vertex neighbors 
			
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			
	
				local_index = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];
    	        V_neighbor = Cell[local_index].measure();
				
            	for (unsigned int i = 0; i < n_quad_points; i++) {
            	    q_point = Cell[local_index].cell_quadrature_point(i);
					j_w = Cell[local_index].jxw(i);
					dx = (q_point(0) - x0)/h;
					dy = (q_point(1) - y0)/h;
					dz = (q_point(2) - z0)/h;

            	    A_R4(index,0) += (1./V_neighbor)*j_w*dx;
            	    A_R4(index,1) += (1./V_neighbor)*j_w*dy;
            	    A_R4(index,2) += (1./V_neighbor)*j_w*dz;
            	    A_R4(index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
            	    A_R4(index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
            	    A_R4(index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
            	    A_R4(index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
            	    A_R4(index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
            	    A_R4(index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
            	    A_R4(index,9) += (1./V_neighbor)*j_w*(dx*dx*dx - WENO_poly_consts[c](9));
            	    A_R4(index,10) += (1./V_neighbor)*j_w*(dy*dy*dy - WENO_poly_consts[c](10));
            	    A_R4(index,11) += (1./V_neighbor)*j_w*(dz*dz*dz - WENO_poly_consts[c](11));
            	    A_R4(index,12) += (1./V_neighbor)*j_w*(dx*dx*dy - WENO_poly_consts[c](12));
            	    A_R4(index,13) += (1./V_neighbor)*j_w*(dx*dx*dz - WENO_poly_consts[c](13));
            	    A_R4(index,14) += (1./V_neighbor)*j_w*(dx*dy*dy - WENO_poly_consts[c](14));
            	    A_R4(index,15) += (1./V_neighbor)*j_w*(dz*dy*dy - WENO_poly_consts[c](15));
            	    A_R4(index,16) += (1./V_neighbor)*j_w*(dx*dz*dz - WENO_poly_consts[c](16));
            	    A_R4(index,17) += (1./V_neighbor)*j_w*(dy*dz*dz - WENO_poly_consts[c](17));
            	    A_R4(index,18) += (1./V_neighbor)*j_w*(dx*dy*dz - WENO_poly_consts[c](18));
            	}
                
                index++; 
			}

	      	for (unsigned int f = 0; f < n_faces; f++) {
				if (n_n_face[f]){

					local_index = n_n_face_index[f];
		            V_neighbor = Cell[local_index].measure();

	            	for (unsigned int i = 0; i < n_quad_points; i++) {
	            	    q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						dx = (q_point(0) - x0)/h;
						dy = (q_point(1) - y0)/h;
						dz = (q_point(2) - z0)/h;

	            	    A_R4(index,0) += (1./V_neighbor)*j_w*dx;
	            	    A_R4(index,1) += (1./V_neighbor)*j_w*dy;
	            	    A_R4(index,2) += (1./V_neighbor)*j_w*dz;
	            	    A_R4(index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
	            	    A_R4(index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
	            	    A_R4(index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
	            	    A_R4(index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
	            	    A_R4(index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
	            	    A_R4(index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
	            	    A_R4(index,9) += (1./V_neighbor)*j_w*(dx*dx*dx - WENO_poly_consts[c](9));
	            	    A_R4(index,10) += (1./V_neighbor)*j_w*(dy*dy*dy - WENO_poly_consts[c](10));
	            	    A_R4(index,11) += (1./V_neighbor)*j_w*(dz*dz*dz - WENO_poly_consts[c](11));
	            	    A_R4(index,12) += (1./V_neighbor)*j_w*(dx*dx*dy - WENO_poly_consts[c](12));
	            	    A_R4(index,13) += (1./V_neighbor)*j_w*(dx*dx*dz - WENO_poly_consts[c](13));
	            	    A_R4(index,14) += (1./V_neighbor)*j_w*(dx*dy*dy - WENO_poly_consts[c](14));
	            	    A_R4(index,15) += (1./V_neighbor)*j_w*(dz*dy*dy - WENO_poly_consts[c](15));
	            	    A_R4(index,16) += (1./V_neighbor)*j_w*(dx*dz*dz - WENO_poly_consts[c](16));
	            	    A_R4(index,17) += (1./V_neighbor)*j_w*(dy*dz*dz - WENO_poly_consts[c](17));
	            	    A_R4(index,18) += (1./V_neighbor)*j_w*(dx*dy*dz - WENO_poly_consts[c](18));
	            	}
                index++; 
				}
            }
            
            CLS_R4[c].initialize(A_R4, C_R4); 
            
            // =====================================================================
            // r = 3 stencil (Centered Stencil)
            // =====================================================================
			C_R3.reinit(6, 9); C_R3 = 0.0;
			
			// Fill C_R3 (Constraint r = 3 matrix)

           	for (unsigned int f = 0; f < n_faces; f++) {			
            	
				neighbor = cell->neighbor(f);
			    neighbor->get_dof_indices(local_neighbor_dof_indices);
				local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
	           	V_neighbor = Cell[local_index].measure();
            	
            	for (unsigned int i = 0; i < n_quad_points; i++) {
            	    q_point = Cell[local_index].cell_quadrature_point(i);
					j_w = Cell[local_index].jxw(i);
					dx = (q_point(0) - x0)/h;
					dy = (q_point(1) - y0)/h;
					dz = (q_point(2) - z0)/h;
            	    C_R3(f,0) += (1./V_neighbor)*j_w*dx;
            	    C_R3(f,1) += (1./V_neighbor)*j_w*dy;
            	    C_R3(f,2) += (1./V_neighbor)*j_w*dz;
            	    C_R3(f,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
            	    C_R3(f,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
            	    C_R3(f,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
            	    C_R3(f,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
            	    C_R3(f,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
            	    C_R3(f,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
            	}

			}
/*
			for(unsigned int i = 0; i < 6; ++i) {
				for(unsigned int j = 0; j < 6; ++j)
					std::cout<<C_R3(i,j)<<"\t";
				std::cout<<std::endl;
			}
*/
			// Fill A_R3 (Least squares r = 3 matrix)

			ROWS = cell_diagonal_neighbor_index[c].size();
//			if(c == 8570) pcout<<"A_R3 size: "<<ROWS<<std::endl; 				
			A_R3.reinit(ROWS, 9); A_R3 = 0.0; index = 0; 
			
			// First fill the vertex neighbors 

			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			
	
				local_index = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];
    	        V_neighbor = Cell[local_index].measure();
				
            	for (unsigned int i = 0; i < n_quad_points; i++) {
            	    q_point = Cell[local_index].cell_quadrature_point(i);
					j_w = Cell[local_index].jxw(i);
					dx = (q_point(0) - x0)/h;
					dy = (q_point(1) - y0)/h;
					dz = (q_point(2) - z0)/h;
            	    A_R3(index,0) += (1./V_neighbor)*j_w*dx;
            	    A_R3(index,1) += (1./V_neighbor)*j_w*dy;
            	    A_R3(index,2) += (1./V_neighbor)*j_w*dz;
            	    A_R3(index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
            	    A_R3(index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
            	    A_R3(index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
            	    A_R3(index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
            	    A_R3(index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
            	    A_R3(index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
            	}
                
                index++; 
			}
			
			CLS_R3[c].initialize(A_R3, C_R3);

            // =====================================================================
            // r = 3 directional stencil 
            // =====================================================================

	      	for (unsigned int f = 0; f < n_faces; ++f) {
				if (n_n_face[f]) {
					if (cell_neighbor_index[c][f].size() >= 3) {
						is_admissible_R3_d[c][f] = true;
					}
					else {
						is_admissible_R3_d[c][f] = false;
					}
				}
	            else {
					is_admissible_R3_d[c][f] = false;
				}

            	if (is_admissible_R3_d[c][f]) {
            
                	// Constraint Matrix 
                	C_R3_d[f].reinit(5,9); 
            
                	C_R3_d[f] = 0.0; index = 0; 

					if(f%2 == 0)
						fff = f+1;            	
					else 
						fff = f-1;

	           		for (unsigned int ff = 0; ff < n_faces; ff++) {

						if(fff != ff) {
							neighbor = cell->neighbor(ff);
						    neighbor->get_dof_indices(local_neighbor_dof_indices);
							local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
	    			       	V_neighbor = Cell[local_index].measure();
	
        		    		for (unsigned int i = 0; i < n_quad_points; i++) {
        		    		    q_point = Cell[local_index].cell_quadrature_point(i);
								j_w = Cell[local_index].jxw(i);
								dx = (q_point(0) - x0)/h;
								dy = (q_point(1) - y0)/h;
								dz = (q_point(2) - z0)/h;
        		    		    C_R3_d[f](index,0) += (1./V_neighbor)*j_w*dx;
        		    		    C_R3_d[f](index,1) += (1./V_neighbor)*j_w*dy;
        		    		    C_R3_d[f](index,2) += (1./V_neighbor)*j_w*dz;
        		    		    C_R3_d[f](index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
        		    		    C_R3_d[f](index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
        		    		    C_R3_d[f](index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
        		    		    C_R3_d[f](index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
        		    		    C_R3_d[f](index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
        		    		    C_R3_d[f](index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
        		    		}
							index++;
						}

					}

					Assert (index == 5, ExcDimensionMismatch(index, 5));                
                          

    	            
	                // Least Squares Matrix 
                
    	            index = 0; 
    	            ROWS = cell_neighbor_index[c][f].size() + 1; 
	
    	            A_R3_d[f].reinit(ROWS, 9); A_R3_d[f] = 0.0; 
//					if(c == 8570) pcout<<"A_R3_d["<<f<<"] size: "<<ROWS<<std::endl; 	
					// vertex neighbor of cell at face f
	
					for (unsigned int d = 0; d < ROWS - 1; ++d) {			
	
						local_index = global_to_local_index_map[cell_neighbor_index[c][f][d] ];
		    	        V_neighbor = Cell[local_index].measure();
				
        		    	for (unsigned int i = 0; i < n_quad_points; i++) {
        		    	    q_point = Cell[local_index].cell_quadrature_point(i);
							j_w = Cell[local_index].jxw(i);
							dx = (q_point(0) - x0)/h;
							dy = (q_point(1) - y0)/h;
							dz = (q_point(2) - z0)/h;
        		    	    A_R3_d[f](index,0) += (1./V_neighbor)*j_w*dx;
        		    	    A_R3_d[f](index,1) += (1./V_neighbor)*j_w*dy;
        		    	    A_R3_d[f](index,2) += (1./V_neighbor)*j_w*dz;
        		    	    A_R3_d[f](index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
        		    	    A_R3_d[f](index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
        		    	    A_R3_d[f](index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
        		    	    A_R3_d[f](index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
        		    	    A_R3_d[f](index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
        		    	    A_R3_d[f](index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
        		    	}
                
						index++; 
					}

					local_index = n_n_face_index[f];
    	       	 	V_neighbor = Cell[local_index].measure();
                
    	    		for (unsigned int i = 0; i < n_quad_points; i++) {
    	    		    q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						dx = (q_point(0) - x0)/h;
						dy = (q_point(1) - y0)/h;
						dz = (q_point(2) - z0)/h;
    	    		    A_R3_d[f](index,0) += (1./V_neighbor)*j_w*dx;
    	    		    A_R3_d[f](index,1) += (1./V_neighbor)*j_w*dy;
    	    		    A_R3_d[f](index,2) += (1./V_neighbor)*j_w*dz;
    	    		    A_R3_d[f](index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
    	    		    A_R3_d[f](index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
    	    		    A_R3_d[f](index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
    	    		    A_R3_d[f](index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
    	    		    A_R3_d[f](index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
    	    		    A_R3_d[f](index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
    	    		}
								
					CLS_R3_d[c][f].initialize(A_R3_d[f], C_R3_d[f]);
            	}

            }

/*
 *         6-------7        6-------7
 *        /|       |       /       /|
 *       / |       |      /   f3  / |
 *      /  |  f5   |     /       /  |
 *     4   |       |    4-------5   |
 *     |f0 2-------3    |       |f1 3			z	y
 *     |  /       /     |       |  /			|  /    
 *     | /   f2  /      |  f4   | /				| /   
 *     |/       /       |       |/				|/      
 *     0-------1        0-------1               ------x
*/
/*
	    	for (unsigned int v = 0; v < 8; v++)
				A_R2_d[v].reinit(3,3);
			            
			// =====================================================================
            // r = 2 stencil 0 (f0 f2 & f4 )	
            // =====================================================================
            
            A_R2_d[0](0,0) = C_R3(0,0); A_R2_d[0](0,1) = C_R3(0,1); A_R2_d[0](0,2) = C_R3(0,2); // Row 1 f0
            A_R2_d[0](1,0) = C_R3(2,0); A_R2_d[0](1,1) = C_R3(2,1); A_R2_d[0](1,2) = C_R3(2,2); // Row 1 f2
            A_R2_d[0](2,0) = C_R3(4,0); A_R2_d[0](2,1) = C_R3(4,1); A_R2_d[0](2,2) = C_R3(4,2); // Row 1 f4

			std::cout<<"0: "<<c<<"\tC_R3: " <<std::endl;
			for(unsigned int i = 0; i < 3; ++i) 
				std::cout<<C_R3(2,i)<<"\t";
			std::cout<<std::endl;
				std::cout<<"0: "<<c<<std::endl;
			for(unsigned int i = 0; i < 3; ++i) {
				for(unsigned int j = 0; j < 3; ++j)
					std::cout<<A_R2_d[0](i,j)<<"\t";
				std::cout<<std::endl;
			}

            LU_R2_d[c][0].initialize(A_R2_d[0]);

//			if(A_R2_d[0](0,0) == 0.0 && A_R2_d[0](0,1) == 0.0 && A_R2_d[0](0,2) == 0.0 )
//				std::cout<<"0: "<<c<<std::endl;

			// =====================================================================
            // r = 2 stencil 1 (f1 f2 & f4 )	
            // =====================================================================
            
            A_R2_d[1](0,0) = C_R3(1,0); A_R2_d[1](0,1) = C_R3(1,1); A_R2_d[1](0,2) = C_R3(1,2); // Row 1 f1
            A_R2_d[1](1,0) = C_R3(2,0); A_R2_d[1](1,1) = C_R3(2,1); A_R2_d[1](1,2) = C_R3(2,2); // Row 1 f2
            A_R2_d[1](2,0) = C_R3(4,0); A_R2_d[1](2,1) = C_R3(4,1); A_R2_d[1](2,2) = C_R3(4,2); // Row 1 f4
            LU_R2_d[c][1].initialize(A_R2_d[1]);
//			std::cout<<"1: "<<std::endl;
//			if(A_R2_d[1](0,0) == 0.0 && A_R2_d[1](0,1) == 0.0 && A_R2_d[1](0,2) == 0.0 )
//				std::cout<<"1: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 2 (f0 f2 & f5 )	
            // =====================================================================
            
            A_R2_d[2](0,0) = C_R3(0,0); A_R2_d[2](0,1) = C_R3(0,1); A_R2_d[2](0,2) = C_R3(0,2); // Row 1 f0
            A_R2_d[2](1,0) = C_R3(2,0); A_R2_d[2](1,1) = C_R3(2,1); A_R2_d[2](1,2) = C_R3(2,2); // Row 1 f2
            A_R2_d[2](2,0) = C_R3(5,0); A_R2_d[2](2,1) = C_R3(5,1); A_R2_d[2](2,2) = C_R3(5,2); // Row 1 f5
            LU_R2_d[c][2].initialize(A_R2_d[2]);
//			std::cout<<"2: "<<std::endl;
//			if(A_R2_d[2](0,0) == 0.0 && A_R2_d[2](0,1) == 0.0 && A_R2_d[2](0,2) == 0.0 )
//				std::cout<<"2: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 3 (f1 f2 & f5 )	
            // =====================================================================
            
            A_R2_d[3](0,0) = C_R3(1,0); A_R2_d[3](0,1) = C_R3(1,1); A_R2_d[3](0,2) = C_R3(1,2); // Row 1 f1
            A_R2_d[3](1,0) = C_R3(2,0); A_R2_d[3](1,1) = C_R3(2,1); A_R2_d[3](1,2) = C_R3(2,2); // Row 1 f2
            A_R2_d[3](2,0) = C_R3(5,0); A_R2_d[3](2,1) = C_R3(5,1); A_R2_d[3](2,2) = C_R3(5,2); // Row 1 f5
            LU_R2_d[c][3].initialize(A_R2_d[3]);
//			std::cout<<"3: "<<std::endl;
//			if(A_R2_d[3](0,0) == 0.0 && A_R2_d[3](0,1) == 0.0 && A_R2_d[3](0,2) == 0.0 )
//				std::cout<<"3: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 4 (f0 f3 & f4 )	
            // =====================================================================
            
            A_R2_d[4](0,0) = C_R3(0,0); A_R2_d[4](0,1) = C_R3(0,1); A_R2_d[4](0,2) = C_R3(0,2); // Row 1 f0
            A_R2_d[4](1,0) = C_R3(3,0); A_R2_d[4](1,1) = C_R3(3,1); A_R2_d[4](1,2) = C_R3(3,2); // Row 1 f3
            A_R2_d[4](2,0) = C_R3(4,0); A_R2_d[4](2,1) = C_R3(4,1); A_R2_d[4](2,2) = C_R3(4,2); // Row 1 f4
            LU_R2_d[c][4].initialize(A_R2_d[4]);
//			std::cout<<"4: "<<std::endl;
//			if(A_R2_d[4](0,0) == 0.0 && A_R2_d[4](0,1) == 0.0 && A_R2_d[4](0,2) == 0.0 )
//				std::cout<<"4: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 5 (f1 f3 & f4 )	
            // =====================================================================
            
            A_R2_d[5](0,0) = C_R3(1,0); A_R2_d[5](0,1) = C_R3(1,1); A_R2_d[5](0,2) = C_R3(1,2); // Row 1 f1
            A_R2_d[5](1,0) = C_R3(3,0); A_R2_d[5](1,1) = C_R3(3,1); A_R2_d[5](1,2) = C_R3(3,2); // Row 1 f3
            A_R2_d[5](2,0) = C_R3(4,0); A_R2_d[5](2,1) = C_R3(4,1); A_R2_d[5](2,2) = C_R3(4,2); // Row 1 f4
            LU_R2_d[c][5].initialize(A_R2_d[5]);
//			std::cout<<"5: "<<std::endl;
//			if(A_R2_d[5](0,0) == 0.0 && A_R2_d[5](0,1) == 0.0 && A_R2_d[5](0,2) == 0.0 )
//				std::cout<<"5: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 6 (f0 f3 & f5 )	
            // =====================================================================
            
            A_R2_d[6](0,0) = C_R3(0,0); A_R2_d[6](0,1) = C_R3(0,1); A_R2_d[6](0,2) = C_R3(0,2); // Row 1 f0
            A_R2_d[6](1,0) = C_R3(3,0); A_R2_d[6](1,1) = C_R3(3,1); A_R2_d[6](1,2) = C_R3(3,2); // Row 1 f3
            A_R2_d[6](2,0) = C_R3(5,0); A_R2_d[6](2,1) = C_R3(5,1); A_R2_d[6](2,2) = C_R3(5,2); // Row 1 f5
            LU_R2_d[c][6].initialize(A_R2_d[6]);
//			std::cout<<"6: "<<std::endl;
//			if(A_R2_d[6](0,0) == 0.0 && A_R2_d[6](0,1) == 0.0 && A_R2_d[6](0,2) == 0.0 )
//				std::cout<<"6: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 7 (f1 f3 & f5 )	
            // =====================================================================
            
            A_R2_d[7](0,0) = C_R3(1,0); A_R2_d[7](0,1) = C_R3(1,1); A_R2_d[7](0,2) = C_R3(1,2); // Row 1 f1
            A_R2_d[7](1,0) = C_R3(3,0); A_R2_d[7](1,1) = C_R3(3,1); A_R2_d[7](1,2) = C_R3(3,2); // Row 1 f3
            A_R2_d[7](2,0) = C_R3(5,0); A_R2_d[7](2,1) = C_R3(5,1); A_R2_d[7](2,2) = C_R3(5,2); // Row 1 f5
            LU_R2_d[c][7].initialize(A_R2_d[7]);
//			if(A_R2_d[7](0,0) == 0.0 && A_R2_d[7](0,1) == 0.0 && A_R2_d[7](0,2) == 0.0 )
//				std::cout<<"7: "<<c<<std::endl;
//			std::cout<<"7: "<<std::endl;
*/


		} // End of interior cell loop 
        
		else {

			unsigned int n_boundary_faces = 0, n_no_slip_boundary = 0;
			bool no_slip_boundary = false, transmissive_boundary = false;
		    for (unsigned int f=0; f<n_faces; ++f) {
				if(cell->face(f)->at_boundary()) {
					if(cell->face(f)->boundary_id() == 2) {	
						n_no_slip_boundary += 1;
						no_slip_boundary = true;
					}
					if(cell->face(f)->boundary_id() == 1 || cell->face(f)->boundary_id() == 0) {	
						transmissive_boundary = true;
					}
					n_boundary_faces += 1;
				}
			}

			// =====================================================================
            // r = 4 stencil 
            // =====================================================================

//			if (n_boundary_faces < 3) {

				unsigned int index_least = 0;

				ROWS = cell_diagonal_neighbor_index[c].size(); 
//				if(c == 36) std::cout<<"cls ROWS: "<<ROWS<<"\tindex_least: "<<index_least<<std::endl;						
				// neighbor of neighbors 

		      	for (unsigned int f = 0; f < n_faces; f++)
					if (n_n_face[f])
						ROWS++; 
				ROWS += n_boundary_faces*4;
//				if(c == 36) std::cout<<"cls ROWS: "<<ROWS<<"\tindex_least: "<<index_least<<std::endl;			
    	        A_R4.reinit(ROWS, 19); A_R4 = 0.0;
			
				for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			
	
					local_index = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];
	    	        V_neighbor = Cell[local_index].measure();
				
	            	for (unsigned int i = 0; i < n_quad_points; i++) {
	            	    q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						dx = (q_point(0) - x0)/h;
						dy = (q_point(1) - y0)/h;
						dz = (q_point(2) - z0)/h;
	
	            	    A_R4(index_least,0) += (1./V_neighbor)*j_w*dx;
	            	    A_R4(index_least,1) += (1./V_neighbor)*j_w*dy;
	            	    A_R4(index_least,2) += (1./V_neighbor)*j_w*dz;
	            	    A_R4(index_least,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
	            	    A_R4(index_least,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
	            	    A_R4(index_least,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
	            	    A_R4(index_least,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
	            	    A_R4(index_least,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
	            	    A_R4(index_least,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
    	        	    A_R4(index_least,9) += (1./V_neighbor)*j_w*(dx*dx*dx - WENO_poly_consts[c](9));
    	        	    A_R4(index_least,10) += (1./V_neighbor)*j_w*(dy*dy*dy - WENO_poly_consts[c](10));
    	        	    A_R4(index_least,11) += (1./V_neighbor)*j_w*(dz*dz*dz - WENO_poly_consts[c](11));
    	        	    A_R4(index_least,12) += (1./V_neighbor)*j_w*(dx*dx*dy - WENO_poly_consts[c](12));
    	        	    A_R4(index_least,13) += (1./V_neighbor)*j_w*(dx*dx*dz - WENO_poly_consts[c](13));
    	        	    A_R4(index_least,14) += (1./V_neighbor)*j_w*(dx*dy*dy - WENO_poly_consts[c](14));
    	        	    A_R4(index_least,15) += (1./V_neighbor)*j_w*(dz*dy*dy - WENO_poly_consts[c](15));
    	        	    A_R4(index_least,16) += (1./V_neighbor)*j_w*(dx*dz*dz - WENO_poly_consts[c](16));
    	        	    A_R4(index_least,17) += (1./V_neighbor)*j_w*(dy*dz*dz - WENO_poly_consts[c](17));
    	        	    A_R4(index_least,18) += (1./V_neighbor)*j_w*(dx*dy*dz - WENO_poly_consts[c](18));
    	        	}
    	            
    	            index_least++; 
				}
	
		      	for (unsigned int f = 0; f < n_faces; f++) {
					if (n_n_face[f]){
	
						local_index = n_n_face_index[f];
			            V_neighbor = Cell[local_index].measure();

		            	for (unsigned int i = 0; i < n_quad_points; i++) {
		            	    q_point = Cell[local_index].cell_quadrature_point(i);
							j_w = Cell[local_index].jxw(i);
							dx = (q_point(0) - x0)/h;
							dy = (q_point(1) - y0)/h;
							dz = (q_point(2) - z0)/h;
	
		            	    A_R4(index_least,0) += (1./V_neighbor)*j_w*dx;
		            	    A_R4(index_least,1) += (1./V_neighbor)*j_w*dy;
		            	    A_R4(index_least,2) += (1./V_neighbor)*j_w*dz;
		            	    A_R4(index_least,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
		            	    A_R4(index_least,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
		            	    A_R4(index_least,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
		            	    A_R4(index_least,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
		            	    A_R4(index_least,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
		            	    A_R4(index_least,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
		            	    A_R4(index_least,9) += (1./V_neighbor)*j_w*(dx*dx*dx - WENO_poly_consts[c](9));
		            	    A_R4(index_least,10) += (1./V_neighbor)*j_w*(dy*dy*dy - WENO_poly_consts[c](10));
		            	    A_R4(index_least,11) += (1./V_neighbor)*j_w*(dz*dz*dz - WENO_poly_consts[c](11));
		            	    A_R4(index_least,12) += (1./V_neighbor)*j_w*(dx*dx*dy - WENO_poly_consts[c](12));
		            	    A_R4(index_least,13) += (1./V_neighbor)*j_w*(dx*dx*dz - WENO_poly_consts[c](13));
		            	    A_R4(index_least,14) += (1./V_neighbor)*j_w*(dx*dy*dy - WENO_poly_consts[c](14));
		            	    A_R4(index_least,15) += (1./V_neighbor)*j_w*(dz*dy*dy - WENO_poly_consts[c](15));
		            	    A_R4(index_least,16) += (1./V_neighbor)*j_w*(dx*dz*dz - WENO_poly_consts[c](16));
		            	    A_R4(index_least,17) += (1./V_neighbor)*j_w*(dy*dz*dz - WENO_poly_consts[c](17));
		            	    A_R4(index_least,18) += (1./V_neighbor)*j_w*(dx*dy*dz - WENO_poly_consts[c](18));
		            	}
    	            index_least++; 
					}
    	        }

			ROWS = 6 - n_boundary_faces;
//			ROWS = 6;
            C_R4.reinit(ROWS, 19); 
			C_R4 = 0.0;
			index = 0; 

           	for (unsigned int f = 0; f < n_faces; f++) {			
				if(!cell->face(f)->at_boundary()) {	

					neighbor = cell->neighbor(f);
				    neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
		           	V_neighbor = Cell[local_index].measure();
	
	            	for (unsigned int i = 0; i < n_quad_points; i++) {
	            	    q_point = Cell[local_index].cell_quadrature_point(i);
						dx = (q_point(0) - x0)/h;
						dy = (q_point(1) - y0)/h;
						dz = (q_point(2) - z0)/h;
	
						j_w = Cell[local_index].jxw(i);
	            	    C_R4(index,0) += (1./V_neighbor)*j_w*dx;
	            	    C_R4(index,1) += (1./V_neighbor)*j_w*dy;
	            	    C_R4(index,2) += (1./V_neighbor)*j_w*dz;
	            	    C_R4(index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
	            	    C_R4(index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
	            	    C_R4(index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
	            	    C_R4(index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
	            	    C_R4(index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
	            	    C_R4(index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
	            	    C_R4(index,9) += (1./V_neighbor)*j_w*(dx*dx*dx - WENO_poly_consts[c](9));
	            	    C_R4(index,10) += (1./V_neighbor)*j_w*(dy*dy*dy - WENO_poly_consts[c](10));
	            	    C_R4(index,11) += (1./V_neighbor)*j_w*(dz*dz*dz - WENO_poly_consts[c](11));
	            	    C_R4(index,12) += (1./V_neighbor)*j_w*(dx*dx*dy - WENO_poly_consts[c](12));
	            	    C_R4(index,13) += (1./V_neighbor)*j_w*(dx*dx*dz - WENO_poly_consts[c](13));
	            	    C_R4(index,14) += (1./V_neighbor)*j_w*(dx*dy*dy - WENO_poly_consts[c](14));
	            	    C_R4(index,15) += (1./V_neighbor)*j_w*(dz*dy*dy - WENO_poly_consts[c](15));
	            	    C_R4(index,16) += (1./V_neighbor)*j_w*(dx*dz*dz - WENO_poly_consts[c](16));
	            	    C_R4(index,17) += (1./V_neighbor)*j_w*(dy*dz*dz - WENO_poly_consts[c](17));
	            	    C_R4(index,18) += (1./V_neighbor)*j_w*(dx*dy*dz - WENO_poly_consts[c](18));
	            	}
//							if(c == 5376) pcout<<"index: "<<index<<"\tneigh: "<<local_index<<std::endl;
    	            index++; 		
				}
	
				else {


	                QGauss<3-1> face_quadrature_formula(n_boundary_quad_points_edge);
	                FEFaceValues<3> fv_face_values (mapping, fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);
					fv_face_values.reinit(cell, f);

	                QGauss<3-1> face_quadrature_formula2(n_boundary_quad_points);
	                FEFaceValues<3> fv_face_values2 (mapping, fv, face_quadrature_formula2, update_quadrature_points | update_normal_vectors);
					fv_face_values2.reinit(cell, f);

					double nx, ny, nz;
/*
	                for (unsigned int i = 0; i < 1; i++) {
							nx = fv_face_values.normal_vector(i)[0];
							ny = fv_face_values.normal_vector(i)[1];
							nz = fv_face_values.normal_vector(i)[2];
							q_point = fv_face_values.quadrature_point(i);
							dx = (q_point(0) - x0)/h;
							dy = (q_point(1) - y0)/h;
							dz = (q_point(2) - z0)/h;
        		            C_R4(index,0) = nx;
        		            C_R4(index,1) = ny;
        		            C_R4(index,2) = nz;
    		        	    C_R4(index,3) = 2.0*dx*nx;
    		        	    C_R4(index,4) = 2.0*dy*ny;
    		        	    C_R4(index,5) = 2.0*dz*nz;
    		        	    C_R4(index,6) = dy*nx+dx*ny;
    		        	    C_R4(index,7) = dz*nx+dx*nz;
    		        	    C_R4(index,8) = dz*ny+dy*nz;
    		        	    C_R4(index,9) = 3.0*dx*dx*nx;
    		        	    C_R4(index,10) = 3.0*dy*dy*ny;
    		        	    C_R4(index,11) = 3.0*dz*dz*nz;
    		        	    C_R4(index,12) = 2.0*dx*dy*nx+dx*dx*ny;
    		        	    C_R4(index,13) = 2.0*dx*dz*nx+dx*dx*nz;
    		        	    C_R4(index,14) = dy*dy*nx+2.0*dx*dy*ny;
    		        	    C_R4(index,15) = 2.0*dy*dz*ny+dy*dy*nz;
    		        	    C_R4(index,16) = dz*dz*nx+2.0*dx*dz*nz;
    		        	    C_R4(index,17) = dz*dz*ny+2.0*dy*dz*nz;
    		        	    C_R4(index,18) = dy*dz*nx+dx*dz*ny+dx*dy*nz;

    		            index++; 		

    		        }
*/
	                for (unsigned int i = 0; i < 4; i++) {
							nx = fv_face_values2.normal_vector(i)[0];
							ny = fv_face_values2.normal_vector(i)[1];
							nz = fv_face_values2.normal_vector(i)[2];
							q_point = fv_face_values2.quadrature_point(i);
							dx = (q_point(0) - x0)/h;
							dy = (q_point(1) - y0)/h;
							dz = (q_point(2) - z0)/h;
        		            A_R4(index_least,0) = weight_ls*(nx);
        		            A_R4(index_least,1) = weight_ls*(ny);
        		            A_R4(index_least,2) = weight_ls*(nz);
    		        	    A_R4(index_least,3) = weight_ls*(2.0*dx*nx);
    		        	    A_R4(index_least,4) = weight_ls*(2.0*dy*ny);
    		        	    A_R4(index_least,5) = weight_ls*(2.0*dz*nz);
    		        	    A_R4(index_least,6) = weight_ls*(dy*nx+dx*ny);
    		        	    A_R4(index_least,7) = weight_ls*(dz*nx+dx*nz);
    		        	    A_R4(index_least,8) = weight_ls*(dz*ny+dy*nz);
    		        	    A_R4(index_least,9) = weight_ls*(3.0*dx*dx*nx);
    		        	    A_R4(index_least,10) = weight_ls*(3.0*dy*dy*ny);
    		        	    A_R4(index_least,11) = weight_ls*(3.0*dz*dz*nz);
    		        	    A_R4(index_least,12) = weight_ls*(2.0*dx*dy*nx+dx*dx*ny);
    		        	    A_R4(index_least,13) = weight_ls*(2.0*dx*dz*nx+dx*dx*nz);
    		        	    A_R4(index_least,14) = weight_ls*(dy*dy*nx+2.0*dx*dy*ny);
    		        	    A_R4(index_least,15) = weight_ls*(2.0*dy*dz*ny+dy*dy*nz);
    		        	    A_R4(index_least,16) = weight_ls*(dz*dz*nx+2.0*dx*dz*nz);
    		        	    A_R4(index_least,17) = weight_ls*(dz*dz*ny+2.0*dy*dz*nz);
    		        	    A_R4(index_least,18) = weight_ls*(dy*dz*nx+dx*dz*ny+dx*dy*nz);

    		            index_least++; 		

    		        }

				}
			}  // face loop for 4rd order contraint matrix ends
//			if(c == 36) std::cout<<"cls ROWS: "<<ROWS<<"\tindex_least: "<<index_least<<std::endl;			
/*
			if( g_i == 960 && Utilities::MPI::this_mpi_process(mpi_communicator) == 3) CLS_R4[c].initialize(A_R4, C_R4, true);
			else CLS_R4[c].initialize(A_R4, C_R4);
*/
			CLS_R4[c].initialize(A_R4, C_R4);
/*
			if(g_i == 9931 && Utilities::MPI::this_mpi_process(mpi_communicator) == 3) {
				pcout<<"c Nuemann: "<<g_i<<"\n";
				print_matrix(A_R4);
		//		pcout<<"c C_R4: "<<c<<"\n";
				print_matrix(C_R4);

			}		
*/
			// =====================================================================
            // r = 3 stencil 
            // =====================================================================

//	                C_R3.reinit(3, 9);
			ROWS = 6 - n_boundary_faces;
            C_R3.reinit(ROWS, 9); 
			C_R3 = 0.0;
			ROWS = cell_diagonal_neighbor_index[c].size() + 4*n_boundary_faces;
			A_R3.reinit(ROWS, 9); A_R3 = 0.0; index = 0; index_least = 0;

//					A_R3.reinit(ROWS + 9, 9); A_R3 = 0.0; index = 0; 
           	for (unsigned int f = 0; f < n_faces; f++) {			
				if(!cell->face(f)->at_boundary()) {	

					neighbor = cell->neighbor(f);
				    neighbor->get_dof_indices(local_neighbor_dof_indices);
					local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
		           	V_neighbor = Cell[local_index].measure();
	
                    for (unsigned int i = 0; i < n_quad_points; i++) {
                        q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						dx = (q_point(0) - x0)/h;
						dy = (q_point(1) - y0)/h;
						dz = (q_point(2) - z0)/h;
		        	    C_R3(index,0) += (1./V_neighbor)*j_w*dx;
		        	    C_R3(index,1) += (1./V_neighbor)*j_w*dy;
		        	    C_R3(index,2) += (1./V_neighbor)*j_w*dz;
		        	    C_R3(index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
		        	    C_R3(index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
		        	    C_R3(index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
		        	    C_R3(index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
		        	    C_R3(index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
		        	    C_R3(index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));

//								if(index == 0 && (g_i == 21376 || g_i == 21408) )
//									std::cout<<"g_i: "<<g_i<<"\tvalue 1: "<<A_R3(index,1)<<"\tj_w: "<<j_w<<"\tdy: "<<dy<<"\tV_neighbor: "<<V_neighbor<<std::endl;
                    }
//							if(c == 5376) pcout<<"index: "<<index<<"\tneigh: "<<local_index<<std::endl;
    	            index++; 		
				}
	
				else {

//			                QGauss<3-1> face_quadrature_formula2(n_boundary_quad_points);
	                QGauss<3-1> face_quadrature_formula2(n_boundary_quad_points);
//			                FEFaceValues<3> fv_face_values2 (mapping, fv, face_quadrature_formula2, update_quadrature_points | update_normal_vectors);
	                FEFaceValues<3> fv_face_values2 (mapping, fv, face_quadrature_formula2, update_quadrature_points | update_normal_vectors);
					fv_face_values2.reinit(cell, f);

					double nx, ny, nz;
		                for (unsigned int i = 0; i < 4; i++) {
							nx = fv_face_values2.normal_vector(i)[0];
							ny = fv_face_values2.normal_vector(i)[1];
							nz = fv_face_values2.normal_vector(i)[2];
							q_point = fv_face_values2.quadrature_point(i);

							dx = (q_point(0) - x0)/h;
							dy = (q_point(1) - y0)/h;
							dz = (q_point(2) - z0)/h;
	    		            A_R3(index_least,0) = weight_ls*(nx);
	    		            A_R3(index_least,1) = weight_ls*(ny);
	    		            A_R3(index_least,2) = weight_ls*(nz);
	    		            A_R3(index_least,3) = weight_ls*(2.0*dx*nx);
	    		            A_R3(index_least,4) = weight_ls*(2.0*dy*ny);
	    		            A_R3(index_least,5) = weight_ls*(2.0*dz*nz);
	    		            A_R3(index_least,6) = weight_ls*(dy*nx+dx*ny);
	    		            A_R3(index_least,7) = weight_ls*(dz*nx+dx*nz);
	    		            A_R3(index_least,8) = weight_ls*(dz*ny+dy*nz);
	    		            index_least++; 		
	
	    		        }
				}
			}  // face loop for 3rd order contraint matrix ends
			
					// First fill the vertex neighbors 
			for (unsigned int d = 0; d < cell_diagonal_neighbor_index[c].size(); ++d) {			
	
				local_index = global_to_local_index_map[cell_diagonal_neighbor_index[c][d] ];
		        V_neighbor = Cell[local_index].measure();
			
  		    	for (unsigned int i = 0; i < n_quad_points; i++) {
   		    	    q_point = Cell[local_index].cell_quadrature_point(i);
					j_w = Cell[local_index].jxw(i);
					dx = (q_point(0) - x0)/h;
					dy = (q_point(1) - y0)/h;
					dz = (q_point(2) - z0)/h;
  		    	    A_R3(index_least,0) += (1./V_neighbor)*j_w*dx;
   		    	    A_R3(index_least,1) += (1./V_neighbor)*j_w*dy;
   		    	    A_R3(index_least,2) += (1./V_neighbor)*j_w*dz;
   		    	    A_R3(index_least,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
   		    	    A_R3(index_least,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
   		    	    A_R3(index_least,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
   		    	    A_R3(index_least,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
   		    	    A_R3(index_least,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
   		    	    A_R3(index_least,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
   		    	}
//						if(c == 5376) pcout<<"index: "<<index<<"\tneigh: "<<local_index<<std::endl;
   		        index_least++; 
			}
					
			CLS_R3[c].initialize(A_R3, C_R3);
//					CLS_R3[c].initialize(A_R3);
/*
					if(g_i == 21376 ) {
						pcout<<"c: "<<c<<"\tg_i: "<<g_i<<"\nC_R3: \n";

						for(unsigned int i = 0; i < 3; ++i) {
							for(unsigned int j = 0; j < 9; ++j)
								pcout<<C_R3(i,j)<<"\t";
							pcout<<std::endl;
						}
						pcout<<"A_R3: sie: "<<index<<"\n";
						for(unsigned int i = 0; i < index; ++i) {
							for(unsigned int j = 0; j < 9; ++j)
								if(std::fabs(A_R3(i,j)) < 1e-14) pcout<<0.0<<"\t";
								else	pcout<<A_R3(i,j)<<"\t";
							pcout<<std::endl;
						}

					}
*/
/*
					if(g_i == 20352 || g_i == 21408) {
						std::cout<<"c: "<<c<<"\tg_i: "<<g_i<<"\nC_R3: \n";

						for(unsigned int i = 0; i < 4; ++i) {
							for(unsigned int j = 0; j < 9; ++j)
								std::cout<<C_R3(i,j)<<"\t";
							std::cout<<std::endl;
						}
						std::cout<<"A_R3: sie: "<<index<<"\n";
						for(unsigned int i = 0; i < index; ++i) {
							for(unsigned int j = 0; j < 9; ++j)
								std::cout<<A_R3(i,j)<<"\t";
							std::cout<<std::endl;
						}

					}
*/

				// =====================================================================
	            // r = 3 directiomal stencil 
	            // =====================================================================

        	for (unsigned int f = 0; f < n_faces; f++) {			

				is_admissible_R3_d[c][f] = false;

				unsigned int f_e = numbers::invalid_unsigned_int;
//					unsigned int n_boundary_faces_d = n_boundary_faces;

				if(f%2 == 0)
					f_e = f + 1;
				else
					f_e = f - 1;
				
//					if(cell->face(f_e)->at_boundary()) n_boundary_faces_d = n_boundary_faces_d - 1;

				unsigned int n_boundary_faces_d = 0, n_no_slip_boundary_d = 0;
				bool no_slip_boundary_d = false, transmissive_boundary_d = false;
			    for (unsigned int ff = 0; ff < n_faces; ++ff) {
					if(ff != f_e)
						if(cell->face(ff)->at_boundary()) {
							if(cell->face(ff)->boundary_id() == 2) {	
								n_no_slip_boundary_d += 1;
								no_slip_boundary_d = true;
							}
							if(cell->face(ff)->boundary_id() == 1 || cell->face(ff)->boundary_id() == 0) {	
								transmissive_boundary_d = true;
							}
							n_boundary_faces_d += 1;
						}
				}

				if(!cell->face(f)->at_boundary() && n_n_face[f]) {	

					if (cell_neighbor_index[c][f].size() >= 3) {
						if( (no_slip_boundary_d == no_slip_boundary) && (transmissive_boundary_d == transmissive_boundary) )
							is_admissible_R3_d[c][f] = true;
					}					
				}

				if(is_admissible_R3_d[c][f]) {
					ROWS = 5 - n_boundary_faces_d;
	                C_R3_d[f].reinit(ROWS, 9);
	
	                C_R3_d[f] = 0.0; 
	
   		            index = 0; index_least = 0; 
   		            ROWS = cell_neighbor_index[c][f].size() + 1 + n_boundary_faces_d*4; 
		
   		            A_R3_d[f].reinit(ROWS, 9); A_R3_d[f] = 0.0; 

	           		for (unsigned int ff = 0; ff < n_faces; ff++) {

						if(ff != f_e) {

							if(!cell->face(ff)->at_boundary()) {

								neighbor = cell->neighbor(ff);
							    neighbor->get_dof_indices(local_neighbor_dof_indices);
								local_index = global_to_local_index_map[local_neighbor_dof_indices[0] ];
				    	       	V_neighbor = Cell[local_index].measure();
		
    			        		for (unsigned int i = 0; i < n_quad_points; i++) {
    			        		    q_point = Cell[local_index].cell_quadrature_point(i);
									j_w = Cell[local_index].jxw(i);
									dx = (q_point(0) - x0)/h;
									dy = (q_point(1) - y0)/h;
									dz = (q_point(2) - z0)/h;
    			        		    C_R3_d[f](index,0) += (1./V_neighbor)*j_w*dx;
    			        		    C_R3_d[f](index,1) += (1./V_neighbor)*j_w*dy;
    			        		    C_R3_d[f](index,2) += (1./V_neighbor)*j_w*dz;
    			        		    C_R3_d[f](index,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
    			        		    C_R3_d[f](index,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
    			    	   		    C_R3_d[f](index,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
    			    	   		    C_R3_d[f](index,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
    			    	   		    C_R3_d[f](index,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
    			        		    C_R3_d[f](index,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
    			        		}
//										pcout<<"inside ff: "<<ff<<"\tindex: "<<index<<"\tsize: "<<std::endl;		
								index++;
		
							}
		
							else {
		
						        QGauss<3-1> face_quadrature_formula2(n_boundary_quad_points);
//				                FEFaceValues<3> fv_face_values2 (fv, face_quadrature_formula2, update_quadrature_points | update_normal_vectors);
					            FEFaceValues<3> fv_face_values2 (mapping, fv, face_quadrature_formula2, update_quadrature_points | update_normal_vectors);
								fv_face_values2.reinit(cell, ff);
  								double nx, ny, nz;
					            for (unsigned int i = 0; i < 4; i++) {
									nx = fv_face_values2.normal_vector(i)[0];
									ny = fv_face_values2.normal_vector(i)[1];
									nz = fv_face_values2.normal_vector(i)[2];
									q_point = fv_face_values2.quadrature_point(i);
									dx = (q_point(0) - x0)/h;
									dy = (q_point(1) - y0)/h;
									dz = (q_point(2) - z0)/h;
			    		            A_R3_d[f](index_least,0) = weight_ls*(nx);
    		    		            A_R3_d[f](index_least,1) = weight_ls*(ny);
    		    		            A_R3_d[f](index_least,2) = weight_ls*(nz);
    				        	    A_R3_d[f](index_least,3) = weight_ls*(2.0*dx*nx);
    				        	    A_R3_d[f](index_least,4) = weight_ls*(2.0*dy*ny);
    				        	    A_R3_d[f](index_least,5) = weight_ls*(2.0*dz*nz);
    				        	    A_R3_d[f](index_least,6) = weight_ls*(dy*nx+dx*ny);
    				        	    A_R3_d[f](index_least,7) = weight_ls*(dz*nx+dx*nz);
    				        	    A_R3_d[f](index_least,8) = weight_ls*(dz*ny+dy*nz);
									index_least++;
    		    		        }
							}
						}
					}

    	            // Least Squares Matrix 
    	        	    
					// vertex neighbor of cell at face f
					for (unsigned int d = 0; d < cell_neighbor_index[c][f].size(); ++d) {			
			
						local_index = global_to_local_index_map[cell_neighbor_index[c][f][d] ];
		    	        V_neighbor = Cell[local_index].measure();
						
    			    	for (unsigned int i = 0; i < n_quad_points; i++) {
	   			    	    q_point = Cell[local_index].cell_quadrature_point(i);
							j_w = Cell[local_index].jxw(i);
							dx = (q_point(0) - x0)/h;
							dy = (q_point(1) - y0)/h;
							dz = (q_point(2) - z0)/h;
    			    	    A_R3_d[f](index_least,0) += (1./V_neighbor)*j_w*dx;
	   			    	    A_R3_d[f](index_least,1) += (1./V_neighbor)*j_w*dy;
    			    	    A_R3_d[f](index_least,2) += (1./V_neighbor)*j_w*dz;
    			    	    A_R3_d[f](index_least,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
    			    	    A_R3_d[f](index_least,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
    			    	    A_R3_d[f](index_least,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
    			    	    A_R3_d[f](index_least,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
    			    	    A_R3_d[f](index_least,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
    			    	    A_R3_d[f](index_least,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
    			    	}
   		            
						index_least++; 
					}

					local_index = n_n_face_index[f];
    			    V_neighbor = Cell[local_index].measure();
    	        	    
    				for (unsigned int i = 0; i < n_quad_points; i++) {
    				    q_point = Cell[local_index].cell_quadrature_point(i);
						j_w = Cell[local_index].jxw(i);
						dx = (q_point(0) - x0)/h;
						dy = (q_point(1) - y0)/h;
						dz = (q_point(2) - z0)/h;
    				    A_R3_d[f](index_least,0) += (1./V_neighbor)*j_w*dx;
    				    A_R3_d[f](index_least,1) += (1./V_neighbor)*j_w*dy;
    				    A_R3_d[f](index_least,2) += (1./V_neighbor)*j_w*dz;
    				    A_R3_d[f](index_least,3) += (1./V_neighbor)*j_w*(dx*dx - WENO_poly_consts[c](3));
    				    A_R3_d[f](index_least,4) += (1./V_neighbor)*j_w*(dy*dy - WENO_poly_consts[c](4));
    				    A_R3_d[f](index_least,5) += (1./V_neighbor)*j_w*(dz*dz - WENO_poly_consts[c](5));
    				    A_R3_d[f](index_least,6) += (1./V_neighbor)*j_w*(dx*dy - WENO_poly_consts[c](6));
    				    A_R3_d[f](index_least,7) += (1./V_neighbor)*j_w*(dx*dz - WENO_poly_consts[c](7));
    				    A_R3_d[f](index_least,8) += (1./V_neighbor)*j_w*(dy*dz - WENO_poly_consts[c](8));
    				}

					CLS_R3_d[c][f].initialize(A_R3_d[f], C_R3_d[f]);		
		
				} // end of admissible loop
			} // end of 3rd order directional stencil

/*
 *         6-------7        6-------7
 *        /|       |       /       /|
 *       / |       |      /   f3  / |
 *      /  |  f5   |     /       /  |
 *     4   |       |    4-------5   |
 *     |f0 2-------3    |       |f1 3			z	y
 *     |  /       /     |       |  /			|  /    
 *     | /   f2  /      |  f4   | /				| /   
 *     |/       /       |       |/				|/      
 *     0-------1        0-------1               ------x
*/
/*
	    	for (unsigned int v = 0; v < 8; v++)
				A_R2_d[v].reinit(3,3);
			            
			// =====================================================================
            // r = 2 stencil 0 (f0 f2 & f4 )	
            // =====================================================================
            
            A_R2_d[0](0,0) = C_R3(0,0); A_R2_d[0](0,1) = C_R3(0,1); A_R2_d[0](0,2) = C_R3(0,2); // Row 1 f0
            A_R2_d[0](1,0) = C_R3(2,0); A_R2_d[0](1,1) = C_R3(2,1); A_R2_d[0](1,2) = C_R3(2,2); // Row 1 f2
            A_R2_d[0](2,0) = C_R3(4,0); A_R2_d[0](2,1) = C_R3(4,1); A_R2_d[0](2,2) = C_R3(4,2); // Row 1 f4

			std::cout<<"0: "<<c<<"\tC_R3: " <<std::endl;
			for(unsigned int i = 0; i < 3; ++i) 
				std::cout<<C_R3(2,i)<<"\t";
			std::cout<<std::endl;
				std::cout<<"0: "<<c<<std::endl;
			for(unsigned int i = 0; i < 3; ++i) {
				for(unsigned int j = 0; j < 3; ++j)
					std::cout<<A_R2_d[0](i,j)<<"\t";
				std::cout<<std::endl;
			}

            LU_R2_d[c][0].initialize(A_R2_d[0]);

//			if(A_R2_d[0](0,0) == 0.0 && A_R2_d[0](0,1) == 0.0 && A_R2_d[0](0,2) == 0.0 )
//				std::cout<<"0: "<<c<<std::endl;

			// =====================================================================
            // r = 2 stencil 1 (f1 f2 & f4 )	
            // =====================================================================
            
            A_R2_d[1](0,0) = C_R3(1,0); A_R2_d[1](0,1) = C_R3(1,1); A_R2_d[1](0,2) = C_R3(1,2); // Row 1 f1
            A_R2_d[1](1,0) = C_R3(2,0); A_R2_d[1](1,1) = C_R3(2,1); A_R2_d[1](1,2) = C_R3(2,2); // Row 1 f2
            A_R2_d[1](2,0) = C_R3(4,0); A_R2_d[1](2,1) = C_R3(4,1); A_R2_d[1](2,2) = C_R3(4,2); // Row 1 f4
            LU_R2_d[c][1].initialize(A_R2_d[1]);
//			std::cout<<"1: "<<std::endl;
//			if(A_R2_d[1](0,0) == 0.0 && A_R2_d[1](0,1) == 0.0 && A_R2_d[1](0,2) == 0.0 )
//				std::cout<<"1: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 2 (f0 f2 & f5 )	
            // =====================================================================
            
            A_R2_d[2](0,0) = C_R3(0,0); A_R2_d[2](0,1) = C_R3(0,1); A_R2_d[2](0,2) = C_R3(0,2); // Row 1 f0
            A_R2_d[2](1,0) = C_R3(2,0); A_R2_d[2](1,1) = C_R3(2,1); A_R2_d[2](1,2) = C_R3(2,2); // Row 1 f2
            A_R2_d[2](2,0) = C_R3(5,0); A_R2_d[2](2,1) = C_R3(5,1); A_R2_d[2](2,2) = C_R3(5,2); // Row 1 f5
            LU_R2_d[c][2].initialize(A_R2_d[2]);
//			std::cout<<"2: "<<std::endl;
//			if(A_R2_d[2](0,0) == 0.0 && A_R2_d[2](0,1) == 0.0 && A_R2_d[2](0,2) == 0.0 )
//				std::cout<<"2: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 3 (f1 f2 & f5 )	
            // =====================================================================
            
            A_R2_d[3](0,0) = C_R3(1,0); A_R2_d[3](0,1) = C_R3(1,1); A_R2_d[3](0,2) = C_R3(1,2); // Row 1 f1
            A_R2_d[3](1,0) = C_R3(2,0); A_R2_d[3](1,1) = C_R3(2,1); A_R2_d[3](1,2) = C_R3(2,2); // Row 1 f2
            A_R2_d[3](2,0) = C_R3(5,0); A_R2_d[3](2,1) = C_R3(5,1); A_R2_d[3](2,2) = C_R3(5,2); // Row 1 f5
            LU_R2_d[c][3].initialize(A_R2_d[3]);
//			std::cout<<"3: "<<std::endl;
//			if(A_R2_d[3](0,0) == 0.0 && A_R2_d[3](0,1) == 0.0 && A_R2_d[3](0,2) == 0.0 )
//				std::cout<<"3: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 4 (f0 f3 & f4 )	
            // =====================================================================
            
            A_R2_d[4](0,0) = C_R3(0,0); A_R2_d[4](0,1) = C_R3(0,1); A_R2_d[4](0,2) = C_R3(0,2); // Row 1 f0
            A_R2_d[4](1,0) = C_R3(3,0); A_R2_d[4](1,1) = C_R3(3,1); A_R2_d[4](1,2) = C_R3(3,2); // Row 1 f3
            A_R2_d[4](2,0) = C_R3(4,0); A_R2_d[4](2,1) = C_R3(4,1); A_R2_d[4](2,2) = C_R3(4,2); // Row 1 f4
            LU_R2_d[c][4].initialize(A_R2_d[4]);
//			std::cout<<"4: "<<std::endl;
//			if(A_R2_d[4](0,0) == 0.0 && A_R2_d[4](0,1) == 0.0 && A_R2_d[4](0,2) == 0.0 )
//				std::cout<<"4: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 5 (f1 f3 & f4 )	
            // =====================================================================
            
            A_R2_d[5](0,0) = C_R3(1,0); A_R2_d[5](0,1) = C_R3(1,1); A_R2_d[5](0,2) = C_R3(1,2); // Row 1 f1
            A_R2_d[5](1,0) = C_R3(3,0); A_R2_d[5](1,1) = C_R3(3,1); A_R2_d[5](1,2) = C_R3(3,2); // Row 1 f3
            A_R2_d[5](2,0) = C_R3(4,0); A_R2_d[5](2,1) = C_R3(4,1); A_R2_d[5](2,2) = C_R3(4,2); // Row 1 f4
            LU_R2_d[c][5].initialize(A_R2_d[5]);
//			std::cout<<"5: "<<std::endl;
//			if(A_R2_d[5](0,0) == 0.0 && A_R2_d[5](0,1) == 0.0 && A_R2_d[5](0,2) == 0.0 )
//				std::cout<<"5: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 6 (f0 f3 & f5 )	
            // =====================================================================
            
            A_R2_d[6](0,0) = C_R3(0,0); A_R2_d[6](0,1) = C_R3(0,1); A_R2_d[6](0,2) = C_R3(0,2); // Row 1 f0
            A_R2_d[6](1,0) = C_R3(3,0); A_R2_d[6](1,1) = C_R3(3,1); A_R2_d[6](1,2) = C_R3(3,2); // Row 1 f3
            A_R2_d[6](2,0) = C_R3(5,0); A_R2_d[6](2,1) = C_R3(5,1); A_R2_d[6](2,2) = C_R3(5,2); // Row 1 f5
            LU_R2_d[c][6].initialize(A_R2_d[6]);
//			std::cout<<"6: "<<std::endl;
//			if(A_R2_d[6](0,0) == 0.0 && A_R2_d[6](0,1) == 0.0 && A_R2_d[6](0,2) == 0.0 )
//				std::cout<<"6: "<<c<<std::endl;
			// =====================================================================
            // r = 2 stencil 7 (f1 f3 & f5 )	
            // =====================================================================
            
            A_R2_d[7](0,0) = C_R3(1,0); A_R2_d[7](0,1) = C_R3(1,1); A_R2_d[7](0,2) = C_R3(1,2); // Row 1 f1
            A_R2_d[7](1,0) = C_R3(3,0); A_R2_d[7](1,1) = C_R3(3,1); A_R2_d[7](1,2) = C_R3(3,2); // Row 1 f3
            A_R2_d[7](2,0) = C_R3(5,0); A_R2_d[7](2,1) = C_R3(5,1); A_R2_d[7](2,2) = C_R3(5,2); // Row 1 f5
            LU_R2_d[c][7].initialize(A_R2_d[7]);
//			if(A_R2_d[7](0,0) == 0.0 && A_R2_d[7](0,1) == 0.0 && A_R2_d[7](0,2) == 0.0 )
//				std::cout<<"7: "<<c<<std::endl;
*/



		} //  End of Boundary cell loop

    } // End of cell loop 
    
    pcout << "Done!" << std::endl;

} // End of function  

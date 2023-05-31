#include "../include/Weno32.h"


// Compute the rhs vectors

void Weno3_2D::compute_rhs() {

    reconstruct();
    
    unsigned int faces_per_cell = GeometryInfo<2>::faces_per_cell;

	// Loop over all the cells
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

	// For obtaining face normal vectors
	QGauss<2-1> face_quadrature_formula(2);
	FEFaceValues<2> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors);

    Point<2> face_quadrature_point_1;
    Point<2> face_quadrature_point_2;

    double V_c, S_f; // Volume of the cell and surface area of the face
	
	double nx1, ny1;   // Face normal vectors
	double nx2, ny2; 
	
	Tensor<1,2> face_normal_vector1; // Face normal vector
	Tensor<1,2> face_normal_vector2; // Face normal vector
	
    Vector<double> UL1(4); Vector<double> UR1(4); // Solving the Riemann Problem
    Vector<double> UL2(4); Vector<double> UR2(4); // Solving the Riemann Problem

	Vector<double> WL1(4); Vector<double> WR1(4); // Solving the Riemann Problem
    Vector<double> WL2(4); Vector<double> WR2(4); // Solving the Riemann Problem
	
    Vector<double> F(4);
    
    bool boundary; unsigned int global_face_index;   
    
    unsigned int n_faces = triangulation.n_active_faces(); 

    std::vector< Vector<double> > Flux1(n_faces); 
    std::vector< Vector<double> > Flux2(n_faces);
    std::vector< bool > did_not_compute_flux_for_the_face(n_faces); 
    
    for (unsigned int f = 0; f < n_faces; f++) {
        Flux1[f].reinit(4);
        Flux2[f].reinit(4);
        did_not_compute_flux_for_the_face[f] = true; 
    }
    
	for (unsigned int c = 0; cell != endc; ++cell, ++c) {
    
		V_c = cell->measure();

		rhs1(c) = 0.0;
		rhs2(c) = 0.0;
		rhs3(c) = 0.0;
		rhs4(c) = 0.0;
        
        
        for (unsigned int f = 0; f < faces_per_cell; ++f) {
            
            global_face_index = cell->face_index(f);
            S_f = cell->face(f)->measure();
            
            if(did_not_compute_flux_for_the_face[global_face_index]) {
                
                // Get some geometry info
        
                
                fv_face_values.reinit(cell, f);
                
                face_normal_vector1 = fv_face_values.normal_vector(0); // Zero corresponds to the presence of only one quadrature point
                nx1 = face_normal_vector1[0]; ny1 = face_normal_vector1[1];
                
                face_normal_vector2 = fv_face_values.normal_vector(1);
                nx2 = face_normal_vector2[0]; ny2 = face_normal_vector2[1];
                
                face_quadrature_point_1 = fv_face_values.quadrature_point(0); 
                face_quadrature_point_2 = fv_face_values.quadrature_point(1); 
                
                // Left face
                UL1(0) = evaluate_weno_polynomial(coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point_1);
                UL1(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point_1);
                UL1(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point_1);
                UL1(3) = evaluate_weno_polynomial(coeffs_E[c], WENO_poly_consts[c], face_quadrature_point_1);
                
                UL2(0) = evaluate_weno_polynomial(coeffs_RHO[c], WENO_poly_consts[c], face_quadrature_point_2);
                UL2(1) = evaluate_weno_polynomial(coeffs_RHO_U[c], WENO_poly_consts[c], face_quadrature_point_2);
                UL2(2) = evaluate_weno_polynomial(coeffs_RHO_V[c], WENO_poly_consts[c], face_quadrature_point_2);
                UL2(3) = evaluate_weno_polynomial(coeffs_E[c], WENO_poly_consts[c], face_quadrature_point_2);
                
                WL1 = conserved_to_primitive(UL1); WL2 = conserved_to_primitive(UL2);
                
                if (cell->face(f)->at_boundary()) {
                    
                    // Apply transmissive boundary conditions
                    
                    WR1 = WL1; WR2 = WL2;
                    
                    boundary = true; 
                }
                
                else {
                    
                    // Get the right state values
                    
                    UR1(0) = evaluate_weno_polynomial(coeffs_RHO[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1);
                    UR1(1) = evaluate_weno_polynomial(coeffs_RHO_U[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1);
                    UR1(2) = evaluate_weno_polynomial(coeffs_RHO_V[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1);
                    UR1(3) = evaluate_weno_polynomial(coeffs_E[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_1);
                    
                    UR2(0) = evaluate_weno_polynomial(coeffs_RHO[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2);
                    UR2(1) = evaluate_weno_polynomial(coeffs_RHO_U[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2);
                    UR2(2) = evaluate_weno_polynomial(coeffs_RHO_V[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2);
                    UR2(3) = evaluate_weno_polynomial(coeffs_E[cell->neighbor_index(f)], WENO_poly_consts[cell->neighbor_index(f)],  face_quadrature_point_2);
                    
                    WR1 = conserved_to_primitive(UR1); WR2 = conserved_to_primitive(UR2);
                    
                    boundary = false; 
                }
        
                
                Flux1[global_face_index] = rotated_HLLC_riemann_solver(WL1, WR1, nx1, ny1, face_quadrature_point_1, boundary);
                Flux2[global_face_index] = rotated_HLLC_riemann_solver(WL2, WR2, nx2, ny2, face_quadrature_point_2, boundary);
                
                
                did_not_compute_flux_for_the_face[global_face_index] = false; 
            }
            
            else {
                
                Flux1[global_face_index] *= -1.0; 
                Flux2[global_face_index] *= -1.0;
            }
            
            
            F(0) = 0.5*(Flux1[global_face_index](0) + Flux2[global_face_index](0)); 
            F(1) = 0.5*(Flux1[global_face_index](1) + Flux2[global_face_index](1)); 
            F(2) = 0.5*(Flux1[global_face_index](2) + Flux2[global_face_index](2));
            F(3) = 0.5*(Flux1[global_face_index](3) + Flux2[global_face_index](3));

            // Add it to the rhs vectors

            rhs1(c) += (-1.0/V_c)*(F(0)*S_f);
            rhs2(c) += (-1.0/V_c)*(F(1)*S_f);
            rhs3(c) += (-1.0/V_c)*(F(2)*S_f);
            rhs4(c) += (-1.0/V_c)*(F(3)*S_f);
        }
    }
}


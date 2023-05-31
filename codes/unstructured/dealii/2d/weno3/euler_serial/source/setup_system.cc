#include "../include/Weno32.h"


// Setup the system - allocate the necessary memory 

void Weno3_2D::setup_system() {

    std::cout << "Setting up system" << std::endl;

	dof_handler.distribute_dofs (fv);

	RHO.reinit(dof_handler.n_dofs());
	RHO_U.reinit(dof_handler.n_dofs());
	RHO_V.reinit(dof_handler.n_dofs());
	E.reinit(dof_handler.n_dofs());

	rhs1.reinit(dof_handler.n_dofs());
	rhs2.reinit(dof_handler.n_dofs());
	rhs3.reinit(dof_handler.n_dofs());
	rhs4.reinit(dof_handler.n_dofs());

    coeffs_RHO.resize(dof_handler.n_dofs());
    coeffs_RHO_U.resize(dof_handler.n_dofs());
    coeffs_RHO_V.resize(dof_handler.n_dofs());
    coeffs_E.resize(dof_handler.n_dofs());
    
    WENO_poly_consts.resize(dof_handler.n_dofs()); 
    IS_constants.resize(dof_handler.n_dofs());
	Cell.resize(dof_handler.n_dofs());
	Stencil.resize(dof_handler.n_dofs());
    
    CLS_R3.resize(dof_handler.n_dofs());
    LU_R21.resize(dof_handler.n_dofs());
    LU_R22.resize(dof_handler.n_dofs());
    LU_R23.resize(dof_handler.n_dofs());
    LU_R24.resize(dof_handler.n_dofs());
    
    is_corner_cell.resize(dof_handler.n_dofs());
    
    
    for (unsigned int i = 0; i < dof_handler.n_dofs(); i++) {
        coeffs_RHO[i].reinit(6);
        coeffs_RHO_U[i].reinit(6);
        coeffs_RHO_V[i].reinit(6);
        coeffs_E[i].reinit(6);
        WENO_poly_consts[i].reinit(5); 
        IS_constants[i].reinit(6);
    }
    
        // Initialize the wall boundary matrices  
    
    // Get the number of cells at the wall boundary 
    
    unsigned int n_wall_boundary_cells = 0; 
    
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
    
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {
        
        if (cell->at_boundary()) {
            
            for (unsigned int f=0; f < GeometryInfo<2>::faces_per_cell; ++f) {
                
                if (cell->face(f)->at_boundary()) {
                    
                    if (cell->face(f)->boundary_id() == 2) { // boundary_id 2 corresponds to the reflecting wall 
                        
                        wall_boundary_global_index_map[c] = n_wall_boundary_cells; 
                        
                        n_wall_boundary_cells ++; 
                        
                    }
                    
                }
            }
        }
    }
            
    std::cout << "Number of wall boundary cells = " <<  n_wall_boundary_cells << std::endl; 
    
    CLS_R3_slip.resize(n_wall_boundary_cells);
    CLS_R2_slip.resize(n_wall_boundary_cells);

    std::cout << "Done!" << std::endl;  
    std::cout << "===========================" << std::endl;
}

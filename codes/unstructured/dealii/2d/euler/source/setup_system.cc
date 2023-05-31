#include "../include/Euler.h"


// Setup the system - allocate the necessary memory 

void Euler_2D::setup_system() {

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

    std::cout << "Done!" << std::endl;  
    std::cout << "===========================" << std::endl;
}

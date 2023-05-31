#include "../include/Weno32.h"

DoFHandler<2>::active_cell_iterator return_cell_pointer(const DoFHandler<2>& dof_handler, unsigned int i) {
	DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();
	
	for (unsigned int c = 0; cell != endc; ++cell, ++c) {
	
		if (c==i) {
			return cell;
		}
		
	}
	
	return cell; 
}

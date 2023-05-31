#include "../include/Weno32.h"


void Weno3_2D::plot_tecplot() const {
	
	
	unsigned int index; 
	
	std::vector<unsigned int> local_dof_index(1); 
	
    std::ofstream fout_RHO("RHO_" + Utilities::int_to_string( Utilities::MPI::this_mpi_process(MPI_COMM_WORLD), 4) + ".dat");
    
    DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

	for (; cell != endc; ++cell) {
	    if (cell->is_locally_owned()){
	        
			index = cell->active_cell_index(); 
			
			cell->get_dof_indices(local_dof_index);

	        fout_RHO << index << " " << local_RHO(local_dof_index[0]) << std::endl; 

		}
	}
    
    fout_RHO.close(); 
} 

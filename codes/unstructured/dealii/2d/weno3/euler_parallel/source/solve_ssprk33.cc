#include "../include/Weno32.h"

// Solve using SSPRK(3,3)

void Weno3_2D::solve_ssprk33() {

    TimerOutput::Scope t(computing_timer, "Solve by ssprk33");

    unsigned int count = 0;
    
    Vector<double> RHO_old(n_locally_cells);
    Vector<double> RHO_U_old(n_locally_cells);
	Vector<double> RHO_V_old(n_locally_cells);
	Vector<double> E_old(n_locally_cells); 

	compute_time_step_based_on_cfl_number();

	unsigned int output_count;
	output_count = 0.25 / dt ;

	while (time < finalTime) {

		compute_time_step_based_on_cfl_number();

		if (count%output_count == 0 ) {
            output_results();
        }

		if (count%100 == 0 ) {
            //restart();
        }

		if ((count+5)%100 == 0 ) {
            //restart_r();
        }

		time += dt;       
        
        pcout << "time = " << time <<"\tdt = "<<dt<< "\tFinal time = " << finalTime<< std::endl;
        count++;

		unsigned int c, g_i ;

		DoFHandler<2>::active_cell_iterator cell = dof_handler.begin_active();
		DoFHandler<2>::active_cell_iterator endc = dof_handler.end();

		for (; cell != endc; ++cell) {
		    if (cell->is_locally_owned()){
		        cell->get_dof_indices(local_dof_indices);
				g_i = local_dof_indices[0] ;
				c = global_to_local_index_map[local_dof_indices[0] ];

				RHO_old(c) = RHO(g_i);
				RHO_U_old(c) = RHO_U(g_i);
				RHO_V_old(c) = RHO_V(g_i);
				E_old(c) = E(g_i);
			}	
		}

		// SSPRK Step 1

        compute_rhs();

		cell = dof_handler.begin_active();

		for (; cell != endc; ++cell) {
		    if (cell->is_locally_owned()){
		        cell->get_dof_indices(local_dof_indices);
				c = global_to_local_index_map[local_dof_indices[0] ];
				g_i = local_dof_indices[0];

            	local_RHO(g_i) = local_RHO(g_i) + dt*rhs1(c);
	            local_RHO_U(g_i) = local_RHO_U(g_i) + dt*rhs2(c);
    	        local_RHO_V(g_i) = local_RHO_V(g_i) + dt*rhs3(c) + dt*local_RHO(g_i);
    	        local_E(g_i) = local_E(g_i) + dt*rhs4(c) + dt*local_RHO_V(g_i);
        
			}
		}

	local_RHO.compress(VectorOperation::insert);
	local_RHO_U.compress(VectorOperation::insert);
	local_RHO_V.compress(VectorOperation::insert);
	local_E.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;

		// SSPRK Step 2
		
        
		compute_rhs();

		cell = dof_handler.begin_active();

		for (; cell != endc; ++cell) {
		    if (cell->is_locally_owned()){
		        cell->get_dof_indices(local_dof_indices);
				c = global_to_local_index_map[local_dof_indices[0] ];
				g_i = local_dof_indices[0];

	            local_RHO(g_i)   = (3./4.)*RHO_old(c)   + (1./4.)*local_RHO(g_i)   + (1./4.)*dt*rhs1(c);
	            local_RHO_U(g_i) = (3./4.)*RHO_U_old(c) + (1./4.)*local_RHO_U(g_i) + (1./4.)*dt*rhs2(c);
	            local_RHO_V(g_i) = (3./4.)*RHO_V_old(c) + (1./4.)*local_RHO_V(g_i) + (1./4.)*dt*rhs3(c) + (1./4.)*dt*local_RHO(g_i);
	            local_E(g_i)     = (3./4.)*E_old(c)     + (1./4.)*local_E(g_i)     + (1./4.)*dt*rhs4(c) + (1./4.)*dt*local_RHO_V(g_i);
        
	        }
		}

	local_RHO.compress(VectorOperation::insert);
	local_RHO_U.compress(VectorOperation::insert);
	local_RHO_V.compress(VectorOperation::insert);
	local_E.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;

		// SSPRK Step 3

		compute_rhs();

		cell = dof_handler.begin_active();

		for (; cell != endc; ++cell) {
		    if (cell->is_locally_owned()){
		        cell->get_dof_indices(local_dof_indices);
				c = global_to_local_index_map[local_dof_indices[0] ];
				g_i = local_dof_indices[0];

	            local_RHO(g_i)   = (1./3.)*RHO_old(c)   + (2./3.)*local_RHO(g_i)   + (2./3.)*dt*rhs1(c);
	            local_RHO_U(g_i) = (1./3.)*RHO_U_old(c) + (2./3.)*local_RHO_U(g_i) + (2./3.)*dt*rhs2(c);
	            local_RHO_V(g_i) = (1./3.)*RHO_V_old(c) + (2./3.)*local_RHO_V(g_i) + (2./3.)*dt*rhs3(c) + (2./3.)*dt*local_RHO(g_i);
	            local_E(g_i)     = (1./3.)*E_old(c)     + (2./3.)*local_E(g_i)     + (2./3.)*dt*rhs4(c) + (2./3.)*dt*local_RHO_V(g_i);
            
			}
		}		

	local_RHO.compress(VectorOperation::insert);
	local_RHO_U.compress(VectorOperation::insert);
	local_RHO_V.compress(VectorOperation::insert);
	local_E.compress(VectorOperation::insert);

		RHO = local_RHO;
		RHO_U = local_RHO_U;
		RHO_V = local_RHO_V;
		E = local_E;

    		
    }
 
    output_results();
    //restart(); 
}
 

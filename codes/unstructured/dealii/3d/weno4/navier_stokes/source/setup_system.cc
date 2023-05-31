#include "../include/Weno432.h"


// Setup the system - allocate the necessary memory 

void Weno4_3D::setup_system() {

//    auto start = std::chrono::system_clock::now();

    pcout << "Setting up system" << std::endl;

	time = 0.0;

	dof_handler.distribute_dofs (fv);

    locally_owned_dofs = dof_handler.locally_owned_dofs ();
	locally_relevant_dofs = locally_owned_dofs;

	dofs_per_cell = fv.dofs_per_cell;
	local_dof_indices.resize(dofs_per_cell);
	local_neighbor_dof_indices.resize(dofs_per_cell);

	n_locally_cells = dof_handler.n_locally_owned_dofs();
//	is_no_slip_cell.resize(2*n_locally_cells, false);
	unsigned int n_locally_cells_min = Utilities::MPI::min (n_locally_cells, MPI_COMM_WORLD);
	unsigned int n_locally_cells_max = Utilities::MPI::max (n_locally_cells, MPI_COMM_WORLD);
   	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
   	    fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
   	    fout_convergence.precision(7) ;
	
   	    const std::string filename = "log.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

   		fout_convergence <<"no of cells on processor "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<" is : "<<n_locally_cells<<"\tn_locally_cells_min: "<<n_locally_cells_min<<"\tn_locally_cells_max: "<<n_locally_cells_max<<std::endl;
   	    fout_convergence.close();
	}

   		pcout <<"no of cells on processor "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<" is : "<<n_locally_cells<<"\tn_locally_cells_min: "<<n_locally_cells_min<<"\tn_locally_cells_max: "<<n_locally_cells_max<<std::endl;

	std::vector<bool> is_extra_ghost_cell, is_post_processing_cell; 
	is_relevant_cell.resize(triangulation.n_active_cells(), false );
	is_extra_ghost_cell.resize(triangulation.n_active_cells(), false );
	is_ghost_cell.resize(triangulation.n_active_cells(), false );
	is_post_processing_cell.resize(triangulation.n_active_cells(), false );

    DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<3>::active_cell_iterator endc = dof_handler.end();

	std::set<unsigned int> vertex_index;	
	std::set<unsigned int> face_index;	

	DoFHandler<3>::active_cell_iterator neighbor = cell, neighbor1 = cell, neighbor2 = cell;
/*
    for (unsigned int c = 0; c < triangulation.n_active_cells(); ++c) {
		is_relevant_cell[c] = false;
		is_extra_ghost_cell[c] = false;
		is_ghost_cell[c] = false;
		is_post_processing_cell[c] = false;
	}
*/
	std::vector< std::vector< std::set<DoFHandler<3>::active_cell_iterator> > >cell_neighbor_iterator;
	std::vector< std::set<DoFHandler<3>::active_cell_iterator> > cell_all_neighbor_iterator;
	std::vector< std::vector< std::set<DoFHandler<3>::active_cell_iterator> > >cell_neighbor_neighbor_iterator;
	std::vector< std::set<DoFHandler<3>::active_cell_iterator> > cell_diagonal_neighbor_iterator;

	cell = dof_handler.begin_active();
	unsigned int local_index = 0;
//	n_wall_boundary_cells = 0;

    for (; cell != endc; ++cell) {
       
		if (cell->is_locally_owned()){

	        cell->get_dof_indices(local_dof_indices);

			is_relevant_cell[local_dof_indices[0] ] = true;
			global_to_local_index_map[local_dof_indices[0] ] = local_index;
			local_to_global_index_map.push_back(local_dof_indices[0] );
			local_index_to_iterator.push_back(cell);

	        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
				vertex_index.insert(cell->vertex_index(i));		

	        for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f){

				face_index.insert(cell->face_index(f));

				if( !cell->face(f)->at_boundary() ) {
					neighbor = cell->neighbor(f);
			        neighbor->get_dof_indices(local_neighbor_dof_indices);
					is_relevant_cell[local_neighbor_dof_indices[0] ] = true;
				}

			}

			local_index++;
		}
	}

	cell = dof_handler.begin_active();

    for (; cell != endc; ++cell) {

        cell->get_dof_indices(local_dof_indices);

		if(is_relevant_cell[local_dof_indices[0] ] && !(cell->is_locally_owned())) {

			global_to_local_index_map[local_dof_indices[0] ] = local_index;
			is_extra_ghost_cell[local_dof_indices[0] ] = true;
			local_to_global_index_map.push_back(local_dof_indices[0] );
			local_index_to_iterator.push_back(cell);

			local_index++;

	        for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f) {
				if( !cell->face(f)->at_boundary() ){
					neighbor = cell->neighbor(f);
			        neighbor->get_dof_indices(local_neighbor_dof_indices);
					if(!neighbor->is_locally_owned()) {
						is_extra_ghost_cell[local_neighbor_dof_indices[0] ] = true;
				        for (unsigned int f1=0; f1<GeometryInfo<3>::faces_per_cell; ++f1) {
							if( !neighbor->face(f1)->at_boundary() ) {
								neighbor1 = neighbor->neighbor(f1);
						        neighbor1->get_dof_indices(local_neighbor_dof_indices);
								if(!neighbor1->is_locally_owned()) {
									is_extra_ghost_cell[local_neighbor_dof_indices[0] ] = true;
							        for (unsigned int f2=0; f2<GeometryInfo<3>::faces_per_cell; ++f2) {
										if( !neighbor1->face(f2)->at_boundary() ) {
											neighbor2 = neighbor1->neighbor(f2);
								    	    neighbor2->get_dof_indices(local_neighbor_dof_indices);
											if(!neighbor2->is_locally_owned())
												is_extra_ghost_cell[local_neighbor_dof_indices[0] ] = true;
										}
									}
								}
							}
						}
					}
				}
			}
		}

	}

	n_relevant_cells = local_index;

	pcout<<"no of relevant cells: "<<n_relevant_cells<<std::endl;

	cell_neighbor_iterator.resize(n_relevant_cells);
	cell_diagonal_neighbor_iterator.resize(n_relevant_cells);
	cell_neighbor_neighbor_iterator.resize(n_relevant_cells);

    
    for (unsigned int i = 0; i < n_relevant_cells; i++) {
		cell_neighbor_iterator[i].resize(6);
		cell_neighbor_neighbor_iterator[i].resize(6);
    }

	cell = dof_handler.begin_active();

    for (; cell != endc; ++cell) {
       
        cell->get_dof_indices(local_dof_indices);

		if (is_extra_ghost_cell[local_dof_indices[0] ])
	        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
				vertex_index.insert(cell->vertex_index(i));
	}

	n_vertices = vertex_index.size();
	n_faces = face_index.size();

	typename std::set<unsigned int>::iterator ver_iter = vertex_index.begin(), ver_end = vertex_index.end();

	unsigned int local_vertex = 0;
	for (; ver_iter != ver_end; ++ver_iter) {
		vertex_map[*ver_iter] = local_vertex;
		local_vertex++;
	}

	typename std::set<unsigned int>::iterator face_iter = face_index.begin(), face_end = face_index.end();

	unsigned int local_face_index = 0;
	for (; face_iter != face_end; ++face_iter) {
		face_index_map[*face_iter] = local_face_index;
		local_face_index++;
	}

	cell = dof_handler.begin_active();

	std::vector< std::set<DoFHandler<3>::active_cell_iterator> > vertex_to_cell_iterator(n_vertices);
	vertex_to_cell_index.resize(n_vertices);

	cell = dof_handler.begin_active();

    for (; cell!=endc; ++cell) {
        cell->get_dof_indices(local_dof_indices);
		if(cell->is_locally_owned() || is_extra_ghost_cell[local_dof_indices[0] ] ) {
	        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i)
				vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].insert(cell);
		}
    }

	typedef typename std::set<DoFHandler<3>::active_cell_iterator>::iterator ver_cell_iter;

    for (unsigned int i = 0; i < n_vertices; ++i) {

		ver_cell_iter adjacent_cell = vertex_to_cell_iterator[i].begin();						
		ver_cell_iter end_cell = vertex_to_cell_iterator[i].end();

		unsigned int n_sharing_cell = vertex_to_cell_iterator[i].size();
		vertex_to_cell_index[i].resize(n_sharing_cell);

		for (unsigned int n = 0; adjacent_cell != end_cell; ++adjacent_cell, ++n) {	
			cell = *adjacent_cell;
			cell->get_dof_indices(local_dof_indices);
			vertex_to_cell_index[i][n] = local_dof_indices[0];	
		}

    }

	cell = dof_handler.begin_active();
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {

		if (cell->is_locally_owned()){
	        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i) {

				ver_cell_iter adjacent_cell = vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].begin();						
				ver_cell_iter end_cell = vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].end();
				for(;adjacent_cell != end_cell; ++adjacent_cell) {
					neighbor = *adjacent_cell;
					neighbor->get_dof_indices(local_dof_indices);
					if(!is_relevant_cell[local_dof_indices[0] ] ) {
						is_post_processing_cell[local_dof_indices[0]] = true;
					}					
				}
			}
		}
	}


	cell = dof_handler.begin_active();
    for (unsigned int c = 0; cell != endc; ++cell, ++c) {


        cell->get_dof_indices(local_dof_indices);
		if(is_post_processing_cell[local_dof_indices[0] ]) {

			global_to_local_index_map[local_dof_indices[0]] = local_index;
			local_to_global_index_map.push_back(local_dof_indices[0]);
			local_index_to_iterator.push_back(cell); 
			locally_relevant_dofs.add_index(local_dof_indices[0]);
			local_index++;
		}
	}

	n_post_process_cell = local_index;

	pcout<<"no of post process cells: "<<n_post_process_cell<<std::endl;

	cell_all_neighbor_iterator.resize(n_post_process_cell);

	cell = dof_handler.begin_active();
    for (; cell != endc; ++cell) {

        cell->get_dof_indices(local_dof_indices);
        
		if (is_relevant_cell[local_dof_indices[0] ]){

			unsigned int local_i = global_to_local_index_map[local_dof_indices[0] ];

	        for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f) {

		        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_face; ++i){

					ver_cell_iter adjacent_cell = vertex_to_cell_iterator[ vertex_map[cell->face(f)->vertex_index(i)] ].begin();						
					ver_cell_iter end_cell = vertex_to_cell_iterator[ vertex_map[cell->face(f)->vertex_index(i)] ].end();

					for(;adjacent_cell != end_cell; ++adjacent_cell) {
						cell_neighbor_iterator[local_i][f].insert(*adjacent_cell);
						cell_all_neighbor_iterator[local_i].insert(*adjacent_cell);
						cell_diagonal_neighbor_iterator[local_i].insert(*adjacent_cell);
					}
				}

				if(!cell->face(f)->at_boundary() ) {

					neighbor = cell->neighbor(f);
					unsigned int neighbor_face_index = cell->neighbor_of_neighbor(f);

					unsigned int neighbor_neighbor_face_index = numbers::invalid_unsigned_int;

					if(neighbor_face_index%2 == 0) 
						neighbor_neighbor_face_index = neighbor_face_index + 1;
					else 
						neighbor_neighbor_face_index = neighbor_face_index - 1;
/*
					if(neighbor_face_index == 0) 
						neighbor_neighbor_face_index = 1;
					else if(neighbor_face_index == 1) 
						neighbor_neighbor_face_index = 0;
					else if(neighbor_face_index == 2) 
						neighbor_neighbor_face_index = 3;
					else if(neighbor_face_index == 3) 
						neighbor_neighbor_face_index = 2;
					else if(neighbor_face_index == 4) 
						neighbor_neighbor_face_index = 5;
					else if(neighbor_face_index == 5) 
						neighbor_neighbor_face_index = 4;
*/
					if(!neighbor->face(neighbor_neighbor_face_index)->at_boundary() )
						cell_neighbor_neighbor_iterator[local_i][f].insert(neighbor->neighbor(neighbor_neighbor_face_index));	

				}

			}

	        for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f) {

				unsigned int ff = numbers::invalid_unsigned_int;

				if(f%2 == 0)
					ff = f + 1;
				else
					ff = f - 1;

		        for (unsigned int f3=0; f3<GeometryInfo<3>::faces_per_cell; ++f3) {
					if(f3 != ff)
						if(!cell->face(f3)->at_boundary() ) 	
							cell_neighbor_iterator[local_i][f].erase(cell->neighbor(f3));
				}

				if(!cell->face(f)->at_boundary() )	
					cell_diagonal_neighbor_iterator[local_i].erase(cell->neighbor(f));
				cell_neighbor_iterator[local_i][f].erase(cell);

			}
			cell_all_neighbor_iterator[local_i].erase(cell);
			cell_diagonal_neighbor_iterator[local_i].erase(cell);
		}

		if (is_post_processing_cell[local_dof_indices[0] ]){
			unsigned int local_i = global_to_local_index_map[local_dof_indices[0] ];
	        for (unsigned int i=0; i<GeometryInfo<3>::vertices_per_cell; ++i){

				ver_cell_iter adjacent_cell = vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].begin();						
				ver_cell_iter end_cell = vertex_to_cell_iterator[ vertex_map[cell->vertex_index(i)] ].end();

				for(;adjacent_cell != end_cell; ++adjacent_cell) {
					cell_all_neighbor_iterator[local_i].insert(*adjacent_cell);
				}
			}
			cell_all_neighbor_iterator[local_i].erase(cell);
		}

	}


	vertex_index.clear();
	face_index.clear();
	vertex_to_cell_iterator.clear();
	vertex_to_cell_iterator.clear();

	cell_neighbor_index.resize(n_relevant_cells);
	cell_diagonal_neighbor_index.resize(n_relevant_cells);
	cell_neighbor_neighbor_index.resize(n_relevant_cells);
	cell_all_neighbor_index.resize(n_post_process_cell);

    for (unsigned int i = 0; i < n_relevant_cells; i++) {
		cell_neighbor_index[i].resize(6);
		cell_neighbor_neighbor_index[i].resize(6);
    }

	std::set<DoFHandler<3>::active_cell_iterator>::iterator iter;
	typedef typename std::set<DoFHandler<3>::active_cell_iterator>::iterator dia_cell_iter;
	dia_cell_iter diagonal_cell, end_cell, adjacent_cell;

    for (unsigned int i = 0; i < n_relevant_cells; ++i) {

	    for (unsigned int f=0; f<GeometryInfo<3>::faces_per_cell; ++f) {

			if (cell_neighbor_neighbor_iterator[i][f].size() > 0) {
				cell_neighbor_neighbor_index[i][f].resize(cell_neighbor_neighbor_iterator[i][f].size() );
				if(cell_neighbor_neighbor_iterator[i][f].size() > 1) {
        			std::cerr << "Error: Multiple neighbor of neighbor are found!"<<std::endl;
		    	    std::exit(EXIT_FAILURE);
				}
				iter = cell_neighbor_neighbor_iterator[i][f].begin();
				cell = *iter;
				cell->get_dof_indices(local_dof_indices);
				cell_neighbor_neighbor_index[i][f][0] = local_dof_indices[0];				
			}

			diagonal_cell = cell_neighbor_iterator[i][f].begin();
			end_cell = cell_neighbor_iterator[i][f].end();
			cell_neighbor_index[i][f].resize(cell_neighbor_iterator[i][f].size() );

			for (unsigned int d = 0; diagonal_cell != end_cell; ++diagonal_cell, ++d) {			
				cell = *diagonal_cell;
				cell->get_dof_indices(local_dof_indices);
				cell_neighbor_index[i][f][d] = local_dof_indices[0];	
			}
		}

		diagonal_cell = cell_diagonal_neighbor_iterator[i].begin();
		end_cell = cell_diagonal_neighbor_iterator[i].end();
		cell_diagonal_neighbor_index[i].resize(cell_diagonal_neighbor_iterator[i].size() );

		for (unsigned int d = 0; diagonal_cell != end_cell; ++diagonal_cell, ++d) {			
			cell = *diagonal_cell;
			cell->get_dof_indices(local_dof_indices);
			cell_diagonal_neighbor_index[i][d] = local_dof_indices[0];	
		}

    }

    for (unsigned int i = 0; i < n_post_process_cell; ++i) {
		adjacent_cell = cell_all_neighbor_iterator[i].begin();
		end_cell = cell_all_neighbor_iterator[i].end();
		cell_all_neighbor_index[i].resize(cell_all_neighbor_iterator[i].size() );
		for (unsigned int n = 0; adjacent_cell != end_cell; ++adjacent_cell, ++n) {	
			cell = *adjacent_cell;
			cell->get_dof_indices(local_dof_indices);
			cell_all_neighbor_index[i][n] = local_dof_indices[0];	
		}
	}
//	std::cout<<"1"<<std::endl;

    for (unsigned int i = 0; i < cell_all_neighbor_index.size(); ++i) {
	    for (unsigned int d = 0; d < cell_all_neighbor_index[i].size(); ++d) {
			locally_relevant_dofs.add_index(cell_all_neighbor_index[i][d] );
			is_ghost_cell[cell_all_neighbor_index[i][d]] = true;
		}

	}

    for (unsigned int i = 0; i < cell_neighbor_neighbor_index.size(); ++i) {
	    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; ++f) {
			if(cell_neighbor_neighbor_index[i][f].size() > 0) {
				locally_relevant_dofs.add_index(cell_neighbor_neighbor_index[i][f][0] );
				is_ghost_cell[cell_neighbor_neighbor_index[i][f][0]] = true;
			}
		}
	}
//	std::cout<<"2"<<std::endl;
	cell = dof_handler.begin_active();
    for (; cell != endc; ++cell) {

        cell->get_dof_indices(local_dof_indices);

		if(is_ghost_cell[local_dof_indices[0] ] && !is_relevant_cell[local_dof_indices[0] ] && !is_post_processing_cell[local_dof_indices[0] ]) {

			local_index_to_iterator.push_back(cell);
			global_to_local_index_map[local_dof_indices[0] ] = local_index;
			local_to_global_index_map.push_back(local_dof_indices[0] );
			local_index++;
		}

	}

	n_store_cell = local_index;

	unsigned int data_cell = locally_relevant_dofs.n_elements();
	unsigned int ghost_cell = data_cell - n_locally_cells;

	pcout<<"data cell: "<<Utilities::MPI::max (data_cell, MPI_COMM_WORLD)<<std::endl;

	pcout<<"max no of ghost cell: "<<Utilities::MPI::max (ghost_cell, MPI_COMM_WORLD)<<std::endl;
	pcout<<"min no of ghost cell: "<<Utilities::MPI::min (ghost_cell, MPI_COMM_WORLD)<<std::endl;

	double ratio_ghost = static_cast<double>(n_locally_cells)/static_cast<double>(ghost_cell);
	double min_ratio = Utilities::MPI::min (ratio_ghost, MPI_COMM_WORLD);
	double max_ratio = Utilities::MPI::max (ratio_ghost, MPI_COMM_WORLD);

	pcout<<"min ration of locally owned cell to ghost cell: "<<min_ratio<<std::endl
		<<"min ration of locally owned cell to ghost cell: "<<max_ratio<<std::endl;

/*	
	Utilities::System::MemoryStats stat;
	Utilities::System::get_memory_stats(stat);
	pcout<<"Setup system memory: "<<std::endl;
	pcout<<stat.VmRSS/std::pow(2,20)<<std::endl;
	pcout<<"total: "<<24.0*stat.VmRSS/std::pow(2,20)<<std::endl;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "log.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Memory consumption in setup system per node = " << stat.VmRSS/std::pow(2,20) << std::endl;
        fout_convergence.close();
	}
*/
	cell_all_neighbor_iterator.clear();
	cell_diagonal_neighbor_iterator.clear();
	cell_neighbor_iterator.clear();
	is_extra_ghost_cell.clear();
/*
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "timer.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Time taken to setup the system = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}
*/
}

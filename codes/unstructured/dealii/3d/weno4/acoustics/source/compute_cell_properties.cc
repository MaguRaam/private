#include "../include/Weno432.h" 

void Weno4_3D::compute_cell_properties() {

    // Iterate over all the cells 
//    auto start = std::chrono::system_clock::now();
    
    DoFHandler<3>::active_cell_iterator
    cell = dof_handler.begin_active(),
    endc = dof_handler.end();
/*
    // For obtaining cell quadrature points and jxws
    QGauss<3> quadrature_formula(2);
    FEValues<3> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);
    
    // For obtaining face normal vectors
	QGauss<3-1> face_quadrature_formula(2);
	FEFaceValues<3> fv_face_values (fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors | update_JxW_values);

    // For obtaining face normal vectors at face center
	QGauss<3-1> face_quadrature_formula_center(1);
	FEFaceValues<3> fv_face_values_center (fv, face_quadrature_formula_center, update_quadrature_points | update_normal_vectors);
*/

    // For obtaining cell quadrature points and jxws
    QGauss<3> quadrature_formula(2);
    FEValues<3> fv_values (mapping, fv, quadrature_formula, update_quadrature_points | update_JxW_values);
    
    // For obtaining face normal vectors
	QGauss<3-1> face_quadrature_formula(2);
	FEFaceValues<3> fv_face_values (mapping, fv, face_quadrature_formula, update_quadrature_points | update_normal_vectors | update_JxW_values);

    // For obtaining face normal vectors at face center
	QGauss<3-1> face_quadrature_formula_center(1);
	FEFaceValues<3> fv_face_values_center (mapping, fv, face_quadrature_formula_center, update_quadrature_points | update_normal_vectors);

    double volume, volume_min_local = 1e6, volume_min;
    std::vector< Point<3> > cell_quadrature_points(8);
    std::vector< double > jxws(8);
    std::vector< std::vector<Point<3> > > face_quadrature_point(6);
	std::vector < std::vector<double> > face_jxw(6);
	std::vector < std::vector<double> > face_normal_x(6);
	std::vector < std::vector<double> > face_normal_y(6);
	std::vector < std::vector<double> > face_normal_z(6);
    std::vector < Point<3> > face_center_quadrature_point(6);
    std::vector <double> face_center_normal_x(6);
    std::vector <double> face_center_normal_y(6);
    std::vector <double> face_center_normal_z(6);
    std::vector <double> surface_area(6);

	Tensor<1,3> face_normal_vector, face_center_normal_vector;

    for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++) {
		face_quadrature_point[f].resize(4);
		face_jxw[f].resize(4);
		face_normal_x[f].resize(4);
		face_normal_y[f].resize(4);
		face_normal_z[f].resize(4);
	}
    
	for (; cell != endc; ++cell) {

        cell->get_dof_indices(local_dof_indices);

		if (is_ghost_cell[local_dof_indices[0] ] ){
           
	        volume = 0.0; 

			fv_values.reinit(cell);
			for (unsigned int i=0; i<fv_values.n_quadrature_points; ++i){
				volume += fv_values.JxW (i);
	            cell_quadrature_points[i] = fv_values.quadrature_point(i);
	            jxws[i] =  fv_values.JxW (i);            
	        }
			if (volume < volume_min_local) volume_min_local = volume;
	        for (unsigned int f = 0; f < GeometryInfo<3>::faces_per_cell; f++) {
                
				fv_face_values.reinit(cell, f);
    	        fv_face_values_center.reinit(cell, f);                
				surface_area[f] = 0.0;        
				for (unsigned int i=0; i<fv_face_values.n_quadrature_points; ++i){    	            

	    	        face_normal_vector = fv_face_values.normal_vector(i);
	    	        face_normal_x[f][i] = face_normal_vector[0]; 
	    	        face_normal_y[f][i] = face_normal_vector[1]; 
	    	        face_normal_z[f][i] = face_normal_vector[2]; 
    	            
	    	        face_quadrature_point[f][i] = fv_face_values.quadrature_point(i);
					face_jxw[f][i] = fv_face_values.JxW (i);
		            surface_area[f] += fv_face_values.JxW (i);
				}

	    	    face_center_normal_vector = fv_face_values_center.normal_vector(0);     	           
    	        face_center_quadrature_point[f] = fv_face_values_center.quadrature_point(0); 
    	        face_center_normal_x[f] = face_center_normal_vector[0]; 
    	        face_center_normal_y[f] = face_center_normal_vector[1]; 
    	        face_center_normal_z[f] = face_center_normal_vector[2]; 
	
			}
				
			unsigned int c = global_to_local_index_map[local_dof_indices[0] ];
            Cell[c].reinit(volume, cell_quadrature_points, jxws, 
					face_quadrature_point, face_jxw, face_normal_x, face_normal_y, face_normal_z,
					face_center_quadrature_point, face_center_normal_x, face_center_normal_y, face_center_normal_z, surface_area); 
		}
	}

	volume_min = Utilities::MPI::min (volume_min_local, MPI_COMM_WORLD);
	h_min = std::pow(volume_min, (1.0/3.0) );
	pcout<<"h_min: "<<h_min<<"\ncompute cell properties"<<std::endl;
/*
    auto end = std::chrono::system_clock::now();
    std::chrono::duration<double> elapsed_seconds = end-start;

    if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0){
	    std::ofstream fout_convergence ;
        fout_convergence.flags( std::ios::dec | std::ios::scientific ) ;
        fout_convergence.precision(7) ;

        const std::string filename = "timer.dat";
	    fout_convergence.open(filename,std::ios::in | std::ios::out | std::ios::app);

    	fout_convergence << "Time taken to compute cell properties = " << elapsed_seconds.count() << std::endl;
        fout_convergence.close();
	}
*/
}

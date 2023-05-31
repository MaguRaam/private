#include "../include/Weno32.h"
// Output the results in data files 

Vector<double> solve_system(FullMatrix<double> A, Vector<double> b) {
    
    Vector<double> x(2);
    
    double det = A(0,0)*A(1,1) - A(0,1)*A(1,0);
    
    assert(det != 0.0);
    
    x(0) = (b(0)*A(1,1) - A(0,1)*b(1))/det;
    x(1) = (A(0,0)*b(1) - b(0)*A(1,0))/det; 
    
    return x; 
}

void Weno3_2D::output_results() {

	DataOut<2> data_out;
	data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (RHO, "RHO",DataOut<2>::type_dof_data);

	data_out.build_patches ();

	std::ostringstream time_;
	time_ << std::setprecision(4);
	time_ << time;

	const std::string filename = ("plots/sol_rho-" +
                                time_.str() +
                                "." +
                                Utilities::int_to_string
                                (triangulation.locally_owned_subdomain(), 4));
	std::ofstream output ((filename + ".vtu").c_str());
	data_out.write_vtu (output);

	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
    {
      std::vector<std::string> filenames;
      for (unsigned int i=0;
           i<Utilities::MPI::n_mpi_processes(mpi_communicator);
           ++i)
        filenames.push_back ("plots/sol_rho-" +
                             time_.str() +
                             "." +
                             Utilities::int_to_string (i, 4) +
                             ".vtu");

      std::ofstream master_output ((filename + ".pvtu").c_str());
      data_out.write_pvtu_record (master_output, filenames);
    }
}


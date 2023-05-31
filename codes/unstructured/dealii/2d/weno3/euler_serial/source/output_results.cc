#include "../include/Weno32.h"

// Output the results in data files 

void Weno3_2D::output_results(unsigned int i) {
	
	DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);
    data_out.add_data_vector (RHO, "RHO", DataOut<2>::type_dof_data);
    data_out.build_patches ();
    const std::string filename = "plots/plot_" + Utilities::int_to_string (i, 4) + ".dat";
    std::ofstream output (filename.c_str());
    data_out.write_tecplot(output);

} 

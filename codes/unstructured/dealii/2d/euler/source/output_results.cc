#include "../include/Euler.h"

// Output the results in data files 

void Euler_2D::output_results(unsigned int i) {
    
    DataOut<2> data_out1;
    data_out1.attach_dof_handler (dof_handler);
    data_out1.add_data_vector (RHO, "RHO", DataOut<2>::type_dof_data);
    data_out1.build_patches ();
    const std::string filename1 = "plots/plot_RHO_" + Utilities::int_to_string (i, 5) + ".dat";
    std::ofstream output1 (filename1.c_str());
    data_out1.write_tecplot(output1);
    /*
    DataOut<2> data_out2;
    data_out2.attach_dof_handler (dof_handler);
    data_out2.add_data_vector (RHO_V, "RHO_V", DataOut<2>::type_dof_data);
    data_out2.build_patches ();
    const std::string filename2 = "plots/plot_RHO_V_" + Utilities::int_to_string (i, 5) + ".dat";
    std::ofstream output2 (filename2.c_str());
    data_out2.write_tecplot(output2);
	*/ 
} 

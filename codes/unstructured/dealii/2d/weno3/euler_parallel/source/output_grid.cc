#include "../include/Weno32.h"

// Output the results in data files 

void Weno3_2D::output_grid() {

    DataOut<2> data_out;
    data_out.attach_dof_handler (dof_handler);

    data_out.build_patches ();

    const std::string filename = "plots/grid_"+
                                  Utilities::int_to_string (triangulation.locally_owned_subdomain(), 4) + ".vtu";
    std::ofstream output (filename.c_str());
    data_out.write_vtu (output);

}

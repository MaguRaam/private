#include "../include/Weno432.h"

// Initialize the solution

void Weno4_3D::initialize()
{

	unsigned int N_gp = 3; // No. of quadrature points
	QGauss<3> quadrature_formula(N_gp);
	//    FEValues<3> fv_values (fv, quadrature_formula, update_quadrature_points | update_JxW_values);
	FEValues<3> fv_values(mapping, fv, quadrature_formula, update_quadrature_points | update_JxW_values);

	DoFHandler<3>::active_cell_iterator cell = dof_handler.begin_active();
	DoFHandler<3>::active_cell_iterator endc = dof_handler.end();

	Point<3> q_point;

	double V0;

	unsigned int g_i;

	Vector<double> U_(2);
	double u, v;

	Point<3> q_point;
	cell = dof_handler.begin_active();
	for (unsigned int c = 0; cell != endc; ++cell, ++c)
	{

		if (cell->is_locally_owned())
		{

			cell->get_dof_indices(local_dof_indices);
			g_i = local_dof_indices[0];

			V0 = 0.0;
			fv_values.reinit(cell);
			for (unsigned int i = 0; i < fv_values.n_quadrature_points; ++i)
				V0 += fv_values.JxW(i);

			u = 0.0, v = 0.0;

			for (unsigned int i = 0; i < fv_values.n_quadrature_points; i++)
			{

				q_point = fv_values.quadrature_point(i);

				U_ = initial_condition(q_point, h_min);

				u += (1. / V0) * fv_values.JxW(i) * U_(0);
				v += (1. / V0) * fv_values.JxW(i) * U_(1);
			}

			local_U(g_i) = u;
			local_V(g_i) = v;
		}
	}

	local_U.compress(VectorOperation::insert);
	local_V.compress(VectorOperation::insert);

	U = local_U;
	V = local_V;

	//	std::cout<<"rank: "<<Utilities::MPI::this_mpi_process(mpi_communicator)<<" RHO: "<<RHO(465653)<<std::endl;

	if (Utilities::MPI::this_mpi_process(mpi_communicator) == 0)
	{
		std::ofstream fout_convergence;
		fout_convergence.flags(std::ios::dec | std::ios::scientific);
		fout_convergence.precision(7);

		const std::string filename = "log.dat";
		fout_convergence.open(filename, std::ios::in | std::ios::out | std::ios::app);

		fout_convergence << "h_min: " << h_min << std::endl;
		fout_convergence.close();
	}
}

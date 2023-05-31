/*
 * polynomials.cc
 *      Author: sunder
 */

#include "polynomials.hh"

Legendre_m_2D::Legendre_m_2D()
{
	order = 0;
	inv_mass_matrix = new double[0];
	dofs = 0;
}


Legendre_m_2D::Legendre_m_2D(unsigned int N)
{

	order = N;

	inv_mass_matrix = new double[10];

	dofs = N*(N+1)/2;



	// Initialize the inverse mass matrices

	inv_mass_matrix[0] = 1.0;    inv_mass_matrix[1] = 12.0;
	inv_mass_matrix[2] = 12.0;   inv_mass_matrix[3] = 180.0;
	inv_mass_matrix[4] = 144.0;  inv_mass_matrix[5] = 180.0;
	inv_mass_matrix[6] = 2800.0; inv_mass_matrix[7] = 2160.0;
	inv_mass_matrix[8] = 2160.0; inv_mass_matrix[9] = 2800.0;
}

void Legendre_m_2D::initialize(unsigned int N) {

	order = N;
	dofs = N*(N+1)/2;


	delete[] inv_mass_matrix;


	inv_mass_matrix = new double[10];

	// Initialize the inverse mass matrices

	inv_mass_matrix[0] = 1.0;    inv_mass_matrix[1] = 12.0;
	inv_mass_matrix[2] = 12.0;   inv_mass_matrix[3] = 180.0;
	inv_mass_matrix[4] = 144.0;  inv_mass_matrix[5] = 180.0;
	inv_mass_matrix[6] = 2800.0; inv_mass_matrix[7] = 2160.0;
	inv_mass_matrix[8] = 2160.0; inv_mass_matrix[9] = 2800.0;

}

Legendre_m_2D::~Legendre_m_2D() {
	//delete[] inv_mass_matrix;
}

double Legendre_m_2D::evaluate_basis_function(const unsigned int& index, const double& x, const double& y) {

	static const double r1_12 = 1./12.;
	static const double r3_20 = 3./20.;

	double retval = 0.0;

	switch (index) {

	case 0:
		retval = 1.0;
		break;
	case 1:
		retval = x;
		break;
	case 2:
		retval = y;
		break;
	case 3:
		retval = x*x - r1_12;
		break;
	case 4:
		retval = y*y - r1_12;
		break;
	case 5:
		retval = x*y;
		break;
	case 6:
		retval = x*(x*x - r3_20);
		break;
	case 7:
		retval = y*(y*y - r3_20);
		break;
	case 8:
		retval = y*(x*x - r1_12);
		break;
	case 9:
		retval = x*(y*y - r1_12);
		break;
	default:
		Assert(index < 10, ErrNotImplemented());
		break;
	}

	return retval;
}


double Legendre_m_2D::evaluate(const std::vector<double>& coeffs, const double& x, const double& y) {

	double retval = 0.0;

	for (unsigned int i = 0; i < dofs; ++i) {
		retval += coeffs[i]*evaluate_basis_function(i, x, y);
	}

	return retval;
}

void Legendre_m_2D::evaluate_basis_function_grad(const unsigned int& index, const double& x, const double& y, double* grad) {

	static const double r1_12 = 1./12.;
	static const double r3_20 = 3./20.;

	switch (index) {

	case 0:
		grad[0] = 0.0;
		grad[1] = 0.0;
		break;
	case 1:
		grad[0] = 1.0;
		grad[1] = 0.0;
		break;
	case 2:
		grad[0] = 0.0;
		grad[1] = 1.0;
		break;
	case 3:
		grad[0] = 2.0*x;
		grad[1] = 0.0;
		break;
	case 4:
		grad[0] = 0.0;
		grad[1] = 2.0*y;
		break;
	case 5:
		grad[0] = y;
		grad[1] = x;
		break;
	case 6:
		grad[0] = 3.0*x*x - r3_20;
		grad[1] = 0.0;
		break;
	case 7:
		grad[0] = 0.0;
		grad[1] = 3.*y*y - r3_20;
		break;
	case 8:
		grad[0] = 2.0*y*x;
		grad[1] = (x*x - r1_12);
		break;
	case 9:
		grad[0] = (y*y - r1_12);
		grad[1] = 2.0*x*y;
		break;
	default:
		Assert(index < 10, ErrNotImplemented());
		break;
	}
}

void Legendre_m_2D::evaluate_grad(const std::vector<double>& coeffs, const double& x, const double& y, double* grad) {

	grad[0] = 0.0; grad[1] = 0.0;

	double basis_grad[2];

	for (unsigned int i = 0; i < dofs; ++i) {

		evaluate_basis_function_grad(i, x, y, basis_grad);

		grad[0] += coeffs[i]*basis_grad[0];

		grad[1] += coeffs[i]*basis_grad[1];
	}
}

unsigned int Legendre_m_2D::no_dofs() const {
	return dofs;
}

double Legendre_m_2D::Inverse_Mass_Matrix(unsigned int i) const {

	Assert(i < 10, ErrNotImplemented());

	return inv_mass_matrix[i];
}




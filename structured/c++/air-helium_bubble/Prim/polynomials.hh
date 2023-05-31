/*
 * polynomials.hh
 *      Author: sunder
 */

#ifndef POLYNOMIALS_HH_
#define POLYNOMIALS_HH_

#include "headers.hh"

// Monic Legendre Polynomials in 1D

class Legendre_m_1D;
class Legendre_m_2D;

class Legendre_m_2D {

	unsigned int order;
	unsigned int dofs;
	double* inv_mass_matrix;


public:
	Legendre_m_2D();
	Legendre_m_2D(unsigned int);
	void initialize(unsigned int);
	~Legendre_m_2D();

	double evaluate_basis_function(const unsigned int&, const double&, const double&);
	double evaluate(const std::vector<double>&, const double&, const double&);
	void evaluate_basis_function_grad(const unsigned int&, const double&, const double&, double*);
	void evaluate_grad(const std::vector<double>&, const double&, const double&, double*);
	unsigned int no_dofs() const;

	double Inverse_Mass_Matrix(unsigned int) const;

};




#endif /* POLYNOMIALS_HH_ */

/*
 * la.cc
 *      Author: sunder
 */

#include "la.hh"

/* Vector class */

//----------------------------------------------------------------------------
// Default constructor
//----------------------------------------------------------------------------

Vector::Vector() : n_(0), vec(NULL) {}

//----------------------------------------------------------------------------
// Constructor taking size of the vector (all elements are initialized to zero
//----------------------------------------------------------------------------

Vector::Vector(int size) : n_(size), vec(n_>0 ? new double[n_] : NULL)
{
	int i;
	if (n_ > 0 ) {
		for (i = 0; i < n_; ++i)
			vec[i] = 0.0;
	}
}

//----------------------------------------------------------------------------
// Copy constructor
//----------------------------------------------------------------------------

Vector::Vector(const Vector &rhs) : n_(rhs.n_), vec(n_>0 ? new double[n_] : NULL)
{
	for(int i=0; i<n_; ++i) vec[i] = rhs.vec[i];
}

//----------------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------------

Vector::~Vector()
{
	if (vec != NULL) delete[] vec;
}

//----------------------------------------------------------------------------
//  Reinitialize to the given size (values are not preserved)
//----------------------------------------------------------------------------

void Vector::reinit(int nn)
{
	int i;

	if (nn != n_) {
		if (vec != NULL) delete[] (vec);
		n_ = nn;
		vec = n_ > 0 ? new double[n_] : NULL;

		if (n_ > 0 ) {
			for (i = 0; i < n_; ++i)
				vec[i] = 0.0;
		}

	}
}

//----------------------------------------------------------------------------
//  Subscript operators
//----------------------------------------------------------------------------

double& Vector::operator[](const int i)
{
#ifdef BOUNDS_CHECK
	if (i < 0 || i >= n_) {
		std::cerr << "Vector index i = " << i << " out of range [0, " << n_ << ")" << std::endl;
		std::exit(EXIT_FAILURE);
	}
#endif
	return vec[i];
}

const double& Vector::operator[](const int i) const
{
#ifdef BOUNDS_CHECK
	if (i < 0 || i >= n_) {
		std::cerr << "Vector index i = " << i << " out of range [0, " << n_ << ")" << std::endl;
		std::exit(EXIT_FAILURE);
	}
#endif
	return vec[i];
}

double& Vector::operator()(const int i)
{
#ifdef BOUNDS_CHECK
	if (i < 0 || i >= n_) {
		std::cerr << "Vector index i = " << i << " out of range [0, " << n_ << ")" << std::endl;
		std::exit(EXIT_FAILURE);
	}
#endif
	return vec[i];
}


const double& Vector::operator()(const int i) const
{
#ifdef BOUNDS_CHECK
	if (i < 0 || i >= n_) {
		std::cerr << "Vector index i = " << i << " out of range [0, " << n_ << ")" << std::endl;
		std::exit(EXIT_FAILURE);
	}
#endif
	return vec[i];
}

/* Matrix class */

//----------------------------------------------------------------------------
//  Default constructor
//----------------------------------------------------------------------------

Matrix::Matrix() : m(0), n(0), nelem(0), mat(NULL) {}

//----------------------------------------------------------------------------
//  Constructor taking size of the array. Initialize all elements to zero
//----------------------------------------------------------------------------

Matrix::Matrix(int mm, int nn) : m(mm), n(nn), nelem(mm*nn), mat(nelem>0 ? new double[nelem] : NULL)
{
	if (nelem > 0) {
		for (int i = 0; i < nelem; ++i)
			mat[i] = 0.0;
	}
}

//----------------------------------------------------------------------------
//  Copy constructor
//----------------------------------------------------------------------------

Matrix::Matrix(const Matrix &rhs) : m(rhs.m), n(rhs.n), nelem(rhs.nelem), mat(nelem>0 ? new double[nelem] : NULL)
{
	for(int i=0; i<nelem; i++) mat[i] = rhs.mat[i];
}

//----------------------------------------------------------------------------
// Destructor
//----------------------------------------------------------------------------

Matrix::~Matrix()
{
	if (mat != NULL) delete[] mat;
}

//----------------------------------------------------------------------------
//  Reinitialize to the given size (values are not preserved)
//----------------------------------------------------------------------------

void Matrix::reinit(int mm, int nn)
{
	if ( mm != m || nn != n) {
		if (mat != NULL) delete[] (mat);
		m = mm;
		n = nn;
		nelem = m*n;
		mat = nelem > 0 ? new double[nelem] : NULL;

		if (nelem > 0) {
			for (int i = 0; i < nelem; ++i)
				mat[i] = 0.0;
		}
	}
}

//----------------------------------------------------------------------------
//  Subscript operators
//----------------------------------------------------------------------------

double& Matrix::operator()(const int i, const int j)
{
#ifdef BOUNDS_CHECK
	if (i < 0 || i >= m) {
		std::cerr << "Matrix index i = " << i << " out of range [0, " << m << ")" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if (j < 0 || j >= n) {
		std::cerr << "Matrix index j = " << j << " out of range [0, " << n << ")" << std::endl;
		std::exit(EXIT_FAILURE);
	}
#endif
	return mat[i*n+j];
}

const double& Matrix::operator()(const int i, const int j) const
{
#ifdef BOUNDS_CHECK
	if (i < 0 || i >= m) {
		std::cerr << "Matrix index i = " << i << " out of range [0, " << m << ")" << std::endl;
		std::exit(EXIT_FAILURE);
	}
	if (j < 0 || j >= n) {
		std::cerr << "Matrix index j = " << j << " out of range [0, " << n << ")" << std::endl;
		std::exit(EXIT_FAILURE);
	}
#endif
	return mat[i*n+j];
}

/* QR decomposition */

//----------------------------------------------------------------------------
// Default constructor (Does nothing)
//----------------------------------------------------------------------------

QRdcmp::QRdcmp()
{}

//----------------------------------------------------------------------------
// Constructor for QRdcmp class. The input is a rectangular (m x n, m >= n)
// matrix
//----------------------------------------------------------------------------

QRdcmp::QRdcmp(const Matrix& A)
{
	initialize(A);
}

void QRdcmp::initialize(const Matrix& A) {

	const int m = A.rows();
	const int n = A.cols();

	double s, nrm;
	int i = 0, j = 0, k = 0;

	Rdiag.reinit(n);

	// Copy A into QR_

	QR_.reinit(m, n);

	for (i = 0; i < m; ++i) {
		for (j = 0; j < n; ++j) {
			QR_(i,j) = A(i,j);
		}
	}

	// Main loop.

	for (k = 0; k < n; k++) {

		// Compute 2-norm of k-th column without under/overflow.

		nrm = 0.0;
		for (i = k; i < m; i++) {
			nrm = std::hypot(nrm,QR_(i,k));
		}

		if (nrm != 0.0) {
		  // Form k-th Householder vector.
		  if (QR_(k,k) < 0.0) {
			 nrm = -nrm;
		  }
		  for (i = k; i < m; i++) {
			 QR_(i,k) /= nrm;
		  }
		  QR_(k,k) += 1.0;

		  // Apply transformation to remaining columns.
		  for (j = k+1; j < n; j++) {
			 s = 0.0;
			 for (i = k; i < m; i++) {
				s += QR_(i,k)*QR_(i,j);
			 }
			 s = -s/QR_(k,k);
			 for (i = k; i < m; i++) {
				QR_(i,j) += s*QR_(i,k);
			 }
		  }
	   }

		Rdiag(k) = -nrm;
	}

}

//----------------------------------------------------------------------------
// Check whether the given system is a full rank system
//----------------------------------------------------------------------------

bool QRdcmp::is_full_rank() const {

	for (int j = 0; j < QR_.cols(); j++) {
       if (std::abs(Rdiag[j]) < 1.0e-12)
          return false;
    }

    return true;
}

//----------------------------------------------------------------------------
// Get the R (n x n) matrix
//----------------------------------------------------------------------------

void QRdcmp::get_R(Matrix& R) const {

	for (int i = 0; i < QR_.cols(); i++) {
		for (int j = 0; j < QR_.cols(); j++) {
			if (i < j) {
				R(i,j) = QR_(i,j);
			} else if (i == j) {
				R(i,j) = Rdiag(i);
			} else {
				R(i,j) = 0.0;
			}
       }
    }
}

//----------------------------------------------------------------------------
// Get the Q (m x n) matrix
//----------------------------------------------------------------------------

void QRdcmp::get_Q(Matrix& Q) const {

	  int i = 0, j = 0, k = 0;
	  double s;

	  for (k = QR_.cols()-1; k >= 0; k--) {

		  for (i = 0; i < QR_.rows(); i++)
			  Q(i,k) = 0.0;

    	Q(k,k) = 1.0;

    	for (j = k; j < QR_.cols(); j++) {
    		if (QR_(k,k) != 0.0) {

    			s = 0.0;

    			for (i = k; i < QR_.rows(); i++)
    				s += QR_(i,k)*Q(i,j);

    			s = -s/QR_(k,k);

    			for (i = k; i < QR_.rows(); i++)
    				Q(i,j) += s*QR_(i,k);

    		}
    	}
    }
}

//----------------------------------------------------------------------------
// Solve least squares system A*x = b. b is m length vector and x is n length
// vector
//----------------------------------------------------------------------------

void QRdcmp::solve(const Vector& b,  Vector& x_) const {

	const int m = QR_.rows();
	const int n = QR_.cols();

	double* x = new double[m];

	for (int i = 0; i < m; ++i)
		x[i] = b[i];

	double s;

    // Compute Y = transpose(Q)*b

	for (int k = 0; k < n; k++) {
		s = 0.0;
		for (int i = k; i < m; i++)
			s += QR_(i,k)*x[i];

		s = -s/QR_(k,k);

		for (int i = k; i < m; i++)
			x[i] += s*QR_(i,k);
    }

    // Solve R*X = Y;

    for (int k = n-1; k >= 0; k--) {

    	x[k] /= Rdiag[k];

       for (int i = 0; i < k; i++)
             x[i] -= x[k]*QR_(i,k);

    }

    // return n x nx portion of X

    for (int i=0; i<n; i++)
    	x_[i] = x[i];

    delete[] x;
}


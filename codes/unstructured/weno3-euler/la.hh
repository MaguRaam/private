/*
 * la.hh
 *      Author: sunder
 */

#ifndef LA_HH_
#define LA_HH_

#include "headers.hh"

/* Functions related to linear algebra */

// Vector class

class Vector {
	int n_;
	double* vec;
public:
	Vector();
    explicit Vector(int);
    Vector(const Vector&);
    ~Vector();
    void reinit(int);
	double & operator[](const int);
	const double & operator[](const int) const;
	double & operator()(const int);
	const double & operator()(const int) const;

	int n() const {return n_;}

};

// Matrix class

class Matrix {
	int m;
	int n;
	int nelem;
    double* mat;
public:
    Matrix();
    explicit Matrix(int,int);
    Matrix(const Matrix&);
    ~Matrix();
    void reinit(int,int);
	double& operator()(const int, const int);
	const double& operator()(const int, const int) const;
	int rows() const {return m;}
	int cols() const {return n;}
};

// QR decomposition for solving least-squares problem

class QRdcmp {
	Matrix QR_;
	Vector Rdiag;

public:
    QRdcmp();
    QRdcmp(const Matrix&);
    void initialize(const Matrix&);
    bool is_full_rank() const;
    void get_R(Matrix&) const;
    void get_Q(Matrix&) const;
    void solve(const Vector&, Vector&) const;
};



#endif /* LA_HH_ */

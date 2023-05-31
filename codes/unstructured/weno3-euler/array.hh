/*
 * array.hh
 *      Author: sunder
 */

#ifndef ARRAY_HH_
#define ARRAY_HH_

#include "headers.hh"

/* Simple arrays classes to hold data */

template<typename T>
void check_index_range(T i, T max) {
	if (i < 0 || i >= max) {
	    std::cerr << "index " << i << " out of range [0, " << max << ")"  << std::endl;
	    std::exit(EXIT_FAILURE);
	}
}

//-----------------------------------------------------------------------
// 1D Array
//-----------------------------------------------------------------------

template <class T>
class Array1D {
	int size_;
    T* data;
public:
    typedef T value_type;
    Array1D();
    explicit Array1D(int);
    Array1D(const Array1D&);
    ~Array1D();
    void reinit(int);
	inline T & operator[](const int);
	inline const T & operator[](const int) const;
	inline T & operator()(const int);
	inline const T & operator()(const int) const;
	int size() const {return size_;}
};

/* Default constructor */

template <class T>
Array1D<T>::Array1D() : size_(0), data(NULL) {}

/* Constructor taking size of the array */

template <class T>
Array1D<T>::Array1D(int n) : size_(n), data(n>0 ? new T[n] : NULL) {}

/* Copy constructor */

template <class T>
Array1D<T>::Array1D(const Array1D<T> &rhs) : size_(rhs.size_), data(size_>0 ? new T[size_] : NULL)
{
	for(int i=0; i<size_; i++) data[i] = rhs.data[i];
}

/* Destructor */

template <class T>
Array1D<T>::~Array1D()
{
	if (data != NULL) delete[] data;
}

/* Reinitialize to the given size */

template <class T>
void Array1D<T>::reinit(int newsize)
{
	if (newsize != size_) {
		if (data != NULL) delete[] (data);
		size_ = newsize;
		data = size_ > 0 ? new T[size_] : NULL;
	}
}

/* Subscript operators */

template <class T>
inline T & Array1D<T>::operator[](const int i)
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size_);
#endif
	return data[i];
}

template <class T>
inline const T& Array1D<T>::operator[](const int i) const
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size_);
#endif
	return data[i];
}

template <class T>
inline T & Array1D<T>::operator()(const int i)
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size_);
#endif
	return data[i];
}

template <class T>
inline const T& Array1D<T>::operator()(const int i) const
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size_);
#endif
	return data[i];
}

//-----------------------------------------------------------------------
// 2D Array
//-----------------------------------------------------------------------

template <class T>
class Array2D {
	int size1_;
	int size2_;
	int nelem;
    T* data;
public:
    typedef T value_type;
    Array2D();
    explicit Array2D(int,int);
    Array2D(const Array2D&);
    ~Array2D();
    void reinit(int,int);
	inline T & operator()(const int, const int);
	inline const T & operator()(const int, const int) const;
	int size1() const {return size1_;}
	int size2() const {return size2_;}
};

/* Default constructor */

template <class T>
Array2D<T>::Array2D() : size1_(0), size2_(0), nelem(0), data(NULL) {}

/* Constructor taking size of the array */

template <class T>
Array2D<T>::Array2D(int m, int n) : size1_(m), size2_(n), nelem(m*n), data(nelem>0 ? new T[nelem] : NULL) {}

/* Copy constructor */

template <class T>
Array2D<T>::Array2D(const Array2D<T> &rhs) : size1_(rhs.size1_), size2_(rhs.size2_), nelem(rhs.nelem), data(nelem>0 ? new T[nelem] : NULL)
{
	for(int i=0; i<nelem; i++) data[i] = rhs.data[i];
}

/* Destructor */

template <class T>
Array2D<T>::~Array2D()
{
	if (data != NULL) delete[] data;
}

/* Reinitialize to the given size */

template <class T>
void Array2D<T>::reinit(int nsize1_, int nsize2_)
{
	if (nsize1_ != size1_ || nsize2_ != size2_) {
		if (data != NULL) delete[] (data);
		size1_ = nsize1_;
		size2_ = nsize2_;
		nelem = size1_*size2_;
		data = nelem > 0 ? new T[nelem] : NULL;
	}
}

/* Subscript operators */

template <class T>
inline T & Array2D<T>::operator()(const int i, const int j)
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size1_);
	check_index_range(j,size2_);
#endif
	return data[i*size2_+j];
}

template <class T>
inline const T& Array2D<T>::operator()(const int i, const int j) const
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size1_);
	check_index_range(j,size2_);
#endif
	return data[i*size2_+j];
}

//-----------------------------------------------------------------------
// 3D Array
//-----------------------------------------------------------------------

template <class T>
class Array3D {
	int size1_;
	int size2_;
	int size3_;
	int C1,C2,C3; // Storage order
	int nelem;
    T* data;
public:
    typedef T value_type;
    Array3D();
    explicit Array3D(int,int,int);
    Array3D(const Array3D&);
    ~Array3D();
    void reinit(int,int,int);
	inline T & operator()(const int,const int,const int);
	inline const T & operator()(const int,const int,const int) const;
	int size1() const {return size1_;}
	int size2() const {return size2_;}
	int size3() const {return size3_;}
};

/* Default constructor */

template <class T>
Array3D<T>::Array3D() : size1_(0), size2_(0), size3_(0), C1(0), C2(0), C3(0), nelem(0), data(NULL) {}

/* Constructor taking size of the array */

template <class T>
Array3D<T>::Array3D(int dim1, int dim2, int dim3) :
	size1_(dim1),
	size2_(dim2),
	size3_(dim3),
	C1(dim2*dim3),
	C2(dim3),
	C3(1),
	nelem(dim1*dim2*dim3),
	data(nelem>0 ? new T[nelem] : NULL)
{}

/* Copy constructor */

template <class T>
Array3D<T>::Array3D(const Array3D<T> &rhs) :
	size1_(rhs.size1_),
	size2_(rhs.size2_),
	size3_(rhs.size3_),
	C1(rhs.C1),
	C2(rhs.C2),
	C3(rhs.C3),
	nelem(rhs.nelem),
	data(nelem>0 ? new T[nelem] : NULL)
{
	for(int i=0; i<nelem; i++) data[i] = rhs.data[i];
}

/* Destructor */

template <class T>
Array3D<T>::~Array3D()
{
	if (data != NULL) delete[] data;
}

/* Reinitialize to the given size */

template <class T>
void Array3D<T>::reinit(int dim1, int dim2, int dim3)
{
	if (dim1 != size1_ || dim2 != size2_ || dim3 != size3_) {
		if (data != NULL) delete[] (data);
		size1_ = dim1;
		size2_ = dim2;
		size3_ = dim3;
		C1 = dim2*dim3;
		C2 = dim3;
		C3 = 1;
		nelem = size1_*size2_*size3_;
		data = nelem > 0 ? new T[nelem] : NULL;
	}
}

/* Subscript operators */

template <class T>
inline T & Array3D<T>::operator()(const int i, const int j, const int k)
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size1_);
	check_index_range(j,size2_);
	check_index_range(k,size3_);
#endif
	return data[i*C1 + j*C2 + k];
}

template <class T>
inline const T& Array3D<T>::operator()(const int i, const int j, const int k) const
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size1_);
	check_index_range(j,size2_);
	check_index_range(k,size3_);
#endif
	return data[i*C1 + j*C2 + k];
}

//-----------------------------------------------------------------------
// 4D Array
//-----------------------------------------------------------------------

template <class T>
class Array4D {
	int size1_;
	int size2_;
	int size3_;
	int size4_;
	int C1,C2,C3,C4; // Storage order
	int nelem;
    T* data;
public:
    typedef T value_type;
    Array4D();
    explicit Array4D(int,int,int,int);
    Array4D(const Array4D&);
    ~Array4D();
    void reinit(int,int,int,int);
	inline T & operator()(const int,const int,const int,const int);
	inline const T & operator()(const int,const int,const int,const int) const;
	int size1() const {return size1_;}
	int size2() const {return size2_;}
	int size3() const {return size3_;}
	int size4() const {return size4_;}
};

/* Default constructor */

template <class T>
Array4D<T>::Array4D() : size1_(0), size2_(0), size3_(0), size4_(0), C1(0), C2(0), C3(0), C4(0), nelem(0), data(NULL) {}

/* Constructor taking size of the array */

template <class T>
Array4D<T>::Array4D(int dim1, int dim2, int dim3, int dim4) :
	size1_(dim1),
	size2_(dim2),
	size3_(dim3),
	size4_(dim4),
	C1(dim2*dim3*dim4),
	C2(dim3*dim4),
	C3(dim4),
	C4(1),
	nelem(dim1*dim2*dim3*dim4),
	data(nelem>0 ? new T[nelem] : NULL)
{}

/* Copy constructor */

template <class T>
Array4D<T>::Array4D(const Array4D<T> &rhs) :
	size1_(rhs.size1_),
	size2_(rhs.size2_),
	size3_(rhs.size3_),
	size4_(rhs.size4_),
	C1(rhs.C1),
	C2(rhs.C2),
	C3(rhs.C3),
	C4(rhs.C4),
	nelem(rhs.nelem),
	data(nelem>0 ? new T[nelem] : NULL)
{
	for(int i=0; i<nelem; i++) data[i] = rhs.data[i];
}

/* Destructor */

template <class T>
Array4D<T>::~Array4D()
{
	if (data != NULL) delete[] data;
}

/* Reinitialize to the given size */

template <class T>
void Array4D<T>::reinit(int dim1, int dim2, int dim3, int dim4)
{
	if (dim1 != size1_ || dim2 != size2_ || dim3 != size3_) {
		if (data != NULL) delete[] (data);
		size1_ = dim1;
		size2_ = dim2;
		size3_ = dim3;
		size4_ = dim4;
		C1 = dim2*dim3*dim4;
		C2 = dim3*dim4;
		C3 = dim4;
		C4 = 1;
		nelem = size1_*size2_*size3_*size4_;
		data = nelem > 0 ? new T[nelem] : NULL;
	}
}

/* Subscript operators */

template <class T>
inline T & Array4D<T>::operator()(const int i, const int j, const int k, const int l)
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size1_);
	check_index_range(j,size2_);
	check_index_range(k,size3_);
	check_index_range(l,size4_);
#endif
	return data[i*C1+j*C2+k*C3+l];
}

template <class T>
inline const T& Array4D<T>::operator()(const int i, const int j, const int k, const int l) const
{
#ifdef BOUNDS_CHECK
	check_index_range(i,size1_);
	check_index_range(j,size2_);
	check_index_range(k,size3_);
	check_index_range(l,size4_);
#endif
	return data[i*C1+j*C2+k*C3+l];
}

#endif /* ARRAY_HH_ */

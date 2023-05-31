/*
 * headers.hh
 *      Author: sunder
 */

#ifndef HEADERS_HH_
#define HEADERS_HH_

//----------------------------------------------------------------------------
// Commonly used C/C++ header files
//----------------------------------------------------------------------------

#include <iostream>
#include <string>
#include <cmath>
#include <iomanip>
#include <cassert>
#include <ctime>
#include <sstream>
#include <vector>
#include <algorithm>
#include <chrono>
#include <fstream>
#include <stdio.h>
#include <stdlib.h>
#include <set>
#include <array>
#include <map>
#include <utility>

//#define BOUNDS_CHECK

const int nVar = 4; // Number of variables in the PDE system

/* Physical constants */

const double GAMMA = 1.4;         // Ratio of specific heats
const double prs_floor = 1.0e-12; // Pressure floor value
const double rho_floor = 1.0e-14; // Density floor value

//----------------------------------------------------------------------------
// Memory allocation and deallocation for arrays
//----------------------------------------------------------------------------

// Allocate and free memory for dynamic arrays

template<typename T>
void allocate_mem_1d(T*&, size_t);

template<typename T>
void allocate_mem_2d(T**&, size_t, size_t);

template<typename T>
void allocate_mem_3d(T***&, size_t, size_t, size_t);

template<typename T>
void allocate_mem_4d(T****&, size_t, size_t, size_t, size_t);

template<typename T>
void free_mem_1d(T*&);

template<typename T>
void free_mem_2d(T**&, size_t);

template<typename T>
void free_mem_3d(T**&, size_t, size_t);

template<typename T>
void free_mem_4d(T****&, size_t, size_t, size_t);

// Allocate memory for 1D array

template<typename T>
void allocate_mem_1d(T*& U, size_t m) {
	U = new T[m];
}


// Allocate memory for 2D array

template<typename T>
void allocate_mem_2d(T**& U, size_t m, size_t n) {
	U = new T*[m];

	for (unsigned int i = 0; i < m; ++i) {
		U[i] = new T[n];
	}
}

// Allocate memory for 3D array

template<typename T>
void allocate_mem_3d(T***& U, size_t m, size_t n, size_t o) {
	U = new T**[m];

	for(unsigned int i = 0 ; i < m; ++i) {
		U[i] = new T*[n];
	    for(unsigned int j = 0; j < n; ++j){
	    	U[i][j] = new T[o];
	    }
	}
}

// Allocate memory for 4D array

template<typename T>
void allocate_mem_4d(T****& U, size_t size1, size_t size2, size_t size3, size_t size4) {

	U = new T***[size1];

	for(size_t i = 0 ; i < size1; ++i) {

		U[i] = new T**[size2];

		for(size_t j = 0; j < size2; ++j){

			U[i][j] = new T*[size3];

			for (size_t k = 0; k < size3; ++k) {

				U[i][j][k] = new T[size4];

			}

		}
	}
}

// Release memory for 1D array

template<typename T>
void free_mem_1d(T*& U) {
	delete[] U;
}

// Release memory for 2D array

template<typename T>
void free_mem_2d(T**& U, size_t m) {

	for (unsigned int i = 0; i < m; ++i) {
		delete[] U[i];
	}

	delete[] U;
}

// Release memory for 3D array

template<typename T>
void free_mem_3d(T***& U, size_t m, size_t n) {
    for (unsigned int i = 0; i < m; ++i) {
        for (unsigned int j = 0; j < n; ++j) {
            delete [] U[i][j];
        }
        delete [] U[i];
    }
    delete [] U;
}

// Release memory for 4D array

template<typename T>
void free_mem_4d(T****& U, size_t size1, size_t size2, size_t size3) {

	for (size_t i = 0; i < size1; ++i) {
		for (size_t j = 0; j < size2; ++j) {
			for (unsigned int k = 0; k < size3; ++k) {
				delete[] U[i][j][k];
			}
			delete[] U[i][j];
		}
		delete[] U[i];
	}

	delete[] U;
}


#endif /* HEADERS_HH_ */


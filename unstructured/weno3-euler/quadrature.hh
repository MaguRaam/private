/*
 * quadrature.hh
 *      Author: sunder
 */

#ifndef QUADRATURE_HH_
#define QUADRATURE_HH_

#include "array.hh"
#include "headers.hh"

/* Quadrature points and weights on triangles and straight lines*/

// 1D
void QGauss(double,double,Array1D<double>&,Array1D<double>&,int);

// 2D (On triangles)
void QCentroid(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QGauss4x4(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QGauss8x8(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QSevenPoint(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang1(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang2(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang3(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang4(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang5(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang6(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang7(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang8(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QStrang9(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QToms584_19(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QToms612_19(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QSToms612_28(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QToms706_37(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QVertex(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);
void QGaussTriangle(Array1D<double>&, Array1D<double>&, Array1D<double>&, int&);


#endif /* QUADRATURE_HH_ */

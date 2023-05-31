#ifndef CELL_PROPERTIES_H_
#define CELL_PROPERTIES_H_
 
#include "Headers.h"

// Class to store properties of a cell 

class cell_properties; 

class cell_properties {
    double area; 
    Point<2>* quadrature_points1; 
    Point<2>* quadrature_points2; 
    double* face_normal_x1; 
    double* face_normal_x2;
    double* face_normal_y1;
    double* face_normal_y2;
    double* surface_area; 
public:
    // Constructors and destructors 
    cell_properties(); 
    cell_properties(double, Point<2>*, Point<2>*, double*, double*, double*, double*, double*); 
    cell_properties (const cell_properties &);
    ~cell_properties(); 
    
    // Reinitialize 
    void reinit(double, Point<2>*, Point<2>*, double*, double*, double*, double*, double*); 
    
    // Properties of cell
    double measure() const; 
    Point<2> face_quadrature_point1(unsigned int) const; 
    Point<2> face_quadrature_point2(unsigned int) const;
    double nx1(unsigned int) const;
    double nx2(unsigned int) const;
    double ny1(unsigned int) const;
    double ny2(unsigned int) const;
    double S_f(unsigned int) const; 
};
 
 
#endif /* CELL_PROPERTIES_H_ */

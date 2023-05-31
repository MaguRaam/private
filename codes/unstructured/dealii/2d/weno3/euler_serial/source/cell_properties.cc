#include "../include/cell_properties.h"

// Class member function definitions 

// Default constructor 

cell_properties::cell_properties() {

    area = 0.0; 
    quadrature_points1 = new Point<2>[0];
    quadrature_points2 = new Point<2>[0];
    face_normal_x1 = new double[0]; 
    face_normal_x2 = new double[0];
    face_normal_y1 = new double[0];
    face_normal_y2 = new double[0];
    surface_area = new double[0]; 
}


// Main constructor taking arguments 

cell_properties::cell_properties(double _area, Point<2> _quadrature_points1[], Point<2> _quadrature_points2[], double* _face_normal_x1, double* _face_normal_x2, double* _face_normal_y1, double* _face_normal_y2, double* _surface_area) {

    area = _area; 
    
    // Initialize all the arrays 
    
    quadrature_points1 = new Point<2>[4];
    quadrature_points2 = new Point<2>[4];
    face_normal_x1 = new double[4]; 
    face_normal_x2 = new double[4];
    face_normal_y1 = new double[4];
    face_normal_y2 = new double[4];
    surface_area = new double[4]; 
    
    // Fill all the arrays 
    
    for (unsigned int f = 0; f < 4; f++) {
    
        quadrature_points1[f] = _quadrature_points1[f];
        quadrature_points2[f] = _quadrature_points2[f];
        face_normal_x1[f]     = _face_normal_x1[f]; 
        face_normal_x2[f]     = _face_normal_x2[f];
        face_normal_y1[f]     = _face_normal_y1[f];
        face_normal_y2[f]     = _face_normal_y2[f];
        surface_area[f]       = _surface_area[f]; 
    }

}

// Copy constructor 

cell_properties::cell_properties(const cell_properties& cp) {

    area = cp.area; 
    
        // Initialize all the arrays 
    
    quadrature_points1 = new Point<2>[4];
    quadrature_points2 = new Point<2>[4];
    face_normal_x1 = new double[4]; 
    face_normal_x2 = new double[4];
    face_normal_y1 = new double[4];
    face_normal_y2 = new double[4];
    surface_area = new double[4]; 
    
    // Fill all the arrays 
    
    for (unsigned int f = 0; f < 4; f++) {
    
        quadrature_points1[f] = cp.quadrature_points1[f];
        quadrature_points2[f] = cp.quadrature_points2[f];
        face_normal_x1[f]     = cp.face_normal_x1[f]; 
        face_normal_x2[f]     = cp.face_normal_x2[f];
        face_normal_y1[f]     = cp.face_normal_y1[f];
        face_normal_y2[f]     = cp.face_normal_y2[f];
        surface_area[f]       = cp.surface_area[f]; 
    }
    
    
}

// Destructor 

cell_properties::~cell_properties() {

    delete[] quadrature_points1;
    delete[] quadrature_points2;
    delete[] face_normal_x1; 
    delete[] face_normal_x2;
    delete[] face_normal_y1;
    delete[] face_normal_y2;
    delete[] surface_area; 
}

// Reinitialize the object 

void cell_properties::reinit(double _area, Point<2>* _quadrature_points1, Point<2>* _quadrature_points2, double* _face_normal_x1, double* _face_normal_x2, double* _face_normal_y1, double* _face_normal_y2, double* _surface_area) {
    
    // Delete the old data 
    
    delete[] quadrature_points1;
    delete[] quadrature_points2;
    delete[] face_normal_x1; 
    delete[] face_normal_x2;
    delete[] face_normal_y1;
    delete[] face_normal_y2;
    delete[] surface_area; 
    
    // Create the data again 
    
    area = _area; 
    
    // Initialize all the arrays 
    
    quadrature_points1 = new Point<2>[4];
    quadrature_points2 = new Point<2>[4];
    face_normal_x1 = new double[4]; 
    face_normal_x2 = new double[4];
    face_normal_y1 = new double[4];
    face_normal_y2 = new double[4];
    surface_area = new double[4]; 
    
    // Fill all the arrays 
    
    for (unsigned int f = 0; f < 4; f++) {
    
        quadrature_points1[f] = _quadrature_points1[f];
        quadrature_points2[f] = _quadrature_points2[f];
        face_normal_x1[f]     = _face_normal_x1[f]; 
        face_normal_x2[f]     = _face_normal_x2[f];
        face_normal_y1[f]     = _face_normal_y1[f];
        face_normal_y2[f]     = _face_normal_y2[f];
        surface_area[f]       = _surface_area[f]; 
    }

}


// Get the propeties of the cell 

// 1. area of the cell 

double cell_properties::measure() const {
    return area;
}

// 2. First quadrature point on face f 

Point<2> cell_properties::face_quadrature_point1(unsigned int face_index) const {

    assert(face_index < 4);
    return quadrature_points1[face_index]; 
}

// 2. Second quadrature point on face f 

Point<2> cell_properties::face_quadrature_point2(unsigned int face_index) const {

    assert(face_index < 4);
    return quadrature_points2[face_index]; 
}

// 3. First normal (x) on face 1 

double cell_properties::nx1(unsigned int face_index) const {

    assert(face_index < 4);
    return face_normal_x1[face_index]; 
}

// 4. Second normal (x) on face 1 

double cell_properties::nx2(unsigned int face_index) const {

    assert(face_index < 4);
    return face_normal_x2[face_index]; 
}

// 5. First normal (y) on face 1 

double cell_properties::ny1(unsigned int face_index) const {

    assert(face_index < 4);
    return face_normal_y1[face_index]; 
}

// 5. Second normal (y) on face 1 

double cell_properties::ny2(unsigned int face_index) const {

    assert(face_index < 4);
    return face_normal_y2[face_index]; 
}

// 6. Area of a face f 

double cell_properties::S_f(unsigned int face_index) const {

    assert(face_index < 4);
    return surface_area[face_index]; 
}
 

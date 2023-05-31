#ifndef CELL_PROPERTIES_H_
#define CELL_PROPERTIES_H_
 
#include "Headers.h"

// Class to store properties of a cell 

class cell_properties; 

class cell_properties {
    double volume;
    std::vector < Point<3> > cell_quadrature_points;     
    std::vector <double> jxws;                      
    std::vector < std::vector< Point<3> > > face_quadrature_points;
	std::vector < std::vector<double> > face_jxws;
	std::vector < std::vector<double> > face_normal_x;
	std::vector < std::vector<double> > face_normal_y;
	std::vector < std::vector<double> > face_normal_z;
    std::vector < Point<3> > face_center_quadrature_points;
    std::vector	< double > face_center_normal_x;
    std::vector	< double > face_center_normal_y;
    std::vector	< double > face_center_normal_z;
    std::vector <double> surface_area;
	double h_size; 
 
public:
    // Constructors and destructors 
    cell_properties(); 
    cell_properties(double, const std::vector < Point<3> >&, const std::vector <double>&, 
			const std::vector < std::vector< Point<3> > >&, const std::vector < std::vector<double> >&, const std::vector < std::vector<double> >&, const std::vector < std::vector<double> >&, const std::vector < std::vector<double> >&, 
			const std::vector < Point<3> >&, const std::vector< double >&, const std::vector< double >&, const std::vector<double >&, const std::vector <double>&);
    cell_properties (const cell_properties &);
    ~cell_properties(); 
    
    // Reinitialize 
    void reinit(double, const std::vector < Point<3> >&, const std::vector <double>&, 
			const std::vector < std::vector< Point<3> > >&, const std::vector < std::vector<double> >&, const std::vector < std::vector<double> >&, const std::vector < std::vector<double> >&, const std::vector < std::vector<double> >&, 
			const std::vector < Point<3> >&, const std::vector< double >&, const std::vector< double >&, const std::vector<double >&, const std::vector <double>&);
    // Properties of cell
    double measure() const;
    Point<3> cell_quadrature_point(unsigned int) const;
    double jxw(unsigned int) const;                        
    Point<3> face_quadrature_point(unsigned int, unsigned int) const;
	double face_jxw(unsigned int, unsigned int) const;
	double nx(unsigned int, unsigned int) const;
	double ny(unsigned int, unsigned int) const;
	double nz(unsigned int, unsigned int) const;
    Point<3> face_center_quadrature_point(unsigned int) const;
    double center_nx(unsigned int) const;
    double center_ny(unsigned int) const; 
    double center_nz(unsigned int) const; 
    double S_f(unsigned int) const;
	double h() const;

};
 
 
#endif /* CELL_PROPERTIES_H_ */

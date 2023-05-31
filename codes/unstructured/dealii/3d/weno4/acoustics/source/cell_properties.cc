#include "../include/cell_properties.h"

// Class member function definitions 

// Default constructor 

cell_properties::cell_properties() {
    volume = 0.0;
	h_size = 0.0;
	// check for std vector constructor
}


// Main constructor taking arguments 

    std::vector < Point<3> > cell_quadrature_points;     
    std::vector <double> jxws;                      
    std::vector < std::vector< Point<3> > > face_quadrature_points;
	std::vector < std::vector<double> > face_jxws;
	std::vector < std::vector<double> > face_normal_x;
	std::vector < std::vector<double> > face_normal_y;
	std::vector < std::vector<double> > face_normal_z;
    std::vector < Point<3> > face_center_quadrature_points;
    std::vector< double > face_center_normal_x;
    std::vector< double > face_center_normal_y;
    std::vector< double > face_center_normal_z;
    std::vector <double> surface_area;

cell_properties::cell_properties(double _volume, const std::vector < Point<3> >& _cell_quadrature_points, const std::vector <double>& _jxws, 
			const std::vector < std::vector< Point<3> > >& _face_quadrature_points, const std::vector < std::vector<double> >& _face_jxws, 
			const std::vector < std::vector<double> >& _face_normal_x, const std::vector < std::vector<double> >& _face_normal_y, 
			const std::vector < std::vector<double> >& _face_normal_z, const std::vector < Point<3> >& _face_center_quadrature_points, 
			const std::vector< double >& _face_center_normal_x, const std::vector< double >& _face_center_normal_y, 
			const std::vector< double >& _face_center_normal_z, const std::vector<double>& _surface_area):

	cell_quadrature_points(_cell_quadrature_points),
	jxws(_jxws),
	face_quadrature_points(_face_quadrature_points),
	face_jxws(_face_jxws),
	face_normal_x(_face_normal_x),
	face_normal_y(_face_normal_y),
	face_normal_z(_face_normal_z),
	face_center_quadrature_points(_face_center_quadrature_points),
	face_center_normal_x(_face_center_normal_x),
	face_center_normal_y(_face_center_normal_y),
	face_center_normal_z(_face_center_normal_z),
	surface_area(_surface_area)
{
    volume = _volume; 
	h_size = std::pow(volume, (1./3.) );
}

// Copy constructor 

cell_properties::cell_properties(const cell_properties& cp) {

    volume = cp.volume; 
	cell_quadrature_points = cp.cell_quadrature_points;
	jxws = cp.jxws;
    face_quadrature_points = cp.face_quadrature_points;
	face_jxws = cp.face_jxws;
	face_normal_x = cp.face_normal_x;
	face_normal_y = cp.face_normal_y;
	face_normal_z = cp.face_normal_z;
	face_center_quadrature_points = cp.face_center_quadrature_points;
	face_center_normal_x = cp.face_center_normal_x;
	face_center_normal_y = cp.face_center_normal_y;
	face_center_normal_z = cp.face_center_normal_z;
	surface_area = cp.surface_area;
	h_size = std::pow(volume, (1./3.) );
    
}

// Destructor 

cell_properties::~cell_properties() {
 
}

// Reinitialize the object 

void cell_properties::reinit(double _volume, const std::vector < Point<3> >& _cell_quadrature_points, const std::vector <double>& _jxws, 
			const std::vector < std::vector< Point<3> > >& _face_quadrature_points, const std::vector < std::vector<double> >& _face_jxws, 
			const std::vector < std::vector<double> >& _face_normal_x, const std::vector < std::vector<double> >& _face_normal_y, 
			const std::vector < std::vector<double> >& _face_normal_z, const std::vector < Point<3> >& _face_center_quadrature_points, 
			const std::vector< double >& _face_center_normal_x, const std::vector< double >& _face_center_normal_y, 
			const std::vector< double >& _face_center_normal_z, const std::vector<double>& _surface_area) {
        
    // Create the data again 
    
    volume = _volume;
	cell_quadrature_points = _cell_quadrature_points;
	jxws = _jxws;
    face_quadrature_points = _face_quadrature_points;
	face_jxws = _face_jxws;
	face_normal_x = _face_normal_x;
	face_normal_y = _face_normal_y;
	face_normal_z = _face_normal_z;
	face_center_quadrature_points = _face_center_quadrature_points;
	face_center_normal_x = _face_center_normal_x;
	face_center_normal_y = _face_center_normal_y;
	face_center_normal_z = _face_center_normal_z;
	surface_area = _surface_area;
	h_size = std::pow(volume, (1./3.) );

}


// Get the propeties of the cell 

// 1. area of the cell 

double cell_properties::measure() const {
    return volume;
}

// 2. cell quadrature point

Point<3> cell_properties::cell_quadrature_point(unsigned int quad_point) const {

    assert(quad_point < 8);
    return cell_quadrature_points[quad_point]; 
}

// 3. cell jxw at quadrature point

double cell_properties::jxw(unsigned int quad_point) const {

    assert(quad_point < 8);
    return jxws[quad_point]; 
}

// 4. jxw at qth quadrature point on face f 

double cell_properties::face_jxw(unsigned int face_index, unsigned int quad_index) const {

    assert(face_index < 6);
    assert(quad_index < 4);
    return face_jxws[face_index][quad_index]; 
}

// 5. nx at qth quadrature point on face f 

double cell_properties::nx(unsigned int face_index, unsigned int quad_index) const {

    assert(face_index < 6);
    assert(quad_index < 4);
    return face_normal_x[face_index][quad_index]; 
}

// 6. nx at qth quadrature point on face f 

double cell_properties::ny(unsigned int face_index, unsigned int quad_index) const {

    assert(face_index < 6);
    assert(quad_index < 4);
    return face_normal_y[face_index][quad_index]; 
}

// 7. nx at qth quadrature point on face f 

double cell_properties::nz(unsigned int face_index, unsigned int quad_index) const {

    assert(face_index < 6);
    assert(quad_index < 4);
    return face_normal_z[face_index][quad_index];
}

// 8. qth quadrature point on face f 

Point<3> cell_properties::face_quadrature_point(unsigned int face_index, unsigned int quad_index) const {

    assert(face_index < 6);
    assert(quad_index < 4);
    return face_quadrature_points[face_index][ quad_index]; 
}

// 9. quadrature point at center of face f 

Point<3> cell_properties::face_center_quadrature_point(unsigned int face_index) const {

    assert(face_index < 6);
    return face_center_quadrature_points[face_index]; 
}

// 10. x componenet of normal at center of face f 

double cell_properties::center_nx(unsigned int face_index) const {

    assert(face_index < 6);
    return face_center_normal_x[face_index]; 
}

// 11. y componenet of normal at center of face f 

double cell_properties::center_ny(unsigned int face_index) const {

    assert(face_index < 6);
    return face_center_normal_y[face_index]; 
}

// 12. z componenet of normal at center of face f 

double cell_properties::center_nz(unsigned int face_index) const {

    assert(face_index < 6);
    return face_center_normal_z[face_index]; 
}

// 13. Area of a face f 

double cell_properties::S_f(unsigned int face_index) const {

    assert(face_index < 6);
    return surface_area[face_index]; 
}

// 14. Size of the cell

double cell_properties::h() const {

	return h_size; 
}

/*
 * mappings.hh
 *      Author: sunder
 */

#ifndef MAPPINGS_HH_
#define MAPPINGS_HH_

#include "headers.hh"

/* Functions related to mappings and transformations */

void reference_to_physical(double, double, double, double, double, double, double, double, double&, double&);
void physical_to_reference(double, double, double, double, double, double, double, double, double&, double&);
void straight_face_to_inclined_face(double, double&, double&);


double basis(double,double,int);


double basis_average_over_triangle(double, double, double, double, double, double,int);


#endif /* MAPPINGS_HH_ */

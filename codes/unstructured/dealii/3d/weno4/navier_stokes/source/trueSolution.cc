#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "../include/trueSolution.h"

Point<3> grid_transform (const Point<3> &in)
{

	double k_x = 0.20, k_y = 0.21, k_z = 0.23, f = 2.0/3.0;
	Point<3> in_star_star(in(0), in(1), in(2) + k_x*std::sin(2*in(0)/f) + k_y*std::sin(2*in(1)/f));
	Point<3> in_star(in_star_star(0), in_star_star(1) + k_x*std::sin(2*in_star_star(0)/f) + k_z*std::sin(2*in_star_star(2)/f), in_star_star(2));
	Point<3> out(in_star(0) + k_y*std::sin(2*in_star(1)/f) + k_z*std::sin(2*in_star(2)/f), in_star(1), in_star(2));
	return out;
}

Point<3>
SineTransform::push_forward(const Point<3> &chart_point) const
{
	double k_x = 0.20, k_y = 0.21, k_z = 0.23, f = 2.0/3.0;
	Point<3> in_star_star(chart_point(0), chart_point(1), chart_point(2) + k_x*std::sin(2*chart_point(0)/f) + k_y*std::sin(2*chart_point(1)/f));
	Point<3> in_star(in_star_star(0), in_star_star(1) + k_x*std::sin(2*in_star_star(0)/f) + k_z*std::sin(2*in_star_star(2)/f), in_star_star(2));
	Point<3> out(in_star(0) + k_y*std::sin(2*in_star(1)/f) + k_z*std::sin(2*in_star(2)/f), in_star(1), in_star(2));
	return out;
}

Point<3>
SineTransform::pull_back(const Point<3> &space_point) const
{

	double k_x = 0.20, k_y = 0.21, k_z = 0.23, f = 2.0/3.0;

	Point<3> in_star_star(space_point(0) - k_y*std::sin(2*space_point(1)/f) - k_z*std::sin(2*space_point(2)/f), space_point(1), space_point(2));
	Point<3> in_star(in_star_star(0), in_star_star(1) - k_x*std::sin(2*in_star_star(0)/f) - k_z*std::sin(2*in_star_star(2)/f), in_star_star(2));
	Point<3> out(in_star(0), in_star(1), in_star(2) - k_x*std::sin(2*in_star(0)/f) - k_y*std::sin(2*in_star(1)/f));
	return out;
}


std::unique_ptr<Manifold<3,3> >
SineTransform::clone() const
{
  return std_cxx14::make_unique<SineTransform>();
}

void
TrueSolution::
soln(const Point<3> &p, Vector<double> &values, const double &t) const
{   
	double x = p(0), y = p(1), z = p(2);

	values(0) = 1.0;
	values(1) = sin(x*2.0*M_PI)*sin(y*2.0*M_PI)*sin(z*2.0*M_PI);
	values(2) = sin(x*2.0*M_PI)*sin(y*2.0*M_PI)*sin(z*2.0*M_PI);
	values(3) = sin(x*2.0*M_PI)*sin(y*2.0*M_PI)*sin(z*2.0*M_PI);
	values(4) = 5.0;  

}


void 
TrueSolution::
soln_list(const std::vector<Point<3> > &points,
             std::vector< Vector<double> > &values, const double &t) const{
  // check whether component is in the valid range is up to the derived
  // class
  Assert (values.size() == points.size(),
          ExcDimensionMismatch(values.size(), points.size()));

  for (unsigned int i=0; i<points.size(); ++i)
    this->soln (points[i], values[i],t);
}

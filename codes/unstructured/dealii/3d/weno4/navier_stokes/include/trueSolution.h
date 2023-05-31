#ifndef TRUESOLUTION_H_
#define TRUESOLUTION_H_

#include <deal.II/base/function.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/lac/vector.h>
#include <cmath>
#include <fstream>
#include <iostream>
#include "Headers.h"

Point<3> grid_transform (const Point<3> &in);

class SineTransform : public ChartManifold<3,3>
{
public:
	virtual
  Point<3>
  pull_back(const Point<3> &space_point) const override;
  virtual
  Point<3>
  push_forward(const Point<3> &chart_point) const override;

  virtual std::unique_ptr<Manifold<3,3> > clone() const override;

};

class TrueSolution : public Function<3>
{
public:
  TrueSolution() : Function<3>(1)
  {}

  void soln(const Point<3> &p, Vector<double> &values, const double &t) const;

  void soln_list(const std::vector<Point<3> > &points,
                             std::vector< Vector<double> > &values, const double &t) const;

};

#endif

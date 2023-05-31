/*
 * weno.hh
 *      Author: sunder
 */

#ifndef WENO_HH_
#define WENO_HH_

#include "headers.hh"

void weno_2d(const std::vector<double>&,
		   const std::vector<double>&,
		   const std::vector<double>&,
		   std::vector<double>&, int);



#endif /* WENO_HH_ */

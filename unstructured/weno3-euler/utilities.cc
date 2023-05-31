/*
 * utilities.cc
 *      Author: sunder
 */

#include "utilities.hh"

//-----------------------------------------------------------------------
// Denote a discontinuity using smooth tanh function.
// f(x) = a for x < x0, f(x) = b otherwise. This function smoothens the
// discontinuity by eps amount
//-----------------------------------------------------------------------

double tanh_profile(double a, double b, double x0, double eps, double x) {
	return (0.5*(a+b) + 0.5*(b-a)*tanh((x-x0)/eps));
}

//-----------------------------------------------------------------------
// Convert integer to string along with the required padding
//-----------------------------------------------------------------------

std::string int_to_string (unsigned int value, const unsigned int digits) {
    std::string lc_string = std::to_string(value);

    if (lc_string.size() < digits) {
        // We have to add the padding zeroes in front of the number
        const unsigned int padding_position = (lc_string[0] == '-')
                                                ?
                                                1
                                                :
                                                0;

        const std::string padding(digits - lc_string.size(), '0');
        lc_string.insert(padding_position, padding);
    }

    return lc_string;
}

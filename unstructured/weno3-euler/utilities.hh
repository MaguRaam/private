/*
 * utilities.hh
 *      Author: sunder
 */

#ifndef UTILITIES_HH_
#define UTILITIES_HH_

#include "headers.hh"

/* Small frequently used functions */

double tanh_profile(double, double, double, double, double);
std::string int_to_string (unsigned int, const unsigned int);


//--------------------------------------------------------------------------------------
// Templated and inline functions

// Check if two floating point numbers are close

template<class T>
bool isclose(const T& a, const T& b, const double rtol=1.0e-05,
		const double atol=1.0e-08) {

	if (std::isnan(a) || std::isnan(b) )
		return false;

	else
		return std::abs(a-b) <= (atol + rtol * std::abs(b)) ? true : false;
}

//
// Small integer powers with minimum number of multiplications.
//
// @see https://en.wikipedia.org/wiki/Addition-chain_exponentiation
//

inline double   pow1(const double& a) {                                                                        return a;               }
inline double   pow2(const double& a) {                                                                        return a * a;           }
inline double   pow3(const double& a) {                                                                        return a * a * a;       }
inline double   pow4(const double& a) { const double a2=a*a;                                                   return a2 * a2;         }
inline double   pow5(const double& a) { const double a2=a*a;                                                   return a2 * a2 * a;     }
inline double   pow6(const double& a) { const double a2=a*a;                                                   return a2 * a2 * a2;    }
inline double   pow7(const double& a) { const double a2=a*a;                                                   return a2 * a2 * a2 * a;}
inline double   pow8(const double& a) { const double a2=a*a;   const double a4=a2*a2;                          return a4 * a4;         }
inline double   pow9(const double& a) { const double a3=a*a*a;                                                 return a3 * a3 * a3;    }
inline double  pow10(const double& a) { const double a2=a*a;   const double a4=a2*a2;                          return a4 * a4 * a2;    }
inline double  pow11(const double& a) { const double a2=a*a;   const double a4=a2*a2;                          return a4*a4*a2*a;      }
inline double  pow12(const double& a) { const double a2=a*a;   const double a4=a2*a2;                          return a4*a4*a4;        }
inline double  pow13(const double& a) { const double a2=a*a;   const double a4=a2*a2;                          return a4*a4*a4*a;      }
inline double  pow14(const double& a) { const double a2=a*a;   const double a4=a2*a2;                          return a4*a4*a4*a2;     }
inline double  pow15(const double& a) { const double a2=a*a;   const double a5=a2*a2*a;                        return a5*a5*a5;        }
inline double  pow16(const double& a) { const double a2=a*a;   const double a4=a2*a2;   const double a8=a4*a4; return a8*a8;           }




#endif /* UTILITIES_HH_ */

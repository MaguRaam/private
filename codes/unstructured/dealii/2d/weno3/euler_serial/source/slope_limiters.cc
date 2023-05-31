#include "../include/Slope_Limiters.h"

// Slope Limiters 

// 1. Superbee Limiter 
double superbee(double r) {
	double min1 = std::min(1.0, 2.0*r);
    double min2 = std::min(2.0, r);
    double max1 = std::max(min1, min2);
    double max = std::max(0.0, max1);
    return max; 
}

// 2. minmod Limiter 
double minmod(double r) {
	return std::max(0.0, std::min(r, 1.0));
}

// 3. Osher Limiter 
double osher(double r) {
    return std::max(0.0, std::min(r, 2.0));
}

// 4. MUSCL Limiter 
double muscl(double r) {
    return (r + std::abs(r))/(1 + std::abs(r));
}

// 5. Monotized-Central Limiter 

double mc(double r) {
    return std::max(0.0, std::min(2.0*r, std::min(0.5*(1.0+r), 2.0))); 
}

// 6. van-Leer Limiter

double van_leer(double r) {

    return ((r + std::abs(r))/(1.0 + std::abs(r))); 
}

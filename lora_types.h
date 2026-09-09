#pragma once
#include <complex>
#include <vector>
#include <cstdint>

namespace lora {

using cd   = std::complex<double>;   // single complex baseband sample (I/Q)
using cvec = std::vector<cd>;

inline double dmod(double a, double b) {
    double r = std::fmod(a, b);
    if (r < 0) r += b;
    return r;
}

inline long imod(long a, long b) {
    long r = a % b;
    if (r < 0) r += b;
    return r;
}

} 

#pragma once

#ifdef NDEBUG
#undef NDEBUG
#endif

#include "namelist.h"

#include <cassert>
#include <cmath>


bool approx(double a, double b) {
    constexpr double epsilon = 1e-14;
    return std::abs(a - b) < epsilon;
}

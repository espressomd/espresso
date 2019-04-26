#ifndef ESPRESSO_ABS_HPP
#define ESPRESSO_ABS_HPP

#include "utils/device_qualifier.hpp"

#ifndef __CUDACC__
#include <cmath>
#endif

namespace Utils {
/**
 * @brief Return the absolute value of x.
 */
inline DEVICE_QUALIFIER double abs(double x) { return fabs(x); }

/**
 * @brief Return the absolute value of x.
 */
inline DEVICE_QUALIFIER float abs(float x) { return fabsf(x); }
} // namespace Utils

#endif // ESPRESSO_DEVICE_MATH_HPP

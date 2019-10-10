#ifndef MATHEVAL_IMPLEMENTATION
#error "Do not include math.hpp directly!"
#endif

#pragma once

#include <boost/math/constants/constants.hpp>

#include <cmath>

namespace matheval {

namespace math {

/// @brief Sign function
template <typename T> T sgn(T x) { return (T{0} < x) - (x < T{0}); }

/// @brief isnan function with adjusted return type
template <typename T> T isnan(T x) { return std::isnan(x); }

/// @brief isinf function with adjusted return type
template <typename T> T isinf(T x) { return std::isinf(x); }

/// @brief Convert radians to degrees
template <typename T> T deg(T x) {
    return x * boost::math::constants::radian<T>();
}

/// @brief Convert degrees to radians
template <typename T> T rad(T x) {
    return x * boost::math::constants::degree<T>();
}

/// @brief unary plus
template <typename T> T plus(T x) { return x; }

/// @brief binary plus
template <typename T> T plus(T x, T y) { return x + y; }

/// @brief unary minus
template <typename T> T minus(T x) { return -x; }

/// @brief binary minus
template <typename T> T minus(T x, T y) { return x - y; }

/// @brief multiply
template <typename T> T multiplies(T x, T y) { return x * y; }

/// @brief divide
template <typename T> T divides(T x, T y) { return x / y; }

} // namespace math

} // namespace matheval

#ifndef UTILS_MATH_GAUSSIAN_HPP
#define UTILS_MATH_GAUSSIAN_HPP

#include <utils/Vector.hpp>

#include <cmath>

namespace Utils {
inline double gaussian(Vector3d x, Vector3d x0, double sigma) {
  return std::exp(-((x - x0).norm2() / (2. * sigma * sigma)));
}

inline Utils::Vector3d del_gaussian(Vector3d x, Vector3d x0, double sigma) {
  return -(x - x0) * gaussian(x, x0, sigma) / (sigma * sigma);
}
} // namespace Utils

#endif // UTILS_MATH_GAUSSIAN_HPP

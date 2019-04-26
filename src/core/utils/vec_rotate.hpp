#ifndef UTILS_VEC_ROTATE_HPP
#define UTILS_VEC_ROTATE_HPP

#include "utils.hpp"
#include "utils/Vector.hpp"
#include "utils/math/sqr.hpp"

#include <cmath>

namespace Utils {
/** rotates vector around axis by alpha */
template <typename T1, typename T2, typename T3>
void vec_rotate(const T1 &axis, double alpha, const T2 &vector, T3 &result) {
  double sina, cosa, absa, a[3];
  sina = std::sin(alpha);
  cosa = std::cos(alpha);
  absa = std::sqrt(scalar(axis, axis));

  a[0] = axis[0] / absa;
  a[1] = axis[1] / absa;
  a[2] = axis[2] / absa;

  result[0] = (cosa + Utils::sqr(a[0]) * (1 - cosa)) * vector[0] +
              (a[0] * a[1] * (1 - cosa) - a[2] * sina) * vector[1] +
              (a[0] * a[2] * (1 - cosa) + a[1] * sina) * vector[2];
  result[1] = (a[0] * a[1] * (1 - cosa) + a[2] * sina) * vector[0] +
              (cosa + Utils::sqr(a[1]) * (1 - cosa)) * vector[1] +
              (a[1] * a[2] * (1 - cosa) - a[0] * sina) * vector[2];
  result[2] = (a[0] * a[2] * (1 - cosa) - a[1] * sina) * vector[0] +
              (a[1] * a[2] * (1 - cosa) + a[0] * sina) * vector[1] +
              (cosa + Utils::sqr(a[2]) * (1 - cosa)) * vector[2];
}

/** rotates vector around axis by alpha */
inline Vector3d vec_rotate(Vector3d axis, double alpha, Vector3d vector) {
  Vector3d result;
  vec_rotate(axis, alpha, vector, result);
  return result;
}

} // namespace Utils

#endif

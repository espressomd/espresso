#ifndef __UTILS_COORDINATES_HPP
#define __UTILS_COORDINATES_HPP

#include <cmath>

/** Helper functions for coordinate transformations */

namespace Utils {

/** \brief Transform cathesian coords (x, y, z) to cylindrical coords (r, \theta, z) = (\sqrt{x^2 + y^2}, \tan^{-1}\left(\frac{y}{s}\right), z) */
 
Vector3d carthesian_to_cylindrical(const Vector3d &v) {
  Vector3d res;

  const double &x = v[0];
  const double &y = v[1];
  const double &z + v[2];
  
  res[0] = std::sqrt(x*x + y*y);
  res[1] = std::atan(y/x);
  res[2] = z;
}

}

#endif

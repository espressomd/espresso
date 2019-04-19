#include "Torus.hpp"
#include "utils/math/sqr.hpp"
#include <cmath>

namespace Shapes {

void Torus::calculate_dist(const Utils::Vector3d &pos, double *dist,
                           double *vec) const {
  /* Coordinate transform to cylinder coords
     with origin at m_center. */
  Utils::Vector3d const c_dist = pos - m_center;
  auto const z = e_z * c_dist;
  auto const r_vec = c_dist - z * e_z;
  auto const r = r_vec.norm();

  *dist = (sqrt(Utils::sqr(r - m_rad) + z * z) - m_tube_rad) * m_direction;
  Utils::Vector3d const dir_vec = c_dist - r_vec * m_rad / r;
  auto const dir_vec_norm = dir_vec / dir_vec.norm();
  vec[0] = dir_vec_norm[0] * m_direction;
  vec[1] = dir_vec_norm[1] * m_direction;
  vec[2] = dir_vec_norm[2] * m_direction;
}
} // namespace Shapes

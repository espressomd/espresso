#include "SimplePore.hpp"

#include <cassert>

namespace Shapes {
/**
 * @brief Calculate the distance function in the coordinates of the cylinder.
 *
 * @param r Distance from the cylinder axis.
 * @param z Distance from the center along the axis.
 *
 * @returns The distance vector from the surface in the cylinder system.
 */
std::pair<double, double> SimplePore::dist_half_pore(double r, double z) const {
  assert(z >= 0.0);
  assert(r >= 0.0);

  if (z <= c_z) {
    /* Cylinder section */
    return {m_rad - r, 0};
  } else if ((r >= c_r) && (z >= m_half_length)) {
    /* Wall section */
    return {0, m_half_length - z};
  } else {
    /* Smoothing area */
    /* Vector to center of torus segment */
    auto const dr = c_r - r;
    auto const dz = c_z - z;

    /* Rescale to surface */
    auto const d = std::sqrt(dr * dr + dz * dz);
    auto const fac = (d - m_smoothing_rad) / d;

    return {fac * dr, fac * dz};
  }
}

int SimplePore::calculate_dist(const double *ppos, double *dist,
                               double *vec) const {
  /* Coordinate transform to cylinder coords
     with origin at m_center. */
  Vector3d const c_dist = Vector3d(ppos, ppos + 3) - m_center;
  auto const z = e_z * c_dist;
  auto const r_vec = c_dist - z * e_z;
  auto const r = r_vec.norm();

  /* If exactly on the axis, chose e_r orthogonal
     to e_z. */
  auto const e_r = (r == 0) ? e_r_axis : r_vec / r;

  /* The pore has mirror symmetry in z with regard to
     the center in the {r,z} system. We calculate always
     for the z > 0 case, and flip the result if appropriate. */
  double dr, dz;
  std::tie(dr, dz) = dist_half_pore(r, std::abs(z));

  if (z <= 0.0) {
    dz *= -1;
  }

  *dist = std::sqrt(dr * dr + dz * dz);
  for (int i = 0; i < 3; i++) {
    vec[i] = -dr * e_r[i] + -dz * e_z[i];
  }

  return 0;
}
}

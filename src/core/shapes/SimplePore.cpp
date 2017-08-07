#include "SimplePore.hpp"
#include "random.hpp"

#include <cassert>

namespace Shapes {
std::pair<double, double> SimplePore::dist_half_pore(double r, double z) const {
  assert(z > 0.0);
  assert(r > 0.0);

  if (z <= c_z) {
    /* Cylinder section */
    return {m_rad - r, 0};
  } else if (r >= c_r) {
    /* Wall section */
    return {0, z - m_half_length};
  } else {
    /* Vector to center */
    auto const dr = c_r - r;
    auto const dz = z - c_z;

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
  auto r_vec = c_dist - z * e_z;
  auto const r = r_vec.norm();

  /* Exactly on axis, chose
   a r_vec orthogonal to e_z. */
  if (r == 0) {
    if ((Vector3d{1., 0., 0} * e_z) < 1.)
      r_vec = Vector3d{1., 0., 0};
    else
      r_vec = Vector3d{0., 1., 0};
  }

  auto const e_r = r_vec / r;

  double dr, dz;
  std::tie(dr, dz) = dist_half_pore(r, std::abs(z));

  if (z <= 0.0) {
    dz *= -1;
  }

  *dist = std::sqrt(dr * dr + dz * dz);
  for (int i = 0; i < 3; i++) {
    vec[i] = dr * e_r[i] + dz * e_z[i];
  }

  return 0;
}
}

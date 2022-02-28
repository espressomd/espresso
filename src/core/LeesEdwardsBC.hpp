#pragma once

#include "utils/Vector.hpp"
#include <iostream>

#include <bitset>

struct LeesEdwardsBC {
  double pos_offset = 0.;
  double shear_velocity = 0.;
  int shear_direction = 0;
  int shear_plane_normal = 0;
  Utils::Vector3d distance(const Utils::Vector3d &d, const Utils::Vector3d &l,
                           const Utils::Vector3d &hal_l,
                           const Utils::Vector3d &l_inv,
                           const std::bitset<3> periodic) const {
    assert(shear_plane_normal != shear_direction);
    assert(shear_direction >= 0 and shear_direction <= 2);
    assert(shear_plane_normal >= 0 and shear_plane_normal <= 2);

    Utils::Vector3d n_jumps = {};
    Utils::Vector3d res = d;

    double n_le_crossings =
        std::round(res[shear_plane_normal] * l_inv[shear_plane_normal]);
    res[shear_direction] += n_le_crossings * pos_offset;

    for (int i : {0, 1, 2}) {
      if (periodic[i])
        n_jumps[i] = std::round(res[i] * l_inv[i]);
      res[i] -= n_jumps[i] * l[i];
    }

    return res;
  }
};

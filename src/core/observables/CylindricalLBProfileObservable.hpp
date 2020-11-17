/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP

#include "CylindricalProfileObservable.hpp"

#include <utils/Vector.hpp>
#include <utils/math/coordinate_transformation.hpp>
#include <utils/math/vec_rotate.hpp>
#include <utils/sampling.hpp>

namespace Observables {

class CylindricalLBProfileObservable : public CylindricalProfileObservable {
public:
  CylindricalLBProfileObservable(Utils::Vector3d const &center, Utils::Vector3d const &axis,
                                 int n_r_bins, int n_phi_bins, int n_z_bins,
                                 double min_r, double min_phi, double min_z,
                                 double max_r, double max_phi, double max_z,
                                 double sampling_density)
      : CylindricalProfileObservable(center, axis, min_r, max_r, min_phi,
                                     max_phi, min_z, max_z, n_r_bins,
                                     n_phi_bins, n_z_bins),
        sampling_density(sampling_density) {
    calculate_sampling_positions();
  }
  void calculate_sampling_positions() {
    sampling_positions = Utils::get_cylindrical_sampling_positions(
        limits[0], limits[1], limits[2], n_bins[0], n_bins[1], n_bins[2],
        sampling_density);
    for (auto &p : sampling_positions) {
      double theta;
      Utils::Vector3d rotation_axis;
      auto p_cart = Utils::transform_coordinate_cylinder_to_cartesian(
          p, Utils::Vector3d{{0.0, 0.0, 1.0}});
      // We have to rotate the coordinates since the utils function assumes
      // z-axis symmetry.
      std::tie(theta, rotation_axis) =
          Utils::rotation_params(Utils::Vector3d{{0.0, 0.0, 1.0}}, axis);
      p_cart = Utils::vec_rotate(rotation_axis, theta, p_cart);
      p = p_cart + center;
    }
  }
  std::vector<Utils::Vector3d> sampling_positions;
  double sampling_density;
};
} // Namespace Observables
#endif

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

#include <utility>
#include <utils/Vector.hpp>
#include <utils/math/coordinate_transformation.hpp>
#include <utils/math/vec_rotate.hpp>
#include <utils/sampling.hpp>

namespace Observables {

class CylindricalLBProfileObservable : public CylindricalProfileObservable {
public:
  CylindricalLBProfileObservable(
      std::shared_ptr<Utils::CylindricalTransformationParameters>
          transform_params,
      int n_r_bins, int n_phi_bins, int n_z_bins, double min_r, double max_r,
      double min_phi, double max_phi, double min_z, double max_z,
      double sampling_density)
      : CylindricalProfileObservable(std::move(transform_params), n_r_bins,
                                     n_phi_bins, n_z_bins, min_r, max_r,
                                     min_phi, max_phi, min_z, max_z),
        sampling_density(sampling_density) {
    calculate_sampling_positions();
  }
  void calculate_sampling_positions() {
    sampling_positions = Utils::get_cylindrical_sampling_positions(
        limits[0], limits[1], limits[2], n_bins[0], n_bins[1], n_bins[2],
        sampling_density);
    for (auto &p : sampling_positions) {
      auto p_cart = Utils::transform_coordinate_cylinder_to_cartesian(p);
      // We have to rotate the coordinates since the utils function assumes
      // z-axis symmetry.
      constexpr Utils::Vector3d z_axis{{0.0, 0.0, 1.0}};
      auto const theta = Utils::angle_between(z_axis, transform_params->axis());
      auto const rot_axis =
          Utils::vector_product(z_axis, transform_params->axis()).normalize();
      if (theta > std::numeric_limits<double>::epsilon())
        p_cart = Utils::vec_rotate(rot_axis, theta, p_cart);
      p = p_cart + transform_params->center();
    }
  }
  std::vector<Utils::Vector3d> sampling_positions;
  double sampling_density;
};
} // Namespace Observables
#endif

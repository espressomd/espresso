/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#pragma once

#include "CylindricalProfileObservable.hpp"

#include "SanityChecksLB.hpp"

#include <utils/Vector.hpp>
#include <utils/math/coordinate_transformation.hpp>
#include <utils/math/vec_rotate.hpp>
#include <utils/sampling.hpp>

#include <limits>
#include <memory>
#include <utility>
#include <vector>

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
        sampling_density(sampling_density), lb_sanity_checks() {}

  double sampling_density;

protected:
  mutable SanityChecksLB lb_sanity_checks;
  mutable std::vector<Utils::Vector3d> sampling_positions;
  mutable std::vector<Utils::Vector3d> sampling_positions_cyl;
  mutable std::vector<Utils::Vector3d> sampling_positions_cart;

public:
  void calculate_sampling_positions(auto const &box_geo, auto const &lb) const {
    lb_sanity_checks = SanityChecksLB(box_geo, lb);
    sampling_positions.clear();
    sampling_positions_cyl.clear();
    sampling_positions_cart.clear();

    auto const lb_position_checker = lb.make_lattice_position_checker(false);
    auto const lim = limits();
    auto const b = n_bins();
    auto const global_positions = Utils::get_cylindrical_sampling_positions(
        lim[0], lim[1], lim[2], b[0], b[1], b[2], sampling_density);
    for (auto const &p : global_positions) {
      auto p_cart = Utils::transform_coordinate_cylinder_to_cartesian(p);
      // We have to rotate the coordinates since the utils function assumes
      // z-axis symmetry.
      auto const z_axis = Utils::Vector3d{{0.0, 0.0, 1.0}};
      auto const theta = Utils::angle_between(z_axis, transform_params->axis());
      auto const rot_axis =
          Utils::vector_product(z_axis, transform_params->axis()).normalize();
      if (theta > std::numeric_limits<double>::epsilon())
        p_cart = Utils::vec_rotate(rot_axis, theta, p_cart);
      auto const pos = p_cart + transform_params->center();
      if (lb_position_checker(pos)) {
        sampling_positions.emplace_back(pos);
        sampling_positions_cart.emplace_back(p_cart);
        auto const pos_cyl = Utils::transform_coordinate_cartesian_to_cylinder(
            p_cart, transform_params->axis(), transform_params->orientation());
        sampling_positions_cyl.emplace_back(pos_cyl);
      }
    }
  }
};

} // Namespace Observables

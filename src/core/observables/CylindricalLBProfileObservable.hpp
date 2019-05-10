/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALLBPROFILEOBSERVABLE_HPP

#include "CylindricalProfileObservable.hpp"
#include <utils/Vector.hpp>
#include <utils/coordinate_transformation.hpp>
#include <utils/sampling.hpp>

using boost::math::constants::pi;
using Utils::Vector3d;

namespace Observables {

class CylindricalLBProfileObservable : public CylindricalProfileObservable {
public:
  void calculate_sampling_positions() {
    sampling_positions = Utils::get_cylindrical_sampling_positions(
        std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
        std::make_pair(min_phi, max_phi), n_r_bins, n_phi_bins, n_z_bins,
        sampling_density);
    for (auto &p : sampling_positions) {
      auto const p_cart =
          Utils::transform_coordinate_cylinder_to_cartesian(p, axis);
      p = p_cart - center;
    }
  }
  std::vector<Vector3d> sampling_positions;
  double sampling_density;
};
} // Namespace Observables
#endif

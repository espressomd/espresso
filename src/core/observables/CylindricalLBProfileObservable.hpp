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

#include <algorithm>
#include <boost/math/constants/constants.hpp>
#include <limits>
#include <math.h>

#include "CylindricalProfileObservable.hpp"
#include <utils/Vector.hpp>

using boost::math::constants::pi;
using Utils::Vector3d;

namespace Observables {

class CylindricalLBProfileObservable : public CylindricalProfileObservable {
public:
  void calculate_sampling_positions() {
    sampling_positions.clear();
    auto const r_range = max_r - min_r;
    auto const z_range = max_z - min_z;
    auto const delta_r = r_bin_size();
    auto const delta_phi = phi_bin_size();

    auto const smallest_bin_volume = pi<double>() * pow(min_r + delta_r, 2.0) *
                                     delta_phi / (2.0 * pi<double>());
    auto const min_n_samples = std::max(
        n_z_bins, static_cast<int>(smallest_bin_volume * sampling_density));
    auto const z_step = z_range / min_n_samples;

    // Take care of the smallest bin
    for (int i = 0; i < min_n_samples; ++i) {
      if (min_r < std::numeric_limits<double>::epsilon()) {
        // Place a sample in the center if r_min is 0.0
        sampling_positions.push_back(
            Vector3d{{0.0, 0.0, min_z + (i + 0.5) * z_step}});
      }
      for (int j = 0; j < n_phi_bins; ++j) {
        sampling_positions.push_back(
            Vector3d{{min_r + 0.5 * delta_r, min_phi + (j + 0.5) * delta_phi,
                      min_z + (i + 0.5) * z_step}});
      }
    }

    // Scale the number of samples for larger bins
    auto arc_length = [&delta_phi](int r_bin) {
      return delta_phi * 2.0 * pi<double>() * (r_bin + 1);
    };
    auto n_phi_samples = [&arc_length](int r_bin) {
      return arc_length(r_bin) / arc_length(0);
    };
    auto phis = [this, &n_phi_samples](int r_bin) {
      std::vector<double> phi_values;
      auto const phi_step =
          (max_phi - min_phi) / (n_phi_bins * n_phi_samples(r_bin));
      for (int i = 0; i < n_phi_bins * n_phi_samples(r_bin); ++i) {
        phi_values.push_back(min_phi + i * phi_step);
      }
      return phi_values;
    };
    // Calculate the sampling positions
    // Along z
    for (int i = 0; i < min_n_samples; ++i) {
      // Along r
      for (int j = 0; j < n_r_bins; ++j) {
        // Along phi
        for (auto const &phi : phis(j)) {
          sampling_positions.push_back(Vector3d{
              {min_r + (j + 1.5) * delta_r, phi, min_z + (i + 0.5) * z_step}});
        }
      }
    }
    // Transform to cartesian coordinates
    std::transform(
        sampling_positions.begin(), sampling_positions.end(),
        sampling_positions.begin(), [this](Vector3d const &p) -> Vector3d {
          return {{center[0] + p[0] * std::cos(p[1]),
                   center[1] + p[0] * std::sin(p[1]), center[2] + p[2]}};
        });
  }
  std::vector<Vector3d> sampling_positions;
  double sampling_density;
};
} // Namespace Observables
#endif

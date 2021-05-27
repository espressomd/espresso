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
#ifndef OBSERVABLES_LBPROFILEOBSERVABLE_HPP
#define OBSERVABLES_LBPROFILEOBSERVABLE_HPP

#include "ProfileObservable.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>

namespace Observables {

class LBProfileObservable : public ProfileObservable {
public:
  LBProfileObservable(double sampling_delta_x, double sampling_delta_y,
                      double sampling_delta_z, double sampling_offset_x,
                      double sampling_offset_y, double sampling_offset_z,
                      int n_x_bins, int n_y_bins, int n_z_bins, double min_x,
                      double max_x, double min_y, double max_y, double min_z,
                      double max_z, bool allow_empty_bins = false)
      : ProfileObservable(n_x_bins, n_y_bins, n_z_bins, min_x, max_x, min_y,
                          max_y, min_z, max_z),
        sampling_delta{sampling_delta_x, sampling_delta_y, sampling_delta_z},
        sampling_offset{sampling_offset_x, sampling_offset_y,
                        sampling_offset_z},
        allow_empty_bins(allow_empty_bins) {
    if (sampling_delta[0] <= 0.)
      throw std::domain_error("sampling_delta_x has to be > 0");
    if (sampling_delta[1] <= 0.)
      throw std::domain_error("sampling_delta_y has to be > 0");
    if (sampling_delta[2] <= 0.)
      throw std::domain_error("sampling_delta_z has to be > 0");
    if (sampling_offset[0] < 0.)
      throw std::domain_error("sampling_offset_x has to be >= 0");
    if (sampling_offset[1] < 0.)
      throw std::domain_error("sampling_offset_y has to be >= 0");
    if (sampling_offset[2] < 0.)
      throw std::domain_error("sampling_offset_z has to be >= 0");
    calculate_sampling_positions();
  }
  std::array<double, 3> sampling_delta;
  std::array<double, 3> sampling_offset;
  bool allow_empty_bins;
  std::vector<Utils::Vector3d> sampling_positions;
  void calculate_sampling_positions() {
    auto const lim = limits();
    sampling_positions.clear();
    assert(sampling_delta[0] > 0. and sampling_delta[1] > 0. and
           sampling_delta[2] > 0.);
    assert(sampling_offset[0] >= 0. and sampling_offset[1] >= 0. and
           sampling_offset[2] >= 0.);
    const auto n_samples_x = static_cast<size_t>(
        std::rint((lim[0].second - lim[0].first) / sampling_delta[0]));
    const auto n_samples_y = static_cast<size_t>(
        std::rint((lim[1].second - lim[1].first) / sampling_delta[1]));
    const auto n_samples_z = static_cast<size_t>(
        std::rint((lim[2].second - lim[2].first) / sampling_delta[2]));
    for (size_t x = 0; x < n_samples_x; ++x) {
      for (size_t y = 0; y < n_samples_y; ++y) {
        for (size_t z = 0; z < n_samples_z; ++z) {
          sampling_positions.push_back(Utils::Vector3d{
              {lim[0].first + sampling_offset[0] +
                   static_cast<double>(x) * sampling_delta[0],
               lim[1].first + sampling_offset[1] +
                   static_cast<double>(y) * sampling_delta[1],
               lim[2].first + sampling_offset[2] +
                   static_cast<double>(z) * sampling_delta[2]}});
        }
      }
    }
  }
};

} // namespace Observables

#endif

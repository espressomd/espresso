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
#ifndef OBSERVABLES_PROFILEOBSERVABLE_HPP
#define OBSERVABLES_PROFILEOBSERVABLE_HPP

#include "Observable.hpp"

#include <utils/math/make_lin_space.hpp>

#include <boost/range/algorithm.hpp>

#include <array>
#include <cstddef>
#include <utility>
#include <vector>

namespace Observables {

/** Cartesian profile observable */
class ProfileObservable : virtual public Observable {
protected:
  /** Range of the profile edges. */
  const std::array<std::pair<double, double>, 3> limits;
  /** Number of bins for each coordinate. */
  const std::array<size_t, 3> n_bins;

public:
  ProfileObservable(int n_x_bins, int n_y_bins, int n_z_bins, double min_x,
                    double max_x, double min_y, double max_y, double min_z,
                    double max_z)
      : limits{{std::make_pair(min_x, max_x), std::make_pair(min_y, max_y),
                std::make_pair(min_z, max_z)}},
        n_bins{{static_cast<size_t>(n_x_bins), static_cast<size_t>(n_y_bins),
                static_cast<size_t>(n_z_bins)}} {
    if (max_x <= min_x)
      throw std::runtime_error("max_x has to be > min_x");
    if (max_y <= min_y)
      throw std::runtime_error("max_y has to be > min_y");
    if (max_z <= min_z)
      throw std::runtime_error("max_z has to be > min_z");
  }

  std::vector<size_t> shape() const override {
    return {n_bins[0], n_bins[1], n_bins[2]};
  }

  auto get_n_bins() const { return n_bins; }

  auto get_limits() const { return limits; }
  /** Calculate the bin edges for each dimension */
  std::array<std::vector<double>, 3> edges() const {
    std::array<std::vector<double>, 3> profile_edges = {
        {std::vector<double>(n_bins[0] + 1), std::vector<double>(n_bins[1] + 1),
         std::vector<double>(n_bins[2] + 1)}};
    boost::copy(
        Utils::make_lin_space(limits[0].first, limits[0].second, n_bins[0] + 1),
        profile_edges[0].begin());
    boost::copy(
        Utils::make_lin_space(limits[1].first, limits[1].second, n_bins[1] + 1),
        profile_edges[1].begin());
    boost::copy(
        Utils::make_lin_space(limits[2].first, limits[2].second, n_bins[2] + 1),
        profile_edges[2].begin());
    return profile_edges;
  }
}; // namespace Observables

} // Namespace Observables
#endif

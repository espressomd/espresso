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
#ifndef OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP
#define OBSERVABLES_CYLINDRICALPROFILEOBSERVABLE_HPP

#include "Observable.hpp"

#include <utils/Vector.hpp>
#include <utils/math/make_lin_space.hpp>

#include <boost/range/algorithm.hpp>

#include <array>
#include <cstddef>
#include <vector>

namespace Observables {

/** Cylindrical profile observable */
class CylindricalProfileObservable : virtual public Observable {
public:
  CylindricalProfileObservable(Utils::Vector3d const &center,
                               Utils::Vector3d const &axis, double min_r,
                               double max_r, double min_phi, double max_phi,
                               double min_z, double max_z, int n_r_bins,
                               int n_phi_bins, int n_z_bins)
      : limits{{std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
         std::make_pair(min_z, max_z)}}, n_bins{{static_cast<size_t>(n_r_bins),
        static_cast<size_t>(n_phi_bins),
        static_cast<size_t>(n_z_bins)}}, center(center), axis(axis) {}
  /** Range of the profile edges. */
  std::array<std::pair<double, double>, 3> limits;
  /** Number of bins for each coordinate. */
  std::array<size_t, 3> n_bins;
  Utils::Vector3d center;
  Utils::Vector3d axis;

  std::vector<size_t> shape() const override {
    return {n_bins[0], n_bins[1], n_bins[2]};
  }

  /** Calculate the bin edges for each dimension */
  std::array<std::vector<double>, 3> edges() {
    std::array<std::vector<double>, 3> profile_edges = {
        {std::vector<double>(n_bins[0] + 1),
         std::vector<double>(n_bins[1] + 1),
         std::vector<double>(n_bins[2] + 1)}};
    boost::copy(Utils::make_lin_space(limits[0].first, limits[0].second, n_bins[0] + 1),
                profile_edges[0].begin());
    boost::copy(Utils::make_lin_space(limits[1].first, limits[1].second, n_bins[1] + 1),
                profile_edges[1].begin());
    boost::copy(Utils::make_lin_space(limits[2].first, limits[2].second, n_bins[2] + 1),
                profile_edges[2].begin());
    return profile_edges;
  }
};

} // Namespace Observables
#endif

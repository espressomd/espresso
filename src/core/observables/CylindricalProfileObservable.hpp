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

#include <boost/range/algorithm.hpp>
#include <vector>

#include "Observable.hpp"

#include <utils/Vector.hpp>
#include <utils/math/make_lin_space.hpp>

namespace Observables {

/** Cylindrical profile observable */
class CylindricalProfileObservable : virtual public Observable {
public:
  CylindricalProfileObservable(Utils::Vector3d const &center,
                               Utils::Vector3d const &axis, double min_r,
                               double max_r, double min_phi, double max_phi,
                               double min_z, double max_z, int n_r_bins,
                               int n_phi_bins, int n_z_bins)
      : center(center), axis(axis), min_r(min_r), max_r(max_r),
        min_phi(min_phi), max_phi(max_phi), min_z(min_z), max_z(max_z),
        n_r_bins(static_cast<size_t>(n_r_bins)),
        n_phi_bins(static_cast<size_t>(n_phi_bins)),
        n_z_bins(static_cast<size_t>(n_z_bins)){};
  Utils::Vector3d center;
  Utils::Vector3d axis;
  // Range of the profile edges.
  double min_r, max_r;
  double min_phi, max_phi;
  double min_z, max_z;
  // Number of bins for each coordinate.
  size_t n_r_bins, n_phi_bins, n_z_bins;

  std::vector<size_t> shape() const override {
    return {n_r_bins, n_phi_bins, n_z_bins};
  }

  /** Calculate the bin edges for each dimension */
  std::array<std::vector<double>, 3> edges() {
    std::array<std::vector<double>, 3> profile_edges = {
        {std::vector<double>(n_r_bins + 1), std::vector<double>(n_phi_bins + 1),
         std::vector<double>(n_z_bins + 1)}};
    boost::copy(Utils::make_lin_space(min_r, max_r, n_r_bins + 1),
                profile_edges[0].begin());
    boost::copy(Utils::make_lin_space(min_phi, max_phi, n_phi_bins + 1),
                profile_edges[1].begin());
    boost::copy(Utils::make_lin_space(min_z, max_z, n_z_bins + 1),
                profile_edges[2].begin());
    return profile_edges;
  }
};

} // Namespace Observables
#endif

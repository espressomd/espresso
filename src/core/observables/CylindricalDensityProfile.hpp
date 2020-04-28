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
#ifndef OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP

#include "CylindricalPidProfileObservable.hpp"
#include <utils/Histogram.hpp>
#include <utils/math/coordinate_transformation.hpp>

namespace Observables {
class CylindricalDensityProfile : public CylindricalPidProfileObservable {
public:
  using CylindricalPidProfileObservable::CylindricalPidProfileObservable;
  std::vector<double>
  evaluate(Utils::Span<const Particle *const> particles) const override {
    std::array<size_t, 3> n_bins{{n_r_bins, n_phi_bins, n_z_bins}};
    std::array<std::pair<double, double>, 3> limits{
        {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
         std::make_pair(min_z, max_z)}};
    Utils::CylindricalHistogram<double, 3> histogram(n_bins, 1, limits);

    for (auto p : particles) {
      histogram.update(Utils::transform_coordinate_cartesian_to_cylinder(
          folded_position(p->r.p, box_geo) - center, axis));
    }

    histogram.normalize();
    return histogram.get_histogram();
  }
  std::vector<size_t> shape() const override {
    return {n_r_bins, n_phi_bins, n_z_bins};
  }
};

} // Namespace Observables

#endif

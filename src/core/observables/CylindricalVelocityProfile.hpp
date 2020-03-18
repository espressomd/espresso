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
#ifndef OBSERVABLES_CYLINDRICALVELOCITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALVELOCITYPROFILE_HPP

#include "CylindricalPidProfileObservable.hpp"
#include <utils/Histogram.hpp>

namespace Observables {
class CylindricalVelocityProfile : public CylindricalPidProfileObservable {
public:
  using CylindricalPidProfileObservable::CylindricalPidProfileObservable;

  std::vector<double>
  evaluate(Utils::Span<const Particle *const> particles) const override {
    std::array<size_t, 3> n_bins{{n_r_bins, n_phi_bins, n_z_bins}};
    std::array<std::pair<double, double>, 3> limits{
        {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
         std::make_pair(min_z, max_z)}};
    Utils::CylindricalHistogram<double, 3> histogram(n_bins, 3, limits);

    for (auto p : particles) {
      auto const pos = folded_position(p->r.p, box_geo) - center;
      histogram.update(
          Utils::transform_coordinate_cartesian_to_cylinder(pos, axis),
          Utils::transform_vector_cartesian_to_cylinder(p->m.v, axis, pos));
    }

    auto hist_tmp = histogram.get_histogram();
    auto tot_count = histogram.get_tot_count();
    for (size_t ind = 0; ind < hist_tmp.size(); ++ind) {
      if (tot_count[ind] > 0) {
        hist_tmp[ind] /= tot_count[ind];
      }
    }
    return hist_tmp;
  }
  std::vector<size_t> shape() const override {
    return {n_r_bins, n_phi_bins, n_z_bins, 3};
  }
};

} // Namespace Observables

#endif

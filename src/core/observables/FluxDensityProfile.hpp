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
#ifndef OBSERVABLES_FLUXDENSITYPROFILE_HPP
#define OBSERVABLES_FLUXDENSITYPROFILE_HPP

#include "PidProfileObservable.hpp"

#include <vector>

namespace Observables {
class FluxDensityProfile : public PidProfileObservable {
public:
  using PidProfileObservable::PidProfileObservable;
  std::vector<size_t> shape() const override {
    return {n_x_bins, n_y_bins, n_z_bins, 3};
  }

  std::vector<double>
  evaluate(Utils::Span<const Particle *const> particles) const override {
    std::array<size_t, 3> n_bins{{n_x_bins, n_y_bins, n_z_bins}};
    std::array<std::pair<double, double>, 3> limits{
        {std::make_pair(min_x, max_x), std::make_pair(min_y, max_y),
         std::make_pair(min_z, max_z)}};
    Utils::Histogram<double, 3> histogram(n_bins, 3, limits);

    for (auto p : particles) {
      auto const ppos = folded_position(p->r.p, box_geo);
      histogram.update(ppos, p->m.v);
    }
    histogram.normalize();
    return histogram.get_histogram();
  }
};

} // Namespace Observables

#endif

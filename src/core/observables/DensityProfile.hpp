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
#ifndef OBSERVABLES_DENSITYPROFILE_HPP
#define OBSERVABLES_DENSITYPROFILE_HPP

#include "PidProfileObservable.hpp"
#include <utils/Histogram.hpp>
#include <vector>

namespace Observables {

class DensityProfile : public PidProfileObservable {
public:
  std::vector<double> operator()(PartCfg &partCfg) const override {
    std::array<size_t, 3> n_bins{{static_cast<size_t>(n_x_bins),
                                  static_cast<size_t>(n_y_bins),
                                  static_cast<size_t>(n_z_bins)}};
    std::array<std::pair<double, double>, 3> limits{
        {std::make_pair(min_x, max_x), std::make_pair(min_y, max_y),
         std::make_pair(min_z, max_z)}};
    Utils::Histogram<double, 3> histogram(n_bins, 1, limits);
    for (auto const &id : ids()) {
      histogram.update(folded_position(partCfg[id]));
    }
    histogram.normalize();
    return histogram.get_histogram();
  }
};
} // Namespace Observables

#endif

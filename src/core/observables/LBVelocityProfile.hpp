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
#ifndef OBSERVABLES_LBVELOCITYPROFILE_HPP
#define OBSERVABLES_LBVELOCITYPROFILE_HPP

#include "LBProfileObservable.hpp"
#include "particle_data.hpp"

#include <vector>

namespace Observables {

class LBVelocityProfile : public LBProfileObservable {
public:
  int n_values() const override { return 3 * n_x_bins * n_y_bins * n_z_bins; }
  std::vector<double> operator()(PartCfg &partCfg) const override;
};

} // Namespace Observables

#endif

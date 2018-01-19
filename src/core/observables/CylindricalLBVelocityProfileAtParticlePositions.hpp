/*
  Copyright (C) 2016,2017 The ESPResSo project

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
#ifndef OBSERVABLES_CYLINDRICALLBVELOCITYPROFILEATPARTICLEPOSITIONS_HPP
#define OBSERVABLES_CYLINDRICALLBVELOCITYPROFILEATPARTICLEPOSITIONS_HPP

#include "CylindricalProfileObservable.hpp"
#include "partCfg_global.hpp"
#include "utils/Histogram.hpp"

namespace Observables {
class CylindricalLBVelocityProfileAtParticlePositions
    : public CylindricalProfileObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override;
  virtual int n_values() const override {
    return 3 * n_r_bins * n_phi_bins * n_z_bins;
  }
};

} // Namespace Observables

#endif

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
#ifndef OBSERVABLES_PARTICLEANGULARVELOCITIES_HPP
#define OBSERVABLES_PARTICLEANGULARVELOCITIES_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"

#include "rotation.hpp"
#include <vector>

namespace Observables {

class ParticleAngularVelocities : public PidObservable {
public:
  using PidObservable::PidObservable;
  std::vector<double> evaluate(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
#ifdef ROTATION
      const Utils::Vector3d omega =
          convert_vector_body_to_space(partCfg[i], partCfg[i].m.omega);
      res[3 * i + 0] = omega[0];
      res[3 * i + 1] = omega[1];
      res[3 * i + 2] = omega[2];
#endif
    }
    return res;
  }
  int n_values() const override { return 3 * ids().size(); }
};

} // Namespace Observables
#endif

/*
Copyright (C) 2019 The ESPResSo project

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
#ifndef OBSERVABLES_PARTICLEDISTANCES_HPP
#define OBSERVABLES_PARTICLEDISTANCES_HPP

#include "PidObservable.hpp"
#include <utils/Vector.hpp>

#include <vector>

namespace Observables {

/** Calculate distances between particles in a polymer.
 *  For @f$ n @f$ bonded particles, return the @f$ n-1 @f$ distances separating
 *  them.
 */
class ParticleDistances : public PidObservable {
public:
  std::vector<double> evaluate(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0, end = n_values(); i < end; i++) {
      auto v = get_mi_vector(partCfg[ids()[i]].r.p, partCfg[ids()[i + 1]].r.p);
      res[i] = v.norm();
    }
    return res;
  }
  int n_values() const override { return ids().size() - 1; }
};

} // Namespace Observables

#endif

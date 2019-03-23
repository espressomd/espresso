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
#ifndef OBSERVABLES_PARTICLEDIHEDRALS_HPP
#define OBSERVABLES_PARTICLEDIHEDRALS_HPP

#include "PidObservable.hpp"
#include "utils/Vector.hpp"

#include <cmath>
#include <vector>

namespace Observables {

/** Calculate dihedral angles between particles in a polymer.
 *  For \f$n\f$ bonded particles, return the \f$n-3\f$ signed dihedral angles
 *  along the chain, in radians.
 */
class ParticleDihedrals : public PidObservable {
public:
  std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    auto r1 = (partCfg[ids()[1]].r.p - partCfg[ids()[0]].r.p);
    auto r2 = (partCfg[ids()[2]].r.p - partCfg[ids()[1]].r.p);
    auto c1 = vector_product(r1, r2);
    for (int i = 0, end = n_values(); i < end; i++) {
      auto r3 = (partCfg[ids()[i + 3]].r.p - partCfg[ids()[i + 2]].r.p);
      auto c2 = vector_product(r2, r3);
      /* the 2-argument arctangent returns an angle in the range [-pi, pi] that
       * allows for an unambiguous determination of the 4th particle position */
      res[i] = atan2(vector_product(c1, c2) * r2 / r2.norm(), c1 * c2);
      r1 = std::move(r2);
      r2 = std::move(r3);
      c1 = std::move(c2);
    }
    return res;
  }
  int n_values() const override { return ids().size() - 3; }
};

} // Namespace Observables

#endif

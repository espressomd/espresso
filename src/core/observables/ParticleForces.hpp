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
#ifndef OBSERVABLES_PARTICLEFORCES_HPP
#define OBSERVABLES_PARTICLEFORCES_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include <vector>

namespace Observables {

/** Extract particle forces.
 *  For \f$n\f$ particles, return \f$3 n\f$ forces ordered as
 *  \f$(f_x^1, f_y^1, f_z^1, \dots, f_x^n, f_y^n, f_z^n)\f$.
 */
class ParticleForces : public PidObservable {
public:
  std::vector<double> evaluate(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
      res[3 * i + 0] = partCfg[ids()[i]].f.f[0];
      res[3 * i + 1] = partCfg[ids()[i]].f.f[1];
      res[3 * i + 2] = partCfg[ids()[i]].f.f[2];
    }
    return res;
  };
  int n_values() const override { return 3 * ids().size(); }
};

} // Namespace Observables
#endif

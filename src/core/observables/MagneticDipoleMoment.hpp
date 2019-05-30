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
#ifndef OBSERVABLES_MAGNETICDIPOLEMOMENT_HPP
#define OBSERVABLES_MAGNETICDIPOLEMOMENT_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class MagneticDipoleMoment : public PidObservable {
public:
  int n_values() const override { return 3; };
  std::vector<double> evaluate(PartCfg &partCfg) const override {
    std::vector<double> res(n_values(), 0.0);
    for (int i : ids()) {
#ifdef DIPOLES
      res[0] += partCfg[i].calc_dip()[0];
      res[1] += partCfg[i].calc_dip()[1];
      res[2] += partCfg[i].calc_dip()[2];
#endif
    }
    return res;
  }
};

} // Namespace Observables

#endif

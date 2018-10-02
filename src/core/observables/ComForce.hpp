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
#ifndef OBSERVABLES_ComForce_HPP
#define OBSERVABLES_ComForce_HPP

#include "PidObservable.hpp"

#include "integrate.hpp"

namespace Observables {

class ComForce : public PidObservable {
public:
  int n_values() const override { return 3; }
  std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
      res[0] += partCfg[ids()[i]].f->f[0] * partCfg[ids()[i]].p->mass;
      res[1] += partCfg[ids()[i]].f->f[1] * partCfg[ids()[i]].p->mass;
      res[2] += partCfg[ids()[i]].f->f[2] * partCfg[ids()[i]].p->mass;
    }
    return res;
  };
};

} // Namespace Observables
#endif

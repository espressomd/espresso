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
#ifndef OBSERVABLES_COMPOSITION_HPP
#define OBSERVABLES_COMPOSITION_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {
class ComPosition : public PidObservable {
public:
  using PidObservable::PidObservable;
  std::vector<size_t> shape() const override { return {3}; }

  std::vector<double>
  evaluate(Utils::Span<const Particle *const> particles) const override {
    std::vector<double> res(n_values());
    double total_mass = 0;
    for (auto p : particles) {
      if (p->p.is_virtual)
        continue;
      double mass = p->p.mass;
      res[0] += mass * p->r.p[0];
      res[1] += mass * p->r.p[1];
      res[2] += mass * p->r.p[2];
      total_mass += mass;
    }
    res[0] /= total_mass;
    res[1] /= total_mass;
    res[2] /= total_mass;
    return res;
  };
};

} // Namespace Observables
#endif

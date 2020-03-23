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
#ifndef OBSERVABLES_PARTICLEFORCES_HPP
#define OBSERVABLES_PARTICLEFORCES_HPP

#include "Particle.hpp"
#include "PidObservable.hpp"
#include "integrate.hpp"
#include <vector>

namespace Observables {

/** Extract particle forces.
 *  For \f$n\f$ particles, return \f$3 n\f$ forces ordered as
 *  \f$(f_x^1, f_y^1, f_z^1, \dots, f_x^n, f_y^n, f_z^n)\f$.
 */
class ParticleForces : public PidObservable {
public:
  using PidObservable::PidObservable;

  std::vector<double>
  evaluate(Utils::Span<const Particle *const> particles) const override {
    std::vector<double> res(n_values());
    for (size_t i = 0; i < particles.size(); i++) {
      res[3 * i + 0] = particles[i]->f.f[0];
      res[3 * i + 1] = particles[i]->f.f[1];
      res[3 * i + 2] = particles[i]->f.f[2];
    }
    return res;
  };
  std::vector<size_t> shape() const override { return {ids().size(), 3}; }
};

} // Namespace Observables
#endif

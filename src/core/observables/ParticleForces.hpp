/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include <cstddef>
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
  evaluate(ParticleReferenceRange particles,
           const ParticleObservables::traits<Particle> &) const override {
    std::vector<double> res(n_values());
    std::size_t i = 0;
    for (auto const &p : particles) {
      auto const &f = p.get().f.f;
      res[3 * i + 0] = f[0];
      res[3 * i + 1] = f[1];
      res[3 * i + 2] = f[2];
      i++;
    }
    return res;
  };
  std::vector<std::size_t> shape() const override { return {ids().size(), 3}; }
};

} // Namespace Observables
#endif

/*
 * Copyright (C) 2019 The ESPResSo project
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
#ifndef OBSERVABLES_PARTICLEDISTANCES_HPP
#define OBSERVABLES_PARTICLEDISTANCES_HPP

#include "BoxGeometry.hpp"
#include "PidObservable.hpp"
#include "grid.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <cstddef>
#include <vector>

namespace Observables {

/** Calculate distances between particles in a polymer.
 *  For @f$ n @f$ bonded particles, return the @f$ n-1 @f$ distances separating
 *  them.
 */
class ParticleDistances : public PidObservable {
public:
  using PidObservable::PidObservable;

  std::vector<double>
  evaluate(Utils::Span<std::reference_wrapper<const Particle>> particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    std::vector<double> res(n_values());

    for (size_t i = 0, end = n_values(); i < end; i++) {
      auto const v = box_geo.get_mi_vector(traits.position(particles[i]),
                                           traits.position(particles[i + 1]));
      res[i] = v.norm();
    }
    return res;
  }
  std::vector<size_t> shape() const override { return {ids().size() - 1}; }
};

} // Namespace Observables

#endif

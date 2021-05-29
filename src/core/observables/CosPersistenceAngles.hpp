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
#ifndef OBSERVABLES_PERSISTENCEANGLES_HPP
#define OBSERVABLES_PERSISTENCEANGLES_HPP

#include "BoxGeometry.hpp"
#include "PidObservable.hpp"
#include "grid.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Observables {

/** Calculate bond angles in a polymer.
 *  The \a ith entry in the result vector corresponds to the
 *  averaged cosine of the angle between bonds that are \a i bonds apart.
 */
class CosPersistenceAngles : public PidObservable {
public:
  using PidObservable::PidObservable;
  explicit CosPersistenceAngles(std::vector<int> ids)
      : PidObservable(std::move(ids)) {
    if (this->ids().size() < 3)
      throw std::runtime_error("At least 3 particles are required");
  }

  std::vector<double>
  evaluate(Utils::Span<std::reference_wrapper<const Particle>> particles,
           const ParticleObservables::traits<Particle> &traits) const override {
    auto const no_of_angles = n_values();
    std::vector<double> angles(no_of_angles);
    auto const no_of_bonds = n_values() + 1;
    std::vector<Utils::Vector3d> bond_vectors(no_of_bonds);
    auto get_bond_vector = [&](auto index) {
      return box_geo.get_mi_vector(traits.position(particles[index + 1]),
                                   traits.position(particles[index]));
    };
    for (size_t i = 0; i < no_of_bonds; ++i) {
      auto const tmp = get_bond_vector(i);
      bond_vectors[i] = tmp / tmp.norm();
    }
    // calculate angles between neighbouring bonds, next neighbours, etc...
    for (size_t i = 0; i < no_of_angles; ++i) {
      auto average = 0.0;
      for (size_t j = 0; j < no_of_angles - i; ++j) {
        average += bond_vectors[j] * bond_vectors[j + i + 1];
      }
      angles[i] = average / static_cast<double>(no_of_angles - i);
    }

    return angles;
  }
  std::vector<size_t> shape() const override {
    assert(ids().size() >= 2);
    return {ids().size() - 2};
  }
};

} // Namespace Observables

#endif

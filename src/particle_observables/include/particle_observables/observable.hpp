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
#ifndef SRC_PARTICLE_OBSERVABLES_OBSERVABLE_HPP
#define SRC_PARTICLE_OBSERVABLES_OBSERVABLE_HPP

#include "algorithms.hpp"
#include "properties.hpp"

namespace ParticleObservables {
/**
 * @brief Meta-Observable that returns the product of two
 *        other observables.
 *
 * The operand observables are stored by privately deriving
 * from them to get the empty case optimization if they do
 * not have state.
 *
 * @tparam Left left operand of the product.
 * @tparam Right right operand of the product.
 */
template <class Left, class Right> struct Product : Left, Right {
  Product(Left left = {}, Right right = {})
      : Left(std::move(left)), Right(std::move(right)) {}

  template <class Particle, class Traits = default_traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits const &traits = {}) const {
    return Left::template operator()<Particle, Traits>(p, traits) *
           Right::template operator()<Particle, Traits>(p, traits);
  }
};

using Momentum = Product<Mass, Velocity>;
using AverageMomentum = Average<Momentum>;
using CenterOfMassPosition = WeightedAverage<Position, Mass>;
using CenterOfMassVelocity = WeightedAverage<Velocity, Mass>;
using Forces = Map<Force>;
using Positions = Map<Position>;
using Velocities = Map<Velocity>;
using Directors = Map<Director>;
using DipoleFields = Map<DipoleField>;
using BodyVelocities = Map<BodyVelocity>;
using AngularVelocities = Map<AngularVelocity>;
using BodyAngularVelocities = Map<BodyAngularVelocity>;
} // namespace ParticleObservables

#endif // SRC_PARTICLE_OBSERVABLES_OBSERVABLE_HPP

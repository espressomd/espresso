#ifndef INCLUDE_OBSERVABLES_OBSERVABLE_HPP
#define INCLUDE_OBSERVABLES_OBSERVABLE_HPP

#include "algorithms.hpp"
#include "properties.hpp"

namespace Observables {
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

  template <class Particle, class Traits = traits<Particle>>
  decltype(auto) operator()(Particle const &p,
                            Traits const &traits = {}) const {
    return Left::template operator()<Particle, Traits>(p, traits) *
           Right::template operator()<Particle, Traits>(p, traits);
  }
};

using Momentum = Product<Mass, Velocity>;
template <class ParticleRange>
using AverageMomentum = Average<ParticleRange, Momentum>;
template <class ParticleRange>
using CenterOfMass = WeightedAverage<ParticleRange, Position, Mass>;
template <class ParticleRange>
using CenterOfMassVelocity = WeightedAverage<ParticleRange, Velocity, Mass>;
} // namespace Observables

#endif // INCLUDE_OBSERVABLES_OBSERVABLE_HPP

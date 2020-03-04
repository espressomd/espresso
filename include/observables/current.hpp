#ifndef INCLUDE_OBSERVABLES_CURRENT_HPP
#define INCLUDE_OBSERVABLES_CURRENT_HPP

#include <algorithms/algorithms.hpp>
#include <traits/particle.hpp>

namespace Observables {

struct Current {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) {
    return Algorithms::weighted_sum(particles, Traits::Particle::Velocity(),
                                    Traits::Particle::Charge());
  }
};

} // namespace Observables

#endif // INCLUDE_OBSERVABLES_CURRENT_HPP

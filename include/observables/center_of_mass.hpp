#ifndef OBSERVABLES_CENTER_OF_MASS_HPP
#define OBSERVABLES_CENTER_OF_MASS_HPP

#include <algorithms/algorithms.hpp>
#include <traits/particle.hpp>

namespace Observables {

struct CenterOfMass {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) {
    return Algorithms::weighted_average(particles, Traits::Particle::Position(),
                                        Traits::Particle::Mass());
  }
};

} // namespace Observables
#endif

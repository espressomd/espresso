#ifndef INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP
#define INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP

#include <algorithms/algorithms.hpp>
#include <traits/particle.hpp>

namespace Observables {

struct DipoleMoment {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) {
    return Algorithms::weighted_sum(particles, Traits::Particle::Position(),
                                    Traits::Particle::Charge());
  }
};

} // namespace Observables

#endif // INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP

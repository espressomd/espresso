#ifndef INCLUDE_OBSERVABLES_CURRENT_HPP
#define INCLUDE_OBSERVABLES_CURRENT_HPP

#include <algorithms/algorithms.hpp>
#include <observables/observable.hpp>
#include <traits/particle.hpp>

namespace Observables {

using Current =
    Observables::Observable<Algorithms::WeightedSum, Traits::Particle::Velocity,
                            Traits::Particle::Charge>;

} // namespace Observables

#endif // INCLUDE_OBSERVABLES_CURRENT_HPP

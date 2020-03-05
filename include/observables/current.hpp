#ifndef INCLUDE_OBSERVABLES_CURRENT_HPP
#define INCLUDE_OBSERVABLES_CURRENT_HPP

#include <observables/algorithms.hpp>
#include <observables/observable.hpp>

namespace Observables {

using Current =
    Observables::Observable<Algorithms::WeightedSum, Traits::Particle::Velocity,
                            Traits::Particle::Charge>;

} // namespace Observables

#endif // INCLUDE_OBSERVABLES_CURRENT_HPP

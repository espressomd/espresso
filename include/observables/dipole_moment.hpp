#ifndef INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP
#define INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP

#include <observables/algorithms.hpp>
#include <observables/observable.hpp>

namespace Observables {

using DipoleMoment =
    Observables::Observable<Algorithms::WeightedSum, Traits::Particle::Position,
                            Traits::Particle::Charge>;

} // namespace Observables

#endif // INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP

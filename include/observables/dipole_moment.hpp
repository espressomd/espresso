#ifndef INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP
#define INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP

#include <algorithms/algorithms.hpp>
#include <observables/observable.hpp>
#include <traits/particle.hpp>

namespace Observables {

using DipoleMoment =
    Observables::Observable<Algorithms::WeightedSum, Traits::Particle::Position,
                            Traits::Particle::Charge>;

} // namespace Observables

#endif // INCLUDE_OBSERVABLES_DIPOLE_MOMENT_HPP

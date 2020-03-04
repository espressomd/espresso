#ifndef OBSERVABLES_CENTER_OF_MASS_HPP
#define OBSERVABLES_CENTER_OF_MASS_HPP

#include <algorithms/algorithms.hpp>
#include <observables/observable.hpp>
#include <traits/particle.hpp>

namespace Observables {

using CenterOfMass =
    Observables::Observable<Algorithms::WeightedAverage,
                            Traits::Particle::Position, Traits::Particle::Mass>;
} // namespace Observables
#endif

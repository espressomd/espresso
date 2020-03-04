#ifndef INCLUDE_OBSERVABLES_OBSERVABLE_HPP
#define INCLUDE_OBSERVABLES_OBSERVABLE_HPP

namespace Observables {

template <class Alg, class ValueOp, class WeightOp> struct Observable {
  template <class ParticleRange>
  auto operator()(ParticleRange const &particles) {
    return Alg()(particles, ValueOp(), WeightOp());
  }
};

} // namespace Observables

#endif // INCLUDE_OBSERVABLES_OBSERVABLE_HPP

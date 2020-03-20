#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP
#include <numeric>
#include <utility>

namespace Observables {
namespace detail {
struct One {
  template <class Particle> auto operator()(Particle const &p) { return 1; }
};

template <class ParticleRange, class ValueOp, class WeightOp>
struct WeightedSum {
  auto operator()(ParticleRange const &particles) {
    using particle_type = typename ParticleRange::value_type;
    using value_op_type = decltype(ValueOp{}(std::declval<particle_type>()));
    using weight_op_type = decltype(WeightOp{}(std::declval<particle_type>()));
    auto func = [](auto sum, auto const &p) {
      auto const w = WeightOp{}(p);
      return std::make_pair(sum.first + ValueOp{}(p)*w, sum.second + w);
    };

    return std::accumulate(std::begin(particles), std::end(particles),
                           std::pair<value_op_type, weight_op_type>(), func);
  }
};
} // namespace detail

template <class ParticleRange, class ValueOp, class WeightOp>
struct WeightedSum {
  auto operator()(ParticleRange const &particles) {
    return detail::WeightedSum<ParticleRange, ValueOp, WeightOp>()(particles)
        .first;
  }
};

template <class ParticleRange, class ValueOp> struct Sum {
  auto operator()(ParticleRange const &particles) {
    return detail::WeightedSum<ParticleRange, ValueOp, detail::One>()(particles)
        .first;
  }
};

template <class ParticleRange, class ValueOp, class WeightOp>
struct WeightedAverage {
  auto operator()(ParticleRange const &particles) {

    auto const ws =
        detail::WeightedSum<ParticleRange, ValueOp, WeightOp>()(particles);
    return ws.first / ws.second;
  }
};

template <class ParticleRange, class ValueOp> struct Average {
  auto operator()(ParticleRange const &particles) {
    return WeightedAverage<ParticleRange, ValueOp, detail::One>()(particles);
  }
};

template <class ParticleRange, class ValueOp> struct Collect {
  template <class OutputIterator>
  void operator()(ParticleRange const &particles, OutputIterator out) {
    std::transform(std::begin(particles), std::end(particles), out,
                   [](auto const &p) { return ValueOp{}(p); });
  }
};
} // namespace Observables
#endif // ALGORITHMS_HPP

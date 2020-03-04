#ifndef ALGORITHMS_HPP
#define ALGORITHMS_HPP
#include <numeric>
#include <utility>

namespace Algorithms {

namespace detail {

struct WeightedSum {
  template <class ParticleRange, class ValueOp, class WeightOp>
  auto operator()(ParticleRange const &particles, ValueOp value_op,
                  WeightOp weight_op) {
    using particle_type = typename ParticleRange::value_type;
    using value_op_type = decltype(value_op(std::declval<particle_type>()));
    using weight_op_type = decltype(weight_op(std::declval<particle_type>()));
    auto func = [&value_op, &weight_op](auto sum, auto const &p) {
      auto const w = weight_op(p);
      return std::make_pair(sum.first + value_op(p) * w, sum.second + w);
    };

    return std::accumulate(std::begin(particles), std::end(particles),
                           std::pair<value_op_type, weight_op_type>(), func);
  }
};

} // namespace detail

struct WeightedSum {
  template <class ParticleRange, class ValueOp, class WeightOp>
  auto operator()(ParticleRange const &particles, ValueOp &&value_op,
                  WeightOp &&weight_op) {
    return detail::WeightedSum()(particles, std::forward<ValueOp>(value_op),
                                 std::forward<WeightOp>(weight_op))
        .first;
  }
};

struct WeightedAverage {
  template <class ParticleRange, class ValueOp, class WeightOp>
  auto operator()(ParticleRange const &particles, ValueOp &&value_op,
                  WeightOp &&weight_op) {

    auto const ws =
        detail::WeightedSum()(particles, std::forward<ValueOp>(value_op),
                              std::forward<WeightOp>(weight_op));
    return ws.first / ws.second;
  }
};

struct Average {
  template <class ParticleRange, class ValueOp>
  auto operator()(ParticleRange const &particles, ValueOp value_op) {
    auto one = [](auto const &p) { return 1; };
    return WeightedAverage()(particles, value_op, one);
  }
};

} // namespace Algorithms
#endif // ALGORITHMS_HPP

#ifndef UTILS_CHECKS_CHARGE_NEUTRALITY_HPP
#define UTILS_CHECKS_CHARGE_NEUTRALITY_HPP

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/sum_kahan.hpp>

#include <cmath>
#include <limits>
#include <stdexcept>

namespace Utils {
template <typename ParticleRange>
bool check_charge_neutrality(
    ParticleRange &prange,
    double relative_tolerance = 1e-12) {
  using namespace boost::accumulators;
  using KahanSum = accumulator_set<double, features<tag::sum_kahan>>;

  KahanSum q_sum;
  auto q_min = std::numeric_limits<double>::infinity();

  for (auto const &p : prange) {
    auto const &q = p.p.q;

    if (q) {
      q_sum(q);
      q_min = std::min(q_min, std::abs(q));
    }
  }

  auto const excess_ratio = std::abs(sum_kahan(q_sum)) / q_min;

  return excess_ratio <= relative_tolerance;
}
}
#endif

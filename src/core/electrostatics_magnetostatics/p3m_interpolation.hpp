#ifndef ESPRESSO_P3M_INTERPOLATION_HPP
#define ESPRESSO_P3M_INTERPOLATION_HPP

#include <utils/Span.hpp>
#include <utils/index.hpp>

#include <boost/range/algorithm/copy.hpp>

#include <cassert>
#include <tuple>
#include <vector>

struct p3m_interpolation_weights {
  size_t m_cao = 0;
  /** Charge fractions for mesh assignment. */
  std::vector<double> ca_frac;
  /** index of first mesh point for charge assignment. */
  std::vector<int> ca_fmp;

  /**
   * @brief Number of points in the cache.
   * @return Number of points currently in the cache.
   */
  auto size() const { return ca_fmp.size(); }

  /**
   * @brief Charge assignment order the weights are for.
   * @return The charge assigment order.
   */
  auto cao() const { return m_cao; }

  /**
   * @brief Push back weights for one point.
   *
   * @param q_ind Mesh index.
   * @param w_x Weights in first direction.
   * @param w_y Weights in second direction.
   * @param w_z Weights in third direction.
   */
  void store(int q_ind, Utils::Span<const double> w_x,
             Utils::Span<const double> w_y, Utils::Span<const double> w_z) {
    ca_fmp.push_back(q_ind);
    auto it = std::back_inserter(ca_frac);
    boost::copy(w_x, it);
    boost::copy(w_y, it);
    boost::copy(w_z, it);
  }

  auto load(size_t i) const {
    using Utils::make_const_span;

    assert(i < size());

    auto const ind = ca_frac.data() + 3 * i * m_cao;
    return std::make_tuple(ca_fmp[i], make_const_span(ind + 0 * m_cao, m_cao),
                           make_const_span(ind + 1 * m_cao, m_cao),
                           make_const_span(ind + 2 * m_cao, m_cao));
  }

  /**
   * @brief Reset the cache.
   *
   * @param cao Interpolation order.
   */
  void reset(int cao) {
    m_cao = cao;
    ca_frac.clear();
    ca_fmp.clear();
  }
};

#endif // ESPRESSO_P3M_INTERPOLATION_HPP

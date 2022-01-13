/*
 * Copyright (C) 2016-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */
#ifndef UTILS_HISTOGRAM_HPP
#define UTILS_HISTOGRAM_HPP

#include "utils/Span.hpp"

#include <boost/multi_array.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

namespace Utils {

/**
 * \brief Histogram in Cartesian coordinates.
 * \tparam T     Histogram data type.
 * \tparam N     Histogram data dimensionality.
 * \tparam M     Coordinates data dimensionality.
 * \tparam U     Coordinates data type.
 */
template <typename T, std::size_t N, std::size_t M = 3, typename U = double>
class Histogram {
  using array_type = boost::multi_array<T, M + 1>;
  using count_type = boost::multi_array<std::size_t, M + 1>;
  using array_index = typename array_type::index;

public:
  /**
   * \brief Histogram constructor.
   * \param n_bins  the number of bins in each histogram dimension.
   * \param limits  the minimum/maximum data values to consider for the
   *        histogram.
   */
  explicit Histogram(std::array<std::size_t, M> n_bins,
                     std::array<std::pair<U, U>, M> limits)
      : m_n_bins(std::move(n_bins)), m_limits(std::move(limits)),
        m_bin_sizes(calc_bin_sizes()), m_array(m_array_dim()),
        m_count(m_array_dim()) {
    m_ones.fill(T{1});
  }

  /** \brief Get the number of bins for each dimension. */
  std::array<std::size_t, M> get_n_bins() const { return m_n_bins; }

  /** \brief Get the histogram data. */
  std::vector<T> get_histogram() const {
    return {m_array.data(), m_array.data() + m_array.num_elements()};
  }

  /** \brief Get the histogram count data. */
  std::vector<std::size_t> get_tot_count() const {
    return {m_count.data(), m_count.data() + m_count.num_elements()};
  }

  /** \brief Get the ranges (min, max) for each dimension. */
  std::array<std::pair<U, U>, M> get_limits() const { return m_limits; }

  /** \brief Get the bin sizes. */
  std::array<U, M> get_bin_sizes() const { return m_bin_sizes; }

  /**
   * \brief Add data to the histogram.
   * \param pos    Position to update.
   */
  void update(Span<const U> pos) { update(pos, m_ones); }

  /**
   * \brief Add data to the histogram.
   * \param pos    Position to update.
   * \param value  Value to add.
   */
  void update(Span<const U> pos, Span<const T> value) {
    if (pos.size() != M) {
      throw std::invalid_argument("Wrong dimensions for the coordinates");
    }
    if (value.size() != N) {
      throw std::invalid_argument("Wrong dimensions for the value");
    }
    if (check_limits(pos)) {
      boost::array<array_index, M + 1> index;
      for (std::size_t i = 0; i < M; ++i) {
        index[i] = calc_bin_index(pos[i], m_limits[i].first, m_bin_sizes[i]);
      }
      for (array_index i = 0; i < N; ++i) {
        index.back() = i;
        m_array(index) += value[i];
        m_count(index)++;
      }
    }
  }

  /** \brief Normalize histogram. */
  virtual void normalize() {
    auto const bin_volume = std::accumulate(
        m_bin_sizes.begin(), m_bin_sizes.end(), U{1}, std::multiplies<U>());
    std::transform(
        m_array.data(), m_array.data() + m_array.num_elements(), m_array.data(),
        [bin_volume](T v) { return static_cast<T>(v / bin_volume); });
  }

private:
  /**
   * \brief Calculate the bin index.
   * \param value  Position on that dimension.
   * \param offset Bin offset on that dimension.
   * \param size   Bin size on that dimension.
   */
  array_index calc_bin_index(double value, double offset, double size) const {
    return static_cast<array_index>(std::floor((value - offset) / size));
  }

  /**
   * \brief Calculate the bin sizes.
   */
  std::array<U, M> calc_bin_sizes() const {
    std::array<U, M> bin_sizes;
    for (std::size_t i = 0; i < M; ++i) {
      bin_sizes[i] = (m_limits[i].second - m_limits[i].first) /
                     static_cast<U>(m_n_bins[i]);
    }
    return bin_sizes;
  }

  /**
   * \brief Check if the position lies within the histogram limits.
   * \param pos     Position to check.
   */
  bool check_limits(Span<const U> pos) const {
    assert(pos.size() == M);
    bool within_range = true;
    for (std::size_t i = 0; i < M; ++i) {
      if (pos[i] < m_limits[i].first or pos[i] >= m_limits[i].second)
        within_range = false;
    }
    return within_range;
  }

  std::array<std::size_t, M + 1> m_array_dim() const {
    std::array<std::size_t, M + 1> dimensions;
    std::copy(m_n_bins.begin(), m_n_bins.end(), dimensions.begin());
    dimensions.back() = N;
    return dimensions;
  }

protected:
  /// Number of bins for each dimension.
  std::array<std::size_t, M> m_n_bins;
  /// Min and max values for each dimension.
  std::array<std::pair<U, U>, M> m_limits;
  /// Bin sizes for each dimension.
  std::array<U, M> m_bin_sizes;
  /// Histogram data.
  array_type m_array;
  /// Track the number of total hits per bin entry.
  count_type m_count;
  std::array<T, N> m_ones;
};

/**
 * \brief Histogram in cylindrical coordinates.
 * \tparam T     Histogram data type.
 * \tparam N     Histogram data dimensionality.
 * \tparam M     Coordinates data dimensionality.
 * \tparam U     Coordinates data type.
 */
template <typename T, std::size_t N, std::size_t M = 3, typename U = double>
class CylindricalHistogram : public Histogram<T, N, M, U> {
  using Histogram<T, N, M, U>::m_n_bins;
  using Histogram<T, N, M, U>::m_limits;
  using Histogram<T, N, M, U>::m_bin_sizes;
  using Histogram<T, N, M, U>::m_array;

public:
  using Histogram<T, N, M, U>::Histogram;

  void normalize() override {
    auto const min_r = m_limits[0].first;
    auto const r_bin_size = m_bin_sizes[0];
    auto const phi_bin_size = m_bin_sizes[1];
    auto const z_bin_size = m_bin_sizes[2];
    auto const n_bins_r = m_n_bins[0];
    for (std::size_t i = 0; i < n_bins_r; i++) {
      auto const r_left = min_r + static_cast<U>(i) * r_bin_size;
      auto const r_right = r_left + r_bin_size;
      auto const bin_volume =
          (r_right * r_right - r_left * r_left) * z_bin_size * phi_bin_size / 2;
      auto *begin = m_array[i].origin();
      std::transform(
          begin, begin + m_array[i].num_elements(), begin,
          [bin_volume](T v) { return static_cast<T>(v / bin_volume); });
    }
  }
};

} // Namespace Utils

#endif

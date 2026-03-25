/*
 * Copyright (C) 2016-2026 The ESPResSo project
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

#pragma once

#include <boost/multi_array.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <functional>
#include <numeric>
#include <span>
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

protected:
  using array_index = typename array_type::index;

public:
  /**
   * \brief Histogram constructor.
   * \param n_bins  the number of bins in each histogram dimension.
   * \param limits  the minimum/maximum data values to consider for the
   *        histogram.
   */
  Histogram(std::array<std::size_t, M> n_bins,
            std::array<std::pair<U, U>, M> limits)
      : m_n_bins(std::move(n_bins)), m_limits(std::move(limits)),
        m_bin_sizes(calc_bin_sizes()), m_array(m_array_dim()),
        m_count(m_array_dim()) {
    m_ones.fill(T{1});
  }

  virtual ~Histogram() = default;

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
  void update(std::span<const U> pos) { update(pos, m_ones); }

  /**
   * \brief Add data to the histogram.
   * \param pos    Position to update.
   * \param value  Value to add.
   */
  void update(std::span<const U> pos, std::span<const T> value) {
    if (pos.size() != M) {
      throw std::invalid_argument("Wrong dimensions for the coordinates");
    }
    if (value.size() != N) {
      throw std::invalid_argument("Wrong dimensions for the value");
    }
    if (check_limits(pos)) {
      auto index = calc_bin_index(pos);
      for (std::size_t i = 0; i < N; ++i) {
        index.back() = static_cast<array_index>(i);
        m_array(index) += value[i];
        m_count(index)++;
      }
    }
  }

  /** \brief Normalize histogram. */
  virtual void normalize() {
    auto const bin_volume = std::accumulate(
        m_bin_sizes.begin(), m_bin_sizes.end(), U{1}, std::multiplies<U>());
    std::ranges::transform(
        std::span(m_array.data(), m_array.num_elements()), m_array.data(),
        [bin_volume](T v) { return static_cast<T>(v / bin_volume); });
  }

private:
  /**
   * \brief Calculate the bin index.
   * \param pos  Position.
   */
  auto calc_bin_index(std::span<const U> const &pos) const {
    boost::array<array_index, M + 1> index;
    for (std::size_t i = 0; i < M; ++i) {
      auto const offset = m_limits[i].first;
      auto const size = m_bin_sizes[i];
      auto const n_bins = static_cast<long>(m_n_bins[i]);
      auto const bin = static_cast<long>(std::floor((pos[i] - offset) / size));
      // handle edge cases when the position is exactly between two bins:
      // due to precision loss in the offset subtraction, the bin index might
      // be off by one, so we fold it here back inside the valid range
      index[i] = static_cast<array_index>(std::clamp(bin, 0l, n_bins - 1l));
    }
    return index;
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
  bool check_limits(std::span<const U> const &pos) const {
    assert(pos.size() == M);
    auto it_limits = m_limits.begin();
    return std::ranges::all_of(pos, [&it_limits](U const value) {
      auto const [lower, upper] = *it_limits;
      ++it_limits;
      return value >= lower and value < upper;
    });
  }

  std::array<std::size_t, M + 1> m_array_dim() const {
    std::array<std::size_t, M + 1> dimensions;
    std::ranges::copy(m_n_bins, dimensions.begin());
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
  using Base = Histogram<T, N, M, U>;
  using Base::m_array;
  using Base::m_bin_sizes;
  using Base::m_limits;
  using Base::m_n_bins;
  using typename Base::array_index;

public:
  using Base::Histogram;

  void normalize() override {
    auto const min_r = m_limits[0].first;
    auto const r_bin_size = m_bin_sizes[0];
    auto const phi_bin_size = m_bin_sizes[1];
    auto const z_bin_size = m_bin_sizes[2];
    auto const n_bins_r = static_cast<array_index>(m_n_bins[0]);
    for (array_index i = 0; i < n_bins_r; i++) {
      auto const r_left = min_r + static_cast<U>(i) * r_bin_size;
      auto const r_right = r_left + r_bin_size;
      auto const bin_volume = (r_right * r_right - r_left * r_left) *
                              z_bin_size * phi_bin_size / U(2);
      auto *begin = m_array[i].origin();
      std::ranges::transform(
          std::span(begin, m_array[i].num_elements()), begin,
          [bin_volume](T v) { return static_cast<T>(v / bin_volume); });
    }
  }
};

} // Namespace Utils

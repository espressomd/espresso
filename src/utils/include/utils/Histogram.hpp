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
#include "utils/constants.hpp"
#include "utils/index.hpp"

#include <algorithm>
#include <array>
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
 * \tparam Dims  Histogram dimensionality.
 * \tparam T     Histogram data type.
 * \tparam U     Coordinates data type.
 */
template <typename T, std::size_t Dims, typename U = double> class Histogram {
public:
  /**
   * \brief Histogram constructor.
   * \param n_bins  the number of bins in each histogram dimension.
   * \param n_dims_data  the number of dimensions the data has (e.g. 3 for
   *        vector field).
   * \param limits  the minimum/maximum data values to consider for the
   *        histogram.
   */
  explicit Histogram(std::array<std::size_t, Dims> n_bins,
                     std::size_t n_dims_data,
                     std::array<std::pair<U, U>, Dims> limits)
      : m_n_bins(n_bins), m_limits(limits), m_n_dims_data(n_dims_data),
        m_ones(n_dims_data, T{1}) {
    if (n_bins.size() != limits.size()) {
      throw std::invalid_argument("Argument for number of bins and limits do "
                                  "not have same number of dimensions!");
    }
    m_bin_sizes = calc_bin_sizes(limits, n_bins);
    std::size_t n_bins_total =
        m_n_dims_data * std::accumulate(std::begin(n_bins), std::end(n_bins), 1,
                                        std::multiplies<std::size_t>());
    m_hist = std::vector<T>(n_bins_total);
    m_tot_count = std::vector<std::size_t>(n_bins_total);
  }

  /** \brief Get the number of bins for each dimension. */
  std::array<std::size_t, Dims> get_n_bins() const { return m_n_bins; }

  /** \brief Get the histogram data. */
  std::vector<T> get_histogram() const { return m_hist; }

  /** \brief Get the histogram count data. */
  std::vector<std::size_t> get_tot_count() const { return m_tot_count; }

  /** \brief Get the ranges (min, max) for each dimension. */
  std::array<std::pair<U, U>, Dims> get_limits() const { return m_limits; }

  /** \brief Get the bin sizes. */
  std::array<U, Dims> get_bin_sizes() const { return m_bin_sizes; }

  /**
   * \brief Add data to the histogram.
   * \param pos  Position to increment.
   */
  void update(Span<const U> pos) {
    if (check_limits(pos, m_limits)) {
      update(pos, m_ones);
    }
  }

  /**
   * \brief Add data to the histogram.
   * \param pos    Position to update.
   * \param value  Value to add.
   */
  void update(Span<const U> pos, Span<const T> value) {
    if (check_limits(pos, m_limits)) {
      Array<std::size_t, Dims> index;
      for (std::size_t dim = 0; dim < m_n_bins.size(); ++dim) {
        index[dim] =
            calc_bin_index(pos[dim], m_limits[dim].first, m_bin_sizes[dim]);
      }
      auto const flat_index =
          m_n_dims_data * ::Utils::ravel_index(index, m_n_bins);
      if (value.size() != m_n_dims_data)
        throw std::invalid_argument("Wrong dimensions for the given value!");
      for (std::size_t ind = 0; ind < m_n_dims_data; ++ind) {
        m_hist[flat_index + ind] += value[ind];
        m_tot_count[flat_index + ind]++;
      }
    }
  }

  /** \brief Histogram normalization. */
  virtual void normalize() {
    auto const bin_volume = std::accumulate(
        m_bin_sizes.begin(), m_bin_sizes.end(), U{1}, std::multiplies<U>());
    std::transform(
        m_hist.begin(), m_hist.end(), m_hist.begin(),
        [bin_volume](T v) { return static_cast<T>(v / bin_volume); });
  }

private:
  /**
   * \brief Calculate the bin index.
   * \param value  Position on that dimension.
   * \param offset Bin offset on that dimension.
   * \param size   Bin size on that dimension.
   */
  std::size_t calc_bin_index(double value, double offset, double size) const {
    return static_cast<std::size_t>(std::floor((value - offset) / size));
  }

  /**
   * \brief Calculate the bin sizes.
   * \param limits  contains min/max values for each dimension.
   * \param n_bins  number of bins for each dimension.
   * \return The bin sizes for each dimension.
   */
  std::array<U, Dims>
  calc_bin_sizes(std::array<std::pair<U, U>, Dims> const &limits,
                 std::array<std::size_t, Dims> const &n_bins) const {
    std::array<U, Dims> tmp;
    for (std::size_t ind = 0; ind < Dims; ++ind) {
      tmp[ind] = static_cast<U>((limits[ind].second - limits[ind].first) /
                                n_bins[ind]);
    }
    return tmp;
  }

  /**
   * \brief Check if data is within limits.
   * \param pos     Position to check.
   * \param limits  the min/max values.
   */
  bool check_limits(Span<const U> pos,
                    std::array<std::pair<U, U>, Dims> limits) const {
    if (pos.size() != limits.size()) {
      throw std::invalid_argument("Dimension of pos and limits not the same!");
    }
    bool within_range = true;
    for (std::size_t i = 0; i < pos.size(); ++i) {
      if (pos[i] < limits[i].first or pos[i] >= limits[i].second)
        within_range = false;
    }
    return within_range;
  }

private:
  /// Number of bins for each dimension.
  std::array<std::size_t, Dims> m_n_bins;
  /// Min and max values for each dimension.
  std::array<std::pair<U, U>, Dims> m_limits;
  /// Bin sizes for each dimension.
  std::array<U, Dims> m_bin_sizes;

protected:
  /// Flat histogram data.
  std::vector<T> m_hist;
  /// Number of dimensions for a single data point.
  std::size_t m_n_dims_data;
  /// Track the number of total hits per bin entry.
  std::vector<std::size_t> m_tot_count;
  std::vector<T> m_ones;
};

/**
 * \brief Histogram in cylindrical coordinates.
 * \tparam Dims  Histogram dimensionality.
 * \tparam T     Histogram data type.
 * \tparam U     Coordinates data type.
 */
template <typename T, std::size_t Dims, typename U = double>
class CylindricalHistogram : public Histogram<T, Dims, U> {
public:
  using Histogram<T, Dims, U>::Histogram;
  using Histogram<T, Dims, U>::get_n_bins;
  using Histogram<T, Dims, U>::get_limits;
  using Histogram<T, Dims, U>::get_bin_sizes;
  using Histogram<T, Dims, U>::m_hist;
  using Histogram<T, Dims, U>::m_n_dims_data;

  void normalize() override {
    std::array<std::size_t, 4> unravelled_index;
    auto const dims = get_n_bins();
    std::array<std::size_t, 4> extended_dims;
    std::copy(dims.begin(), dims.end(), extended_dims.begin());
    extended_dims[3] = m_n_dims_data;
    for (std::size_t ind = 0; ind < m_hist.size(); ind += m_n_dims_data) {
      // Get the unraveled indices and calculate the bin volume.
      ::Utils::unravel_index(extended_dims.begin(), extended_dims.end(),
                             unravelled_index.begin(), ind);
      auto const r_bin = static_cast<U>(unravelled_index[0]);
      auto const min_r = get_limits()[0].first;
      auto const r_bin_size = get_bin_sizes()[0];
      auto const phi_bin_size = get_bin_sizes()[1];
      auto const z_bin_size = get_bin_sizes()[2];
      auto const bin_volume =
          Utils::pi<U>() *
          ((min_r + (r_bin + 1) * r_bin_size) *
               (min_r + (r_bin + 1) * r_bin_size) -
           (min_r + r_bin * r_bin_size) * (min_r + r_bin * r_bin_size)) *
          z_bin_size * phi_bin_size / (2 * Utils::pi<U>());
      for (std::size_t dim = 0; dim < m_n_dims_data; ++dim) {
        m_hist[ind + dim] /= bin_volume;
      }
    }
  }
};

} // Namespace Utils

#endif

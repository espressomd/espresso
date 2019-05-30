/*
  Copyright (C) 2016-2018 The ESPResSo project

  This file is part of ESPResSo.


  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#ifndef UTILS_HISTOGRAM_HPP
#define UTILS_HISTOGRAM_HPP

#include <algorithm>
#include <cmath>
#include <functional>
#include <numeric>
#include <vector>

#include "Span.hpp"
#include "constants.hpp"
#include "utils/index.hpp"

namespace Utils {

inline size_t calculate_bin_index(double value, double bin_size,
                                  double offset) {
  return std::floor((value - offset) / bin_size);
}

/**
 * \brief Calculate the bin sizes.
 * \param limits  contains min/max values for each dimension.
 * \param n_bins  number of bins for each dimension.
 * \return The bin sizes for each dimension.
 */
template <typename T, size_t Dims>
std::array<T, Dims>
calc_bin_sizes(std::array<std::pair<T, T>, Dims> const &limits,
               std::array<size_t, Dims> const &n_bins) {
  std::array<T, Dims> tmp;
  for (size_t ind = 0; ind < Dims; ++ind) {
    tmp[ind] = (limits[ind].second - limits[ind].first) / n_bins[ind];
  }
  return tmp;
}

/*
 * \brief Check if data is within limits.
 * \param data  data value to check.
 * \param limits  the min/max values.
 */
template <typename T, size_t Dims>
inline bool check_limits(Span<const T> data,
                         std::array<std::pair<T, T>, Dims> limits) {
  if (data.size() != limits.size()) {
    throw std::invalid_argument("Dimension of data and limits not the same!");
  }
  bool within_range = true;
  for (size_t i = 0; i < data.size(); ++i) {
    if (data[i] < limits[i].first or data[i] >= limits[i].second)
      within_range = false;
  }
  return within_range;
}

template <typename T, size_t Dims> class Histogram {
public:
  explicit Histogram(std::array<size_t, Dims> n_bins, size_t n_dims_data,
                     std::array<std::pair<T, T>, Dims> limits);
  std::array<size_t, Dims> get_n_bins() const;
  std::vector<T> get_histogram() const;
  std::vector<size_t> get_tot_count() const;
  std::array<std::pair<T, T>, Dims> get_limits() const;
  std::array<T, Dims> get_bin_sizes() const;
  void update(Span<const T> data);
  void update(Span<const T> data, Span<const T> weights);
  void normalize();

private:
  // Number of bins for each dimension.
  std::array<size_t, Dims> m_n_bins;
  // Min and max values for each dimension.
  std::array<std::pair<T, T>, Dims> m_limits;
  // Bin sizes for each dimension.
  std::array<T, Dims> m_bin_sizes;
  virtual void do_normalize();

protected:
  // Flat histogram data.
  std::vector<T> m_hist;
  // Number of dimensions for a single data point.
  size_t m_n_dims_data;
  // Track the number of total hits per bin entry.
  std::vector<size_t> m_tot_count;
  std::vector<T> m_ones;
};

/**
 * \brief Histogram constructor.
 * \param n_bins  the number of bins in each histogram dimension.
 * \param n_dims_data  the number of dimensions the data has (e.g. 3 for
 *        vector field).
 * \param limits  the minimum/maximum data values to consider for the
 *        histogram.
 */
template <typename T, size_t Dims>
Histogram<T, Dims>::Histogram(std::array<size_t, Dims> n_bins,
                              size_t n_dims_data,
                              std::array<std::pair<T, T>, Dims> limits)
    : m_n_bins(n_bins), m_limits(limits), m_n_dims_data(n_dims_data),
      m_ones(n_dims_data, T{1.}) {
  if (n_bins.size() != limits.size()) {
    throw std::invalid_argument("Argument for number of bins and limits do "
                                "not have same number of dimensions!");
  }
  m_bin_sizes = calc_bin_sizes(limits, n_bins);
  size_t n_bins_total =
      m_n_dims_data * std::accumulate(std::begin(n_bins), std::end(n_bins), 1,
                                      std::multiplies<size_t>());
  m_hist = std::vector<T>(n_bins_total);
  m_tot_count = std::vector<size_t>(n_bins_total);
}

/**
 * \brief Add data to the histogram.
 * \param data  vector of single data value with type T.
 *              The size of the given vector has to match the number
 *              of dimensions of the histogram.
 */
template <typename T, size_t Dims>
void Histogram<T, Dims>::update(Span<const T> data) {
  if (check_limits(data, m_limits)) {
    update(data, m_ones);
  }
}

/**
 * \brief Add data to the histogram.
 * \param data  vector of single data value with type T.
 *              The size of the given vector has to match the number
 *              of dimensions of the histogram.
 * \param weights  m_n_dims_data dimensional weights.
 */
template <typename T, size_t Dims>
void Histogram<T, Dims>::update(Span<const T> data, Span<const T> weights) {
  if (check_limits(data, m_limits)) {
    Array<size_t, Dims> index;
    for (size_t dim = 0; dim < m_n_bins.size(); ++dim) {
      index[dim] =
          calculate_bin_index(data[dim], m_bin_sizes[dim], m_limits[dim].first);
    }
    auto const flat_index =
        m_n_dims_data * ::Utils::ravel_index(index, m_n_bins);
    if (weights.size() != m_n_dims_data)
      throw std::invalid_argument("Wrong dimensions of given weights!");
    for (size_t ind = 0; ind < m_n_dims_data; ++ind) {
      m_hist[flat_index + ind] += weights[ind];
      m_tot_count[flat_index + ind] += 1;
    }
  }
}

/**
 * \brief Get the bin sizes.
 */
template <typename T, size_t Dims>
std::array<T, Dims> Histogram<T, Dims>::get_bin_sizes() const {
  return m_bin_sizes;
}

/**
 * \brief Get the number of bins for each dimension.
 */
template <typename T, size_t Dims>
std::array<size_t, Dims> Histogram<T, Dims>::get_n_bins() const {
  return m_n_bins;
}

/**
 * \brief Get the ranges (min, max) for each dimension.
 */
template <typename T, size_t Dims>
std::array<std::pair<T, T>, Dims> Histogram<T, Dims>::get_limits() const {
  return m_limits;
}

/**
 * \brief Get the histogram data.
 */
template <typename T, size_t Dims>
std::vector<T> Histogram<T, Dims>::get_histogram() const {
  return m_hist;
}

/**
 * \brief Get the histogram count data.
 */
template <typename T, size_t Dims>
std::vector<size_t> Histogram<T, Dims>::get_tot_count() const {
  return m_tot_count;
}
/**
 * \brief Histogram normalization. (private member function can be overridden by
 * subclasses).
 */
template <typename T, size_t Dims> void Histogram<T, Dims>::normalize() {
  do_normalize();
}

/**
 * \brief Histogram normalization.
 *        Divide by the bin volume.
 */
template <typename T, size_t Dims> void Histogram<T, Dims>::do_normalize() {
  T bin_volume = std::accumulate(m_bin_sizes.begin(), m_bin_sizes.end(),
                                 static_cast<T>(1.0), std::multiplies<T>());
  std::transform(m_hist.begin(), m_hist.end(), m_hist.begin(),
                 [bin_volume](T v) { return v / bin_volume; });
}

template <typename T, size_t Dims>
class CylindricalHistogram : public Histogram<T, Dims> {
public:
  using Histogram<T, Dims>::Histogram;
  using Histogram<T, Dims>::get_n_bins;
  using Histogram<T, Dims>::get_limits;
  using Histogram<T, Dims>::get_bin_sizes;
  using Histogram<T, Dims>::m_hist;
  using Histogram<T, Dims>::m_n_dims_data;

private:
  void do_normalize() override {
    std::array<std::size_t, 4> unravelled_index;
    int r_bin;
    double min_r, r_bin_size, phi_bin_size, z_bin_size, bin_volume;
    auto const dims = get_n_bins();
    std::array<std::size_t, 4> extended_dims;
    std::copy(dims.begin(), dims.end(), extended_dims.begin());
    extended_dims[3] = m_n_dims_data;
    for (size_t ind = 0; ind < m_hist.size(); ind += m_n_dims_data) {
      // Get the unraveled indices and calculate the bin volume.
      ::Utils::unravel_index(extended_dims.begin(), extended_dims.end(),
                             unravelled_index.begin(), ind);
      r_bin = unravelled_index[0];
      min_r = get_limits()[0].first;
      r_bin_size = get_bin_sizes()[0];
      phi_bin_size = get_bin_sizes()[1];
      z_bin_size = get_bin_sizes()[2];
      bin_volume =
          Utils::pi() *
          ((min_r + (r_bin + 1) * r_bin_size) *
               (min_r + (r_bin + 1) * r_bin_size) -
           (min_r + r_bin * r_bin_size) * (min_r + r_bin * r_bin_size)) *
          z_bin_size * phi_bin_size / (2 * Utils::pi());
      for (size_t dim = 0; dim < m_n_dims_data; ++dim) {
        m_hist[ind + dim] /= bin_volume;
      }
    }
  }
};

} // Namespace Utils

#endif

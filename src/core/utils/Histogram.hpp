/*
  Copyright (C) 2016,2017 The ESPResSo project

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

#include "constants.hpp"

#include <algorithm>
#include <functional>
#include <numeric>
#include <vector>
#include <cmath>

namespace Utils {

inline size_t calculate_bin_index(double value, double bin_size, double offset) {
  return std::floor((value - offset) / bin_size);
}

inline size_t ravel_index(std::vector<size_t> unravelled_indices,
                          std::vector<size_t> n_bins) {
  // index calculation: using the following formula for N dimensions:
  //   ind = ind_{N-1} + sum_{j=0}^{N-2} (ind_j * prod_{k=j+1}^{N-1} n_k)
  size_t res = unravelled_indices.back();
  for (size_t j = 0; j < unravelled_indices.size() - 1; ++j) {
    res += unravelled_indices[j] * std::accumulate(n_bins.begin() + j + 1,
                                                   n_bins.end(), 1,
                                                   std::multiplies<size_t>());
  }
  return res;
}

/**
 * \brief Returns the unravelled index of the provided flat index.
 *        Therefore is the inversion of flattening an ndims dimensional index.
 * \param len_dims an int array of length ndims containing the lengths of the
 * dimensions. (Input) \param ndims int denoting the number of dimensions.
 * (Input) \flattened_index an int denoting the flat index. (Input)
 * \unravelled_index_out an int array with length ndims where the unflat indices
 * are written to. (Output)
 */
inline void unravel_index(const int *const len_dims, const int ndims,
                          const int flattened_index,
                          int *unravelled_index_out) {
  // idea taken from
  // http://codinghighway.com/2014/02/22/c-multi-dimensional-arrays-part-2-flattened-to-unflattened-index/
  std::vector<int> mul(ndims);
  mul[ndims - 1] = 1;
  for (int j = ndims - 2; j >= 0; j--)
    mul[j] = mul[j + 1] * len_dims[j + 1];
  for (int j = 0; j < ndims; j++)
    unravelled_index_out[j] = (flattened_index / mul[j]) % len_dims[j];
}

/**
 * \brief Calculate the bin sizes.
 * \param limits: contains min/max values for each dimension.
 * \param n_bins: number of bins for each dimension.
 * \return The bin sizes for each dimension.
 */
template <typename T>
std::vector<T> calc_bin_sizes(std::vector<std::pair<T, T>> const &limits,
                              std::vector<size_t> const &n_bins) {
  std::vector<T> tmp;
  for (size_t ind = 0; ind < limits.size(); ++ind) {
    tmp.push_back((limits[ind].second - limits[ind].first) / n_bins[ind]);
  }
  return tmp;
}

/*
 * \brief Check if data is within limits.
 * \param data: data value to check.
 * \param limits: the min/max values.
 */
template <typename T>
inline bool check_limits(std::vector<T> const &data,
                         std::vector<std::pair<T, T>> limits) {
  if (data.size() != limits.size()) {
    throw std::invalid_argument("Dimension of data and limits not the same!");
  }
  bool within_range = true;
  for (size_t i = 0; i < data.size(); ++i) {
    if (data[i] < limits[i].first or data[i] > limits[i].second)
      within_range = false;
  }
  return within_range;
}

template <typename T> class Histogram {
public:
  explicit Histogram(std::vector<size_t> n_bins, size_t n_dims_data,
                     std::vector<std::pair<T, T>> limits);
  std::vector<size_t> get_n_bins() const;
  std::vector<T> get_histogram() const;
  std::vector<size_t> get_tot_count() const;
  std::vector<std::pair<T, T>> get_limits() const;
  std::vector<T> get_bin_sizes() const;
  void update(std::vector<T> const &data);
  void update(std::vector<T> const &data, std::vector<T> const &weights);
  void normalize();

private:
  // Number of bins for each dimension.
  std::vector<size_t> m_n_bins;
  // Min and max values for each dimension.
  std::vector<std::pair<T, T>> m_limits;
  // Bin sizes for each dimension.
  std::vector<T> m_bin_sizes;
  virtual void do_normalize();

protected:
  // Flat histogram data.
  std::vector<T> m_hist;
  // Number of dimensions for a single data point.
  size_t m_n_dims_data;
  // Track the number of total hits per bin entry.
  std::vector<size_t> m_tot_count;
};

/**
 * \brief Histogram constructor.
 * \param n_bins: the number of bins in each histogram dimension.
 * \param n_dims_data: the number of dimensions the data has (e.g. 3 for
 *        vector field).
 * \param limits: the minimum/maximum data values to consider for the
 *        histogram.
 */
template <typename T>
Histogram<T>::Histogram(std::vector<size_t> n_bins, size_t n_dims_data,
                        std::vector<std::pair<T, T>> limits)
    : m_n_bins(n_bins), m_limits(limits), m_n_dims_data(n_dims_data) {
  if (n_bins.size() != limits.size()) {
    throw std::invalid_argument("Argument for number of bins and limits do "
                                "not have same number of dimensions!");
  }
  m_bin_sizes = calc_bin_sizes<T>(limits, n_bins);
  size_t n_bins_total =
      m_n_dims_data * std::accumulate(std::begin(n_bins), std::end(n_bins), 1,
                                      std::multiplies<size_t>());
  m_hist = std::vector<T>(n_bins_total);
  m_tot_count = std::vector<size_t>(n_bins_total);
}

/**
 * \brief Add data to the histogram.
 * \param data: vector of single data value with type T.
 *              The size of the given vector has to match the number
 *              of dimensions of the histogram.
 */
template <typename T> void Histogram<T>::update(std::vector<T> const &data) {
  if (check_limits(data, m_limits)) {
    std::vector<T> weights(m_n_dims_data, static_cast<T>(1.0));
    update(data, weights);
  }
}

/**
 * \brief Add data to the histogram.
 * \param data: vector of single data value with type T.
 *              The size of the given vector has to match the number
 *              of dimensions of the histogram.
 * \param weights: m_n_dims_data dimensional weights.
 */
template <typename T>
void Histogram<T>::update(std::vector<T> const &data,
                          std::vector<T> const &weights) {
  if (check_limits(data, m_limits)) {
    std::vector<size_t> index;
    for (size_t dim = 0; dim < m_n_bins.size(); ++dim) {
      index.push_back(calculate_bin_index(data[dim], m_bin_sizes[dim],
                                          m_limits[dim].first));
    }
    size_t flat_index = m_n_dims_data * ::Utils::ravel_index(index, m_n_bins);
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
template <typename T> std::vector<T> Histogram<T>::get_bin_sizes() const {
  return m_bin_sizes;
}

/**
 * \brief Get the number of bins for each dimension.
 */
template <typename T> std::vector<size_t> Histogram<T>::get_n_bins() const {
  return m_n_bins;
}

/**
 * \brief Get the ranges (min, max) for each dimension.
 */
template <typename T>
std::vector<std::pair<T, T>> Histogram<T>::get_limits() const {
  return m_limits;
}

/**
 * \brief Get the histogram data.
 */
template <typename T> std::vector<T> Histogram<T>::get_histogram() const {
  return m_hist;
}

/**
 * \brief Get the histogram count data.
 */
template <typename T> std::vector<size_t> Histogram<T>::get_tot_count() const {
  return m_tot_count;
}
/**
 * \brief Histogram normalization. (private member function can be overridden by
 * subclasses).
 */
template <typename T> void Histogram<T>::normalize() { do_normalize(); }

/**
 * \brief Histogram normalization.
 *        Divide by the bin volume.
 */
template <typename T> void Histogram<T>::do_normalize() {
  T bin_volume = std::accumulate(m_bin_sizes.begin(), m_bin_sizes.end(),
                                 static_cast<T>(1.0), std::multiplies<T>());
  std::transform(m_hist.begin(), m_hist.end(), m_hist.begin(),
                 [this, bin_volume](T v) { return v / bin_volume;});
}


template <typename T>
class CylindricalHistogram : public Histogram<T> {
public:
  using Histogram<T>::Histogram;
  using Histogram<T>::get_n_bins;
  using Histogram<T>::get_limits;
  using Histogram<T>::get_bin_sizes;
  using Histogram<T>::m_hist;
  using Histogram<T>::m_n_dims_data;

private:
  void do_normalize() override {
    int unravelled_index[4];
    int r_bin;
    double min_r, r_bin_size, phi_bin_size, z_bin_size, bin_volume;
    // Ugly vector cast due to "unravel_index" function.
    std::vector<size_t> len_bins_u = get_n_bins();
    std::vector<int> len_bins(len_bins_u.begin(), len_bins_u.end());
    len_bins.push_back(m_n_dims_data);
    for (size_t ind = 0; ind < m_hist.size(); ind += m_n_dims_data) {
      // Get the unravelled indices and calculate the bin volume.
      ::Utils::unravel_index(len_bins.data(), 4, ind, unravelled_index);
      r_bin = unravelled_index[0];
      min_r = get_limits()[0].first;
      r_bin_size = get_bin_sizes()[0];
      phi_bin_size = get_bin_sizes()[1];
      z_bin_size = get_bin_sizes()[2];
      bin_volume =
          PI * ((min_r + (r_bin + 1) * r_bin_size) *
                    (min_r + (r_bin + 1) * r_bin_size) -
                (min_r + r_bin * r_bin_size) * (min_r + r_bin * r_bin_size)) *
          z_bin_size * phi_bin_size / (2 * PI);
      for (size_t dim = 0; dim < m_n_dims_data; ++dim) {
        m_hist[ind + dim] /= bin_volume;
      }
    }
  }
};

} // Namespace Utils

#endif

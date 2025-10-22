/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef CORE_UTILS_ACCUMULATOR
#define CORE_UTILS_ACCUMULATOR

#include <boost/serialization/access.hpp>
#include <boost/serialization/vector.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Utils {

template <typename T> struct AccumulatorData {
  T mean = T{};
  T m = T{};

private:
  // Allow serialization to access non-public data members.
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned /*version*/) {
    ar & mean & m;
  }
};

class Accumulator {
public:
  explicit Accumulator(std::size_t N) : m_n(0u), m_acc_data(N) {}
  void operator()(const std::vector<double> &);
  std::vector<double> mean() const;
  std::vector<double> variance() const;
  std::vector<double> std_error() const;

private:
  std::size_t m_n;
  std::vector<AccumulatorData<double>> m_acc_data;
  // Allow serialization to access non-public data members.
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned /*version*/) {
    ar & m_n & m_acc_data;
  }
};

inline void Accumulator::operator()(const std::vector<double> &data) {
  if (data.size() != m_acc_data.size())
    throw std::runtime_error(
        "The given data size does not fit the initialized size!");
  ++m_n;
  if (m_n == 1u) {
    std::ranges::transform(
        data, m_acc_data.begin(),
        [](double d) -> AccumulatorData<double> { return {d, 0.0}; });
  } else {
    auto const denominator = static_cast<double>(m_n);
    std::ranges::transform(m_acc_data, data, m_acc_data.begin(),
                           [denominator](auto const &a, double d) {
                             auto const old_mean = a.mean;
                             auto const new_mean =
                                 old_mean + (d - old_mean) / denominator;
                             auto const new_m =
                                 a.m + (d - old_mean) * (d - new_mean);
                             return AccumulatorData<double>{new_mean, new_m};
                           });
  }
}

/**
 * @brief Compute the sample mean.
 */
inline std::vector<double> Accumulator::mean() const {
  std::vector<double> res{};
  std::ranges::transform(m_acc_data, std::back_inserter(res),
                         [](auto const &acc_data) { return acc_data.mean; });
  return res;
}

/**
 * @brief Compute the Bessel-corrected sample variance,
 * assuming uncorrelated data.
 */
inline std::vector<double> Accumulator::variance() const {
  std::vector<double> res{};
  if (m_n == 1u) {
    res = std::vector<double>(m_acc_data.size(),
                              std::numeric_limits<double>::max());
  } else {
    auto const denominator = static_cast<double>(m_n) - 1.;
    std::ranges::transform(m_acc_data, std::back_inserter(res),
                           [denominator](auto const &acc_data) {
                             return acc_data.m / denominator;
                           });
  }
  return res;
}

/**
 * @brief Compute the standard error of the mean, assuming uncorrelated data.
 */
inline std::vector<double> Accumulator::std_error() const {
  auto const var = variance();
  std::vector<double> err{};
  auto const denominator = static_cast<double>(m_n);
  std::ranges::transform(var, std::back_inserter(err), [denominator](double d) {
    return std::sqrt(d / denominator);
  });
  return err;
}

} // namespace Utils

#endif

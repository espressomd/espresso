/*
Copyright (C) 2010-2018 The ESPResSo project

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
#ifndef CORE_UTILS_ACCUMULATOR
#define CORE_UTILS_ACCUMULATOR

#include <boost/serialization/access.hpp>

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

namespace Utils {

template <typename T> struct AccumulatorData {
  AccumulatorData() = default;
  T mean;
  T m;

private:
  // Allow serialization to access non-public data members.
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned /*version*/) {
    ar &mean &m;
  }
};

class Accumulator {
public:
  explicit Accumulator(std::size_t N) : m_n(0), m_acc_data(N) {}
  void operator()(const std::vector<double> &);
  std::vector<double> get_mean() const;
  std::vector<double> get_variance() const;
  std::vector<double> get_std_error() const;

private:
  std::size_t m_n;
  std::vector<AccumulatorData<double>> m_acc_data;
  // Allow serialization to access non-public data members.
  friend class boost::serialization::access;

  template <typename Archive>
  void serialize(Archive &ar, const unsigned /*version*/) {
    ar &m_n &m_acc_data;
  }
};

inline void Accumulator::operator()(const std::vector<double> &data) {
  if (data.size() != m_acc_data.size())
    throw std::runtime_error(
        "The given data size does not fit the initialized size!");
  ++m_n;
  if (m_n == 1) {
    std::transform(data.begin(), data.end(), m_acc_data.begin(),
                   [](double d) -> AccumulatorData<double> {
                     return {d, 0.0};
                   });
  } else {
    std::transform(m_acc_data.begin(), m_acc_data.end(), data.begin(),
                   m_acc_data.begin(),
                   [this](AccumulatorData<double> &a,
                          double d) -> AccumulatorData<double> {
                     auto const old_mean = a.mean;
                     auto const new_mean = old_mean + (d - old_mean) / m_n;
                     auto const new_m = a.m + (d - old_mean) * (d - new_mean);
                     return {new_mean, new_m};
                   });
  }
}

inline std::vector<double> Accumulator::get_mean() const {
  std::vector<double> res;
  std::transform(
      m_acc_data.begin(), m_acc_data.end(), std::back_inserter(res),
      [](const AccumulatorData<double> &acc_data) { return acc_data.mean; });
  return res;
}

inline std::vector<double> Accumulator::get_variance() const {
  std::vector<double> res;
  if (m_n == 1) {
    res = std::vector<double>(m_acc_data.size(),
                              std::numeric_limits<double>::max());
  } else {
    std::transform(
        m_acc_data.begin(), m_acc_data.end(), std::back_inserter(res),
        [this](const AccumulatorData<double> &acc_data) {
          return acc_data.m /
                 (static_cast<double>(m_n) -
                  1); // numerically stable sample variance, see
                      // https://en.wikipedia.org/wiki/Algorithms_for_calculating_variance
        });
  }
  return res;
}

/**
returns the standard error of the mean of uncorrelated data. if data are
correlated the correlation time needs to be known...
*/
inline std::vector<double> Accumulator::get_std_error() const {
  auto const variance = get_variance();
  std::vector<double> std_error(variance.size());
  std::transform(variance.begin(), variance.end(), std_error.begin(),
                 [this](double d) { return std::sqrt(d / m_n); });
  return std_error;
}

} // namespace Utils

#endif

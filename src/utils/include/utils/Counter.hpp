/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef UTILS_COUNTER_HPP
#define UTILS_COUNTER_HPP

#include <boost/serialization/access.hpp>

namespace Utils {
template <typename T> class Counter {
private:
  T m_val;
  T m_initial;
  friend class boost::serialization::access;
  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &m_val;
    ar &m_initial;
  }

public:
  explicit Counter(T initial_value = T(0)) noexcept
      : m_val(initial_value), m_initial(initial_value) {}
  Counter(T initial_value, T value) noexcept
      : m_val(value), m_initial(initial_value) {}

  void increment() { ++m_val; }

  T value() const { return m_val; }
  T initial_value() const { return m_initial; }
};
} // namespace Utils
#endif // UTILS_COUNTER_HPP

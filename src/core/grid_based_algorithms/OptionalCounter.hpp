/*
 * Copyright (C) 2020 The ESPResSo project
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
#ifndef OPTIONAL_COUNTER_HPP
#define OPTIONAL_COUNTER_HPP

#include <cstdint>
#include <utility>

#include <utils/Counter.hpp>

/** Re-implementation of a boost::optional for a RNG counter.
 *
 *  Workaround for a compiler error with Clang 9.0, boost 1.71
 *  and CUDA 10.1 (see espressomd/espresso#3650).
 */
class OptionalCounter {
private:
  Utils::Counter<uint64_t> m_counter;
  bool m_initialized;

public:
  OptionalCounter() : m_counter{}, m_initialized(false) {}
  OptionalCounter(Utils::Counter<uint64_t> const &counter)
      : m_counter(counter), m_initialized(true) {}
  OptionalCounter &operator=(Utils::Counter<uint64_t> counter) {
    m_counter = std::move(counter);
    m_initialized = true;
    return *this;
  }
  template <class Archive>
  void serialize(Archive &ar, const unsigned int /* version */) {
    ar &m_counter;
    ar &m_initialized;
  }
  bool is_initialized() noexcept { return m_initialized; }
  explicit operator bool() const noexcept { return m_initialized; }
  bool operator!() const noexcept { return !m_initialized; }
  Utils::Counter<uint64_t> &operator*() { return m_counter; }
  Utils::Counter<uint64_t> *operator->() { return &m_counter; }
};

#endif

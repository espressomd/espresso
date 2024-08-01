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

#pragma once

#include <cassert>
#include <cstddef>
#include <string>
#include <vector>

// Forward declarations
namespace boost::mpi {
class communicator;
}
namespace System {
class System;
}

namespace Accumulators {

class AccumulatorBase {
public:
  AccumulatorBase(::System::System const *system, int delta_N)
      : m_system(reinterpret_cast<void const *>(system)), m_delta_N(delta_N) {}
  virtual ~AccumulatorBase() = default;

  int &delta_N() { return m_delta_N; }
  bool has_same_system_handle(::System::System const *system) const {
    return reinterpret_cast<void const *>(system) == m_system;
  }
  void override_system_handle(::System::System const *system) {
    assert(m_system == nullptr);
    m_system = reinterpret_cast<void const *>(system);
  }

  virtual void update(boost::mpi::communicator const &comm) = 0;
  /** Dimensions needed to reshape the flat array returned by the accumulator */
  virtual std::vector<std::size_t> shape() const = 0;
  /** Serialization of private members. */
  virtual std::string get_internal_state() const = 0;
  virtual void set_internal_state(std::string const &) = 0;

protected:
  void const *m_system; ///< for bookkeeping purposes
private:
  /// Number of time steps between automatic updates.
  int m_delta_N;
};
} // namespace Accumulators

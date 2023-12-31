/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#include <utils/Vector.hpp>

#include <cmath>
#include <variant>

namespace LeesEdwards {

// Protocols determining shear rate and positional offset as a function of time

/** Lees-Edwards protocol for un-altered periodic boundary conditions */
struct Off {
  double shear_velocity(double time) const { return 0.; }
  double pos_offset(double time) const { return 0.; }
};

/** Lees-Edwards protocol for linear shearing */
struct LinearShear {
  LinearShear()
      : m_initial_pos_offset{0.}, m_shear_velocity{0.}, m_time_0{0.} {}
  LinearShear(double initial_offset, double shear_velocity, double time_0)
      : m_initial_pos_offset{initial_offset},
        m_shear_velocity{shear_velocity}, m_time_0{time_0} {}
  double shear_velocity(double time) const { return m_shear_velocity; }
  double pos_offset(double time) const {
    return m_initial_pos_offset + (time - m_time_0) * m_shear_velocity;
  }
  double m_initial_pos_offset;
  double m_shear_velocity;
  double m_time_0;
};

/** Lees-Edwards protocol for oscillatory shearing */
struct OscillatoryShear {
  OscillatoryShear()
      : m_initial_pos_offset{0.}, m_amplitude{0.}, m_omega{0.}, m_time_0{0.} {}
  OscillatoryShear(double initial_offset, double amplitude, double omega,
                   double time_0)
      : m_initial_pos_offset{initial_offset},
        m_amplitude{amplitude}, m_omega{omega}, m_time_0{time_0} {}
  double pos_offset(double time) const {
    return m_initial_pos_offset +
           m_amplitude * std::sin(m_omega * (time - m_time_0));
  }
  double shear_velocity(double time) const {
    return m_omega * m_amplitude * std::cos(m_omega * (time - m_time_0));
  }
  double m_initial_pos_offset;
  double m_amplitude;
  double m_omega;
  double m_time_0;
};

/** Type which holds the currently active protocol */
using ActiveProtocol = std::variant<Off, LinearShear, OscillatoryShear>;

class PosOffsetGetter {
public:
  PosOffsetGetter(double time) : m_time{time} {}
  template <typename T> double operator()(T const &protocol) const {
    return protocol.pos_offset(m_time);
  }

private:
  double m_time;
};

inline double get_pos_offset(double time, ActiveProtocol const &protocol) {
  return std::visit(PosOffsetGetter(time), protocol);
}

/** Visitor to get shear velocity from the Lees-Edwards protocol */
class ShearVelocityGetter {
public:
  ShearVelocityGetter(double time) : m_time{time} {}
  template <typename T> double operator()(T const &protocol) const {
    return protocol.shear_velocity(m_time);
  }

private:
  double m_time;
};

/** Calculation of current velocity */
inline double get_shear_velocity(double time, ActiveProtocol const &protocol) {
  return std::visit(ShearVelocityGetter(time), protocol);
}

} // namespace LeesEdwards

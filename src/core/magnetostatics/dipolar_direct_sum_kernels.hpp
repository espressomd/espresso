/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include <config/config.hpp>

#ifdef ESPRESSO_DIPOLES

#include <utils/Vector.hpp>

#include <cmath>

/** @brief Force and torque of one pair interaction. */
struct PairForce {
  Utils::Vector3d f{};
  Utils::Vector3d torque{};
  PairForce &operator+=(PairForce const &o) {
    f += o.f;
    torque += o.torque;
    return *this;
  }
};

/** @brief Pair force of two interacting dipoles (see dipolar_direct_sum.cpp).
 */
inline PairForce pair_force(Utils::Vector3d const &d, Utils::Vector3d const &m1,
                            Utils::Vector3d const &m2) {
  auto const pe2 = m1 * d;
  auto const pe3 = m2 * d;

  auto const r2 = d.norm2();
  auto const r = std::sqrt(r2);
  auto const r5 = r2 * r2 * r;
  auto const r7 = r5 * r2;

  auto const a = 3.0 * (m1 * m2) / r5;
  auto const b = -15.0 * pe2 * pe3 / r7;

  auto const f = (a + b) * d + 3.0 * (pe3 * m1 + pe2 * m2) / r5;
  auto const r3 = r2 * r;
  auto const t =
      -vector_product(m1, m2) / r3 + 3.0 * pe3 * vector_product(m1, d) / r5;

  return PairForce{f, t};
}

/** @brief Pair potential for two interacting dipoles. */
inline double pair_potential(Utils::Vector3d const &d,
                             Utils::Vector3d const &m1,
                             Utils::Vector3d const &m2) {
  auto const r2 = d * d;
  auto const r = std::sqrt(r2);
  auto const r3 = r2 * r;
  auto const r5 = r3 * r2;

  auto const pe1 = m1 * m2;
  auto const pe2 = m1 * d;
  auto const pe3 = m2 * d;

  return pe1 / r3 - 3.0 * pe2 * pe3 / r5;
}

#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
/** @brief Dipole field contribution from a dipole @c m1 at distance @c d. */
inline Utils::Vector3d dipole_field(Utils::Vector3d const &d,
                                    Utils::Vector3d const &m1) {
  auto const r2 = d * d;
  auto const r = std::sqrt(r2);
  auto const r3 = r2 * r;
  auto const r5 = r3 * r2;
  auto const pe2 = m1 * d;

  return 3.0 * pe2 * d / r5 - m1 / r3;
}
#endif // ESPRESSO_DIPOLE_FIELD_TRACKING

#endif // ESPRESSO_DIPOLES

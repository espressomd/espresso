/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef _BONDED_COULOMB_HPP
#define _BONDED_COULOMB_HPP
/** \file
 *  Routines to calculate the bonded Coulomb potential between
 *  particle pairs.
 */

#include "config.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cmath>

/** Parameters for %Coulomb bond Potential */
struct BondedCoulomb {
  /** %Coulomb prefactor */
  double prefactor;

  double cutoff() const { return 0.; }

  static constexpr int num = 1;

  BondedCoulomb(double prefactor) { this->prefactor = prefactor; }

  boost::optional<Utils::Vector3d> force(double q1q2,
                                         Utils::Vector3d const &dx) const;
  boost::optional<double> energy(double q1q2, Utils::Vector3d const &dx) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &prefactor;
  }
};

/** Compute the bonded Coulomb pair force.
 *  @param[in]  q1q2      Product of the particle charges.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<Utils::Vector3d>
BondedCoulomb::force(double const q1q2, Utils::Vector3d const &dx) const {
#ifdef ELECTROSTATICS
  auto const dist2 = dx.norm2();
  auto const dist3 = dist2 * std::sqrt(dist2);
  auto const fac = prefactor * q1q2 / dist3;
  return fac * dx;
#else
  return Utils::Vector3d{};
#endif
}

/** Compute the bonded Coulomb pair energy.
 *  @param[in]  q1q2      Product of the particle charges.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<double>
BondedCoulomb::energy(double const q1q2, Utils::Vector3d const &dx) const {
#ifdef ELECTROSTATICS
  auto const dist = dx.norm();
  return prefactor * q1q2 / dist;
#else
  return .0;
#endif
}

#endif

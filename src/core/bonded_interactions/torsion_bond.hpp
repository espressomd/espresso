/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#ifndef TORSION_BOND_HPP
#define TORSION_BOND_HPP
/** \file
 *  Routines to calculate the torsion_bond potential between particle pairs.
 *
 *  Implementation in \ref torsion_bond.cpp.
 */

#include "config/config.hpp"

#include "Particle.hpp"
#include <utils/Vector.hpp>

#include <boost/optional.hpp>

/** Parameters for torsion_bond Potential. */
struct TorsionBond {
  /** spring constant */
  double k;

  static constexpr int num = 1;

  double cutoff() const { return 0; }

  TorsionBond(double k);

  boost::optional<Utils::Vector3d> torque(Particle const &p1,
                                          Particle const &p2) const;
  boost::optional<double> energy(Particle const &p1, Particle const &p2) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &k;
  }
};

/** Compute the torque resulting from the torsion_bond between p1 and p2.
 *  @param[in]  p1        %First particle.
 *  @param[in]  p2        %Second particle.
 */
inline boost::optional<Utils::Vector3d>
TorsionBond::torque(Particle const &p1, Particle const &p2) const {
#ifdef ROTATION
  return -k * vector_product(p1.calc_director(), p2.calc_director());
#else
  return {};
#endif
}

/** Compute the torsion_bond energy.
 *  @param[in]  p1        %First particle.
 *  @param[in]  p2        %Second particle.
 */
inline boost::optional<double> TorsionBond::energy(Particle const &p1,
                                                   Particle const &p2) const {
#ifdef ROTATION
  double const angle = std::acos(p1.calc_director() * p2.calc_director());
  return 0.5 * k * angle * angle;
#else
  return {};
#endif
}

#endif // TORSION_BOND_HPP

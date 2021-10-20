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
#ifndef ANGLE_HARMONIC_H
#define ANGLE_HARMONIC_H
/** \file
 *  Routines to calculate the angle energy or/and and force
 *  for a particle triple using the potential described in
 *  @ref bondedIA_angle_harmonic.
 */

#include "angle_common.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <tuple>

/** Parameters for three-body angular potential (harmonic). */
struct AngleHarmonicBond {
  /** bending constant */
  double bend;
  /** equilibrium angle (default is 180 degrees) */
  double phi0;

  double cutoff() const { return 0.; }

  static constexpr int num = 2;

  AngleHarmonicBond(double bend, double phi0) {
    this->bend = bend;
    this->phi0 = phi0;
  }

  std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
  forces(Utils::Vector3d const &r_mid, Utils::Vector3d const &r_left,
         Utils::Vector3d const &r_right) const;
  double energy(Utils::Vector3d const &r_mid, Utils::Vector3d const &r_left,
                Utils::Vector3d const &r_right) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar &bend;
    ar &phi0;
  }
};

/** Compute the three-body angle interaction force.
 *  @param[in]  r_mid     Position of second/middle particle.
 *  @param[in]  r_left    Position of first/left particle.
 *  @param[in]  r_right   Position of third/right particle.
 *  @return Forces on the second, first and third particles, in that order.
 */
inline std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
AngleHarmonicBond::forces(Utils::Vector3d const &r_mid,
                          Utils::Vector3d const &r_left,
                          Utils::Vector3d const &r_right) const {

  auto forceFactor = [this](double const cos_phi) {
    auto const sin_phi = sqrt(1 - Utils::sqr(cos_phi));
    auto const phi = acos(cos_phi);
    return -bend * (phi - phi0) / sin_phi;
  };

  return angle_generic_force(r_mid, r_left, r_right, forceFactor, true);
}

/** Compute the three-body angle interaction energy.
 *  @param[in]  r_mid     Position of second/middle particle.
 *  @param[in]  r_left    Position of first/left particle.
 *  @param[in]  r_right   Position of third/right particle.
 */
inline double AngleHarmonicBond::energy(Utils::Vector3d const &r_mid,
                                        Utils::Vector3d const &r_left,
                                        Utils::Vector3d const &r_right) const {
  auto const vectors = calc_vectors_and_cosine(r_mid, r_left, r_right, true);
  auto const cos_phi = std::get<4>(vectors);
  auto const phi = acos(cos_phi);
  return 0.5 * bend * Utils::sqr(phi - phi0);
}

#endif /* ANGLE_HARMONIC_H */

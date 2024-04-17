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

#pragma once

/** \file
 *  Routines to calculate the dihedral energy or/and
 *  force for a particle quadruple. Note that usage of dihedrals
 *  increases the interaction range of bonded interactions to 2 times
 *  the maximal bond length!
 */

#include "config/config.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>

#include <boost/optional.hpp>

#include <cmath>
#include <tuple>

/** Parameters for four-body angular potential (dihedral-angle potentials). */
struct DihedralBond {
  double mult;
  double bend;
  double phase;

  double cutoff() const { return 0.; }

  static constexpr int num = 3;

  DihedralBond(int mult, double bend, double phase) {
    this->mult = mult;
    this->bend = bend;
    this->phase = phase;
  }

  boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d,
                             Utils::Vector3d>>
  forces(Utils::Vector3d const &v12, Utils::Vector3d const &v23,
         Utils::Vector3d const &v34) const;

  boost::optional<double> energy(Utils::Vector3d const &v12,
                                 Utils::Vector3d const &v23,
                                 Utils::Vector3d const &v34) const;

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {
    ar & mult;
    ar & bend;
    ar & phase;
  }
};

/**
 * @brief Calculates the dihedral angle between particle quadruple p1, p2, p3
 * and p4.
 *
 * The dihedral angle is the angle between the planes
 * specified by the particle triples (p1,p2,p3) and (p2,p3,p4).
 * Vectors a, b and c are the bond vectors between consecutive particles.
 * If the a,b or b,c are parallel the dihedral angle is not defined in which
 * case the function returns true. Calling functions should check for that.
 *
 * @param[in]  a Vector from @p p1 to @p p2
 * @param[in]  b Vector from @p p2 to @p p3
 * @param[in]  c Vector from @p p3 to @p p4
 * @param[out] aXb Vector product of a and b
 * @param[out] l_aXb |aXB|
 * @param[out] bXc Vector product of b and c
 * @param[out] l_bXc |bXc|
 * @param[out] cosphi Cosine of the dihedral angle
 * @param[out] phi Dihedral angle in the range [0, pi]
 * @return Whether the angle is undefined.
 */
inline bool calc_dihedral_angle(Utils::Vector3d const &a,
                                Utils::Vector3d const &b,
                                Utils::Vector3d const &c, Utils::Vector3d &aXb,
                                double &l_aXb, Utils::Vector3d &bXc,
                                double &l_bXc, double &cosphi, double &phi) {

  /* calculate vector product a X b and b X c */
  aXb = vector_product(a, b);
  bXc = vector_product(b, c);

  /* calculate the unit vectors */
  l_aXb = aXb.norm();
  l_bXc = bXc.norm();

  /* catch case of undefined dihedral angle */
  if (l_aXb <= TINY_LENGTH_VALUE || l_bXc <= TINY_LENGTH_VALUE) {
    phi = -1.0;
    cosphi = 0.0;
    return true;
  }

  aXb /= l_aXb;
  bXc /= l_bXc;

  cosphi = aXb * bXc;

  if (fabs(fabs(cosphi) - 1) < TINY_SIN_VALUE)
    cosphi = std::round(cosphi);

  /* Calculate dihedral angle */
  phi = acos(cosphi);
  if ((aXb * c) < 0.0)
    phi = (2.0 * Utils::pi()) - phi;
  return false;
}

/** Compute the four-body dihedral interaction force.
 *  The forces have a singularity at @f$ \phi = 0 @f$ and @f$ \phi = \pi @f$
 *  (see @cite swope92a page 592).
 *
 *  @param[in] v12  Vector from @p p1 to @p p2
 *  @param[in] v23  Vector from @p p2 to @p p3
 *  @param[in] v34  Vector from @p p3 to @p p4
 *  @return the forces on @p p2, @p p1, @p p3
 */
inline boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                                  Utils::Vector3d, Utils::Vector3d>>
DihedralBond::forces(Utils::Vector3d const &v12, Utils::Vector3d const &v23,
                     Utils::Vector3d const &v34) const {
  /* vectors for dihedral angle calculation */
  Utils::Vector3d v12Xv23, v23Xv34;
  double l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cos_phi, sin_mphi_over_sin_phi;

  /* dihedral angle */
  auto const angle_is_undefined = calc_dihedral_angle(
      v12, v23, v34, v12Xv23, l_v12Xv23, v23Xv34, l_v23Xv34, cos_phi, phi);
  /* dihedral angle not defined - force zero */
  if (angle_is_undefined) {
    return {};
  }

  auto const f1 = (v23Xv34 - cos_phi * v12Xv23) / l_v12Xv23;
  auto const f4 = (v12Xv23 - cos_phi * v23Xv34) / l_v23Xv34;

  auto const v23Xf1 = vector_product(v23, f1);
  auto const v23Xf4 = vector_product(v23, f4);
  auto const v34Xf4 = vector_product(v34, f4);
  auto const v12Xf1 = vector_product(v12, f1);

  /* calculate force magnitude */
  auto fac = -bend * mult;

  if (fabs(sin(phi)) < TINY_SIN_VALUE) {
    /* comes from taking the first term of the MacLaurin expansion of
     * sin(n * phi - phi0) and sin(phi) and then making the division */
    sin_mphi_over_sin_phi = mult * cos(mult * phi - phase) / cos_phi;
  } else {
    sin_mphi_over_sin_phi = sin(mult * phi - phase) / sin(phi);
  }

  fac *= sin_mphi_over_sin_phi;

  /* store dihedral forces */
  auto const force1 = fac * v23Xf1;
  auto const force2 = fac * (v34Xf4 - v12Xf1 - v23Xf1);
  auto const force3 = fac * (v12Xf1 - v23Xf4 - v34Xf4);

  return std::make_tuple(force2, force1, force3, -(force1 + force2 + force3));
}

/** Compute the four-body dihedral interaction energy.
 *  The energy doesn't have any singularity if the angle phi is well-defined.
 *
 *  @param[in] v12  Vector from @p p1 to @p p2
 *  @param[in] v23  Vector from @p p2 to @p p3
 *  @param[in] v34  Vector from @p p3 to @p p4
 */
inline boost::optional<double>
DihedralBond::energy(Utils::Vector3d const &v12, Utils::Vector3d const &v23,
                     Utils::Vector3d const &v34) const {
  /* vectors for dihedral calculations. */
  Utils::Vector3d v12Xv23, v23Xv34;
  double l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cos_phi;

  auto const angle_is_undefined = calc_dihedral_angle(
      v12, v23, v34, v12Xv23, l_v12Xv23, v23Xv34, l_v23Xv34, cos_phi, phi);
  /* dihedral angle not defined - energy zero */
  if (angle_is_undefined) {
    return {};
  }

  return bend * (1. - cos(mult * phi - phase));
}

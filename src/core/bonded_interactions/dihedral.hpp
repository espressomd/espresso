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

/** \file
 *  Routines to calculate the dihedral energy or/and
 *  force for a particle quadruple. Note that usage of dihedrals
 *  increases the interaction range of bonded interactions to 2 times
 *  the maximal bond length!
 */

#include "config/config.hpp"

#include "errorhandling.hpp"

#include <utils/Vector.hpp>

#include <algorithm>
#include <cmath>
#include <numbers>
#include <optional>
#include <tuple>

/** @brief Tiny length cutoff. */
inline constexpr auto dihe_tiny_length_value{0.0001};

/** Parameters for four-body angular potential (dihedral-angle potentials). */
struct DihedralBond {
  int mult;
  double bend;
  double phase;

  double cutoff() const { return 0.; }

  static constexpr int num = 3;

  DihedralBond(int mult, double bend, double phase) {
    this->mult = mult;
    this->bend = bend;
    this->phase = phase;
  }

  std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d,
                           Utils::Vector3d>>
  forces(Utils::Vector3d const &v12, Utils::Vector3d const &v23,
         Utils::Vector3d const &v34) const;

  std::optional<double> energy(Utils::Vector3d const &v12,
                               Utils::Vector3d const &v23,
                               Utils::Vector3d const &v34) const;
};

/**
 * @brief Calculates the dihedral angle between particle quadruple p1, p2, p3
 * and p4.
 *
 * The dihedral angle is the angle between the planes
 * specified by the particle triples (p1,p2,p3) and (p2,p3,p4).
 * Vectors a, b and c are the bond vectors between consecutive particles.
 * The angle is undefined when any three consecutive particles are collinear,
 * in which case the function returns true. Calling functions should check for
 * that.
 *
 * @param[in]  a Vector from @p p1 to @p p2
 * @param[in]  b Vector from @p p2 to @p p3
 * @param[in]  c Vector from @p p3 to @p p4
 * @param[out] phi Dihedral angle in the range [0, 2 pi)
 * @return Whether the angle is undefined.
 */
inline bool calc_dihedral_angle(Utils::Vector3d const &a,
                                Utils::Vector3d const &b,
                                Utils::Vector3d const &c, double &phi) {

  auto const aXb = vector_product(a, b);
  auto const bXc = vector_product(b, c);
  auto const l_aXb = aXb.norm();
  auto const l_bXc = bXc.norm();

  /* catch case of undefined dihedral angle (three collinear particles) */
  if (l_aXb <= dihe_tiny_length_value or l_bXc <= dihe_tiny_length_value) {
    phi = -1.;
    return true;
  }

  auto const n1 = aXb / l_aXb;
  auto const n2 = bXc / l_bXc;

  /* acos() only accepts values in [-1, 1]. cosphi should mathematically
   * always be in that range, but floating-point round-off can push it a
   * hair outside, so clip it back in to be safe. This only fires for that
   * kind of tiny round-off overshoot -- it leaves alone angles that are
   * just genuinely close to 0 or pi. */
  auto const cosphi = std::clamp(n1 * n2, -1., 1.);
  phi = std::acos(cosphi);
  if ((n1 * c) < 0.)
    phi = 2. * std::numbers::pi - phi;
  return false;
}

/**
 * @brief Gradients of the dihedral angle with respect to the four particle
 * positions.
 *
 * Follows @cite blondel96a eq. 27. The traditional expressions obtained by
 * differentiating @f$ \cos\phi @f$ (@cite swope92a eq. 30) carry a factor
 * @f$ -1/\sin\phi @f$ which is unbounded at @f$ \phi = 0 @f$ and
 * @f$ \phi = \pi @f$. Differentiating @f$ \phi @f$ directly, as done here,
 * avoids that factor entirely: the only denominators are the squared plane
 * normals and the length of the central bond, so the result is well-behaved
 * for every angle. The gradients are undefined only when three consecutive
 * particles are collinear, which is also the only case in which the dihedral
 * angle itself is undefined.
 *
 * @param[in]  a Vector from @p p1 to @p p2
 * @param[in]  b Vector from @p p2 to @p p3
 * @param[in]  c Vector from @p p3 to @p p4
 * @param[out] phi Dihedral angle in the range [0, 2 pi)
 * @param[out] grad1 Gradient of @p phi with respect to the position of @p p1
 * @param[out] grad2 Gradient of @p phi with respect to the position of @p p2
 * @param[out] grad3 Gradient of @p phi with respect to the position of @p p3
 * @param[out] grad4 Gradient of @p phi with respect to the position of @p p4
 * @return Whether the angle is undefined.
 */
inline bool calc_dihedral_angle_gradients(
    Utils::Vector3d const &a, Utils::Vector3d const &b,
    Utils::Vector3d const &c, double &phi, Utils::Vector3d &grad1,
    Utils::Vector3d &grad2, Utils::Vector3d &grad3, Utils::Vector3d &grad4) {

  /* Plane normals. In the notation of @cite blondel96a, F = -a, G = -b and
   * H = c, so that A = F x G = a x b and B = H x G = b x c. */
  auto const A = vector_product(a, b);
  auto const B = vector_product(b, c);
  auto const A_sqr = A.norm2();
  auto const B_sqr = B.norm2();
  auto const l_A = std::sqrt(A_sqr);
  auto const l_B = std::sqrt(B_sqr);
  auto const l_b = b.norm();

  /* catch case of undefined dihedral angle (three collinear particles) */
  if (l_A <= dihe_tiny_length_value or l_B <= dihe_tiny_length_value or
      l_b <= dihe_tiny_length_value) {
    phi = -1.;
    return true;
  }

  auto const n1 = A / l_A;
  auto const n2 = B / l_B;

  /* acos() only accepts values in [-1, 1]; clip a tiny round-off overshoot
   * back in, see @ref calc_dihedral_angle */
  auto const cosphi = std::clamp(n1 * n2, -1., 1.);
  phi = std::acos(cosphi);
  if ((n1 * c) < 0.)
    phi = 2. * std::numbers::pi - phi;

  /* @cite blondel96a eq. 27, with (F.G) = a.b and (H.G) = -(b.c) */
  auto const cA = l_b / A_sqr;              /* |G| / A^2         */
  auto const cB = l_b / B_sqr;              /* |G| / B^2         */
  auto const dA = (a * b) / (A_sqr * l_b);  /* (F.G) / (A^2 |G|) */
  auto const dB = -(b * c) / (B_sqr * l_b); /* (H.G) / (B^2 |G|) */

  grad1 = -cA * A;
  grad2 = (cA + dA) * A - dB * B;
  grad3 = (dB - cB) * B - dA * A;
  grad4 = cB * B;
  return false;
}

/** Compute the four-body dihedral interaction force.
 *  The force is assembled as
 *  @f$ -(\mathrm{d}V/\mathrm{d}\phi)(\partial\phi/\partial r) @f$
 *  (@cite blondel96a eq. 6), which has no singularity at @f$ \phi = 0 @f$ or
 *  @f$ \phi = \pi @f$. See @ref calc_dihedral_angle_gradients.
 *
 *  If three consecutive particles are collinear, the dihedral angle and its
 *  gradients are undefined (@ref calc_dihedral_angle_gradients returns
 *  @c true); a runtime warning is raised and the force is set to zero.
 *
 *  @param[in] v12  Vector from @p p1 to @p p2
 *  @param[in] v23  Vector from @p p2 to @p p3
 *  @param[in] v34  Vector from @p p3 to @p p4
 *  @return the forces on @p p2, @p p1, @p p3, @p p4
 */
inline std::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                                Utils::Vector3d, Utils::Vector3d>>
DihedralBond::forces(Utils::Vector3d const &v12, Utils::Vector3d const &v23,
                     Utils::Vector3d const &v34) const {
  double phi;
  Utils::Vector3d grad1, grad2, grad3, grad4;

  auto const angle_is_undefined = calc_dihedral_angle_gradients(
      v12, v23, v34, phi, grad1, grad2, grad3, grad4);
  if (angle_is_undefined) {
    runtimeWarningMsg() << "Dihedral angle is undefined because three "
                           "consecutive particles are collinear; setting "
                           "the dihedral force to zero";
    return std::make_tuple(Utils::Vector3d{}, Utils::Vector3d{},
                           Utils::Vector3d{}, Utils::Vector3d{});
  }

  auto const mult_ = static_cast<double>(mult);
  auto const dV_dphi = bend * mult_ * std::sin(mult_ * phi - phase);

  return std::make_tuple(-dV_dphi * grad2, -dV_dphi * grad1, -dV_dphi * grad3,
                         -dV_dphi * grad4);
}

/** Compute the four-body dihedral interaction energy.
 *  The energy doesn't have any singularity if the angle phi is well-defined.
 *
 *  If three consecutive particles are collinear, the dihedral angle is
 *  undefined (@ref calc_dihedral_angle returns @c true); a runtime warning is
 *  raised and the energy is set to zero.
 *
 *  @param[in] v12  Vector from @p p1 to @p p2
 *  @param[in] v23  Vector from @p p2 to @p p3
 *  @param[in] v34  Vector from @p p3 to @p p4
 */
inline std::optional<double>
DihedralBond::energy(Utils::Vector3d const &v12, Utils::Vector3d const &v23,
                     Utils::Vector3d const &v34) const {
  double phi;

  auto const angle_is_undefined = calc_dihedral_angle(v12, v23, v34, phi);
  if (angle_is_undefined) {
    runtimeWarningMsg() << "Dihedral angle is undefined because three "
                           "consecutive particles are collinear; setting "
                           "the dihedral energy to zero";
    return 0.;
  }

  auto const mphi = static_cast<double>(mult) * phi - phase;
  return bend * (1. - std::cos(mphi));
}

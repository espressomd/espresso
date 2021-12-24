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

#include "immersed_boundary/ibm_triel.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "grid.hpp"
#include "ibm_common.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/optional.hpp>

#include <tuple>

namespace {
/** Rotate calculated trielastic forces in the 2d plane back to the 3d plane.
 *  Use knowledge that the x-axis in rotated system is parallel to r(p1->p2) in
 *  original system; To find the corresponding unit vector to y in the rotated
 *  system, construct vector perpendicular to r(p1->p2); note that f3 is not
 *  calculated here but is implicitly calculated by f3 = -(f1+f2) which is
 *  consistent with the literature
 *  @return forces on particles 1, 2, 3
 */
std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
RotateForces(Utils::Vector2d const &f1_rot, Utils::Vector2d const &f2_rot,
             Utils::Vector3d const &v12, Utils::Vector3d const &v13) {
  // fRot is in the rotated system, i.e. in a system where the side lPrime of
  // the triangle (i.e. v12) is parallel to the x-axis, and the y-axis is
  // perpendicular to the x-axis (cf. Krueger, Fig. 7.1c).
  // I.e. fRot[XX] is the component parallel to the x-axis, fRot[YY]] the
  // component parallel to the y-axis. Now, the idea is to get these unit
  // vectors for the x- and y-axis in the real coordinate system. They are named
  // xu and yu below. The x-component of the force in the real coordinate system
  // is therefore: fRot[XX]*xu

  // xu is simple: The x-axis in the rotated system is parallel to v12 --> xu =
  // v12 (+ normalization)
  auto const xu = Utils::Vector3d(v12).normalize();

  // yu needs to be orthogonal to xu, and point in the direction of node 3 in
  // Krueger, Fig. 7.1b. Therefore: First get the projection of v13 onto v12:
  // The direction is defined by xu, the length by the scalar product (scalar
  // product can be interpreted as a projection, after all). --> sca * xu Then:
  // v13 - sca * xu gives the component of v13 orthogonal to v12, i.e.
  // perpendicular to the x-axis --> yu Last: Normalize yu.
  auto const yu = (v13 - (v13 * xu) * xu).normalize();

  // Calculate forces in 3D
  auto const force1 = f1_rot[0] * xu + f1_rot[1] * yu;
  auto const force2 = f2_rot[0] * xu + f2_rot[1] * yu;
  auto const force3 = -(force1 + force2);
  return std::make_tuple(force1, force2, force3);
}
} // namespace

boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>>
IBMTriel::calc_forces(Particle const &p1, Particle const &p2,
                      Particle const &p3) const {

  // Calculate the current shape of the triangle (l,lp,cos(phi),sin(phi));
  // l = length between 1 and 3
  // get_mi_vector is an ESPResSo function which considers PBC
  auto const vec2 = box_geo.get_mi_vector(p3.r.p, p1.r.p);
  auto const l = vec2.norm();

  // lp = length between 1 and 2
  auto const vec1 = box_geo.get_mi_vector(p2.r.p, p1.r.p);
  auto const lp = vec1.norm();

  // Check for sanity
  if ((lp - lp0 > maxDist) || (l - l0 > maxDist)) {
    return {};
  }

  // angles between these vectors; calculated directly via the products
  auto const cosPhi = (vec1 * vec2) / (lp * l);
  auto const vecpro = vector_product(vec1, vec2);
  auto const sinPhi = vecpro.norm() / (l * lp);

  // Displacement gradient tensor D: eq. (C.9) in @cite kruger12a
  const double Dxx = lp / lp0;
  const double Dxy = ((l / l0 * cosPhi) - (lp / lp0 * cosPhi0)) / sinPhi0;
  const double Dyx = 0.0;
  const double Dyy = l / l0 * sinPhi / sinPhi0;

  // Tensor G: eq. (C.12)
  const double Gxx = Utils::sqr(Dxx) + Utils::sqr(Dyx);
  const double Gxy = Dxx * Dxy + Dyx * Dyy;
  const double Gyx = Dxx * Dxy + Dyy * Dyx; // = Gxy because of symmetry
  const double Gyy = Utils::sqr(Dxy) + Utils::sqr(Dyy);

  // Strain invariants: eq. (C.11) and (C.12)
  const double i1 = (Gxx + Gyy) - 2;
  const double i2 = ((Gxx * Gyy) - (Gxy * Gyx)) - 1;

  // Derivatives of energy density E used in chain rule below: eq. (C.14)
  double dEdI1;
  double dEdI2;
  if (elasticLaw == tElasticLaw::NeoHookean) {
    // Neo-Hookean
    dEdI1 = k1 / 6.0;
    dEdI2 = -k1 / (6.0 * (i2 + 1.0) * (i2 + 1.0));
  } else {
    // Skalak
    dEdI1 = k1 * (i1 + 1) / 6.0;
    dEdI2 = -k1 / 6.0 + k2 * i2 / 6.0;
  }

  // Derivatives of Is: eq. (C.15)
  const double dI1dGxx = 1;
  const double dI1dGxy = 0;
  const double dI1dGyx = 0;
  const double dI1dGyy = 1;

  const double dI2dGxx = Gyy;
  const double dI2dGxy = -Gyx; // Note: Krueger has a factor 2 here, because he
                               // uses the symmetry of the G-matrix.
  const double dI2dGyx = -Gxy; // But we don't use it. So, Krueger is missing
                               // the yx term, whereas we have it.
  const double dI2dGyy = Gxx;

  // Derivatives of G: eq. (C.16)
  const double dGxxdV1x = 2 * a1 * Dxx;
  const double dGxxdV1y = 0;
  const double dGxxdV2x = 2 * a2 * Dxx;
  const double dGxxdV2y = 0;

  const double dGxydV1x = a1 * Dxy + b1 * Dxx;
  const double dGxydV1y = a1 * Dyy;
  const double dGxydV2x = a2 * Dxy + b2 * Dxx;
  const double dGxydV2y = a2 * Dyy;

  const double dGyxdV1x = a1 * Dxy + b1 * Dxx;
  const double dGyxdV1y = a1 * Dyy;
  const double dGyxdV2x = a2 * Dxy + b2 * Dxx;
  const double dGyxdV2y = a2 * Dyy;

  const double dGyydV1x = 2 * b1 * Dxy;
  const double dGyydV1y = 2 * b1 * Dyy;
  const double dGyydV2x = 2 * b2 * Dxy;
  const double dGyydV2y = 2 * b2 * Dyy;

  // Calculate forces per area in rotated system: chain rule as in appendix C of
  // Krüger (chain rule applied in eq. (C.13), but for the energy density). Only
  // two nodes are needed, third one is calculated from momentum conservation
  // Note: If you calculate the derivatives in a straightforward manner, you get
  // 8 terms (done here). Krueger exploits the symmetry of the G-matrix, which
  // results in 6 elements, but with an additional factor 2 for the xy-elements
  // (see also above at the definition of dI2dGxy).
  Utils::Vector2d f1_rot{};
  Utils::Vector2d f2_rot{};
  f1_rot[0] = -(dEdI1 * dI1dGxx * dGxxdV1x) - (dEdI1 * dI1dGxy * dGxydV1x) -
              (dEdI1 * dI1dGyx * dGyxdV1x) - (dEdI1 * dI1dGyy * dGyydV1x) -
              (dEdI2 * dI2dGxx * dGxxdV1x) - (dEdI2 * dI2dGxy * dGxydV1x) -
              (dEdI2 * dI2dGyx * dGyxdV1x) - (dEdI2 * dI2dGyy * dGyydV1x);
  f1_rot[1] = -(dEdI1 * dI1dGxx * dGxxdV1y) - (dEdI1 * dI1dGxy * dGxydV1y) -
              (dEdI1 * dI1dGyx * dGyxdV1y) - (dEdI1 * dI1dGyy * dGyydV1y) -
              (dEdI2 * dI2dGxx * dGxxdV1y) - (dEdI2 * dI2dGxy * dGxydV1y) -
              (dEdI2 * dI2dGyx * dGyxdV1y) - (dEdI2 * dI2dGyy * dGyydV1y);
  f2_rot[0] = -(dEdI1 * dI1dGxx * dGxxdV2x) - (dEdI1 * dI1dGxy * dGxydV2x) -
              (dEdI1 * dI1dGyx * dGyxdV2x) - (dEdI1 * dI1dGyy * dGyydV2x) -
              (dEdI2 * dI2dGxx * dGxxdV2x) - (dEdI2 * dI2dGxy * dGxydV2x) -
              (dEdI2 * dI2dGyx * dGyxdV2x) - (dEdI2 * dI2dGyy * dGyydV2x);
  f2_rot[1] = -(dEdI1 * dI1dGxx * dGxxdV2y) - (dEdI1 * dI1dGxy * dGxydV2y) -
              (dEdI1 * dI1dGyx * dGyxdV2y) - (dEdI1 * dI1dGyy * dGyydV2y) -
              (dEdI2 * dI2dGxx * dGxxdV2y) - (dEdI2 * dI2dGxy * dGxydV2y) -
              (dEdI2 * dI2dGyx * dGyxdV2y) - (dEdI2 * dI2dGyy * dGyydV2y);

  // Multiply by undeformed area
  f1_rot *= area0;
  f2_rot *= area0;

  // Rotate forces back into original position of triangle
  auto forces = RotateForces(f1_rot, f2_rot, vec1, vec2);

  return forces;
}

IBMTriel::IBMTriel(const int ind1, const int ind2, const int ind3,
                   const double maxDist, const tElasticLaw elasticLaw,
                   const double k1, const double k2) {

  // collect particles from nodes
  auto const pos1 = get_ibm_particle_position(ind1);
  auto const pos2 = get_ibm_particle_position(ind2);
  auto const pos3 = get_ibm_particle_position(ind3);

  // Calculate equilibrium lengths and angle; Note the sequence of the points!
  // l0 = length between 1 and 3
  auto const temp_l0 = box_geo.get_mi_vector(pos3, pos1);
  l0 = temp_l0.norm();
  // lp0 = length between 1 and 2
  auto const temp_lp0 = box_geo.get_mi_vector(pos2, pos1);
  lp0 = temp_lp0.norm();

  // cospo / sinpo angle functions between these vectors; calculated directly
  // via the products
  cosPhi0 = (temp_l0 * temp_lp0) / (l0 * lp0);
  auto const vecpro = vector_product(temp_l0, temp_lp0);
  sinPhi0 = vecpro.norm() / (l0 * lp0);

  // Use the values determined above for further constants of the stretch-force
  // calculation
  const double area2 = l0 * lp0 * sinPhi0;
  a1 = -(l0 * sinPhi0) / area2;
  a2 = -a1;
  b1 = (l0 * cosPhi0 - lp0) / area2;
  b2 = -(l0 * cosPhi0) / area2;
  area0 = 0.5 * area2;
  this->maxDist = maxDist;
  this->elasticLaw = elasticLaw;
  // Always store two constants, for NeoHookean only k1 is used
  this->k1 = k1;
  this->k2 = k2;
}

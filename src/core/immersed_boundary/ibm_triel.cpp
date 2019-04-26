/*
Copyright (C) 2010-2018 The ESPResSo project

This file is part of ESPResSo.

ESPResSo is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

ESPResSo is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "immersed_boundary/ibm_triel.hpp"

#ifdef IMMERSED_BOUNDARY
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

namespace {
/** Rotate calculated trielastic forces in the 2d plane back to the 3d plane
 *Use knowledge that the x-axis in rotated system is parallel to r(p1->p2) in
 *original system; To find the corresponding unit vector to y in the rotated
 *system, construct vector perpendicular to r(p1->p2); note that f3 is not
 *calculated here but is implicitly calculated by f3 = -(f1+f2) which is
 *consistent with the literature
 */
void RotateForces(const double f1_rot[2], const double f2_rot[2], double f1[3],
                  double f2[3], const Utils::Vector3d &v12,
                  const Utils::Vector3d &v13) {
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
  // The direction is definied by xu, the length by the scalar product (scalar
  // product can be interpreted as a projection, after all). --> sca * xu Then:
  // v13 - sca * xu gives the component of v13 orthogonal to v12, i..e.
  // perpendicular to the x-axis --> yu Last: Normalize yu.
  auto const yu = (v13 - (v13 * xu) * xu).normalize();

  // Calculate forces in 3D
  f1[0] = (f1_rot[0] * xu[0]) + (f1_rot[1] * yu[0]);
  f1[1] = (f1_rot[0] * xu[1]) + (f1_rot[1] * yu[1]);
  f1[2] = (f1_rot[0] * xu[2]) + (f1_rot[1] * yu[2]);

  f2[0] = (f2_rot[0] * xu[0]) + (f2_rot[1] * yu[0]);
  f2[1] = (f2_rot[0] * xu[1]) + (f2_rot[1] * yu[1]);
  f2[2] = (f2_rot[0] * xu[2]) + (f2_rot[1] * yu[2]);
}
} // namespace

/*************
   IBM_Triel_CalcForce
Calculate the repulsion and add it to the particle
 **************/

int IBM_Triel_CalcForce(Particle *p1, Particle *p2, Particle *p3,
                        Bonded_ia_parameters *iaparams) {

  // Calculate the current shape of the triangle (l,lp,cos(phi),sin(phi));
  // l = length between 1 and 3
  // get_mi_vector is an Espresso function which considers PBC
  auto const vec2 = get_mi_vector(p3->r.p, p1->r.p);
  auto const l = vec2.norm();

  // lp = lenght between 1 and 2
  auto const vec1 = get_mi_vector(p2->r.p, p1->r.p);
  auto const lp = vec1.norm();

  // angles between these vectors; calculated directly via the products
  double const cosPhi = (vec1 * vec2) / (lp * l);
  auto const vecpro = vector_product(vec1, vec2);
  double const sinPhi = vecpro.norm() / (l * lp);

  // Check for sanity
  if ((lp - iaparams->p.ibm_triel.lp0 > iaparams->p.ibm_triel.maxDist) ||
      (l - iaparams->p.ibm_triel.l0 > iaparams->p.ibm_triel.maxDist)) {
    return 1;
  }

  // Variables in the reference state
  double const l0 = iaparams->p.ibm_triel.l0;
  double const lp0 = iaparams->p.ibm_triel.lp0;
  double const cosPhi0 = iaparams->p.ibm_triel.cosPhi0;
  double const sinPhi0 = iaparams->p.ibm_triel.sinPhi0;
  double const a1 = iaparams->p.ibm_triel.a1;
  double const a2 = iaparams->p.ibm_triel.a2;
  double const b1 = iaparams->p.ibm_triel.b1;
  double const b2 = iaparams->p.ibm_triel.b2;
  double const A0 = iaparams->p.ibm_triel.area0;

  // Displacement gradient tensor D: Eq. (C.9) Krüger thesis
  double const Dxx = lp / lp0;
  double const Dxy = ((l / l0 * cosPhi) - (lp / lp0 * cosPhi0)) / sinPhi0;
  double const Dyx = 0.0;
  double const Dyy = l / l0 * sinPhi / sinPhi0;

  // Tensor G: (C.12)
  double const Gxx = Utils::sqr(Dxx) + Utils::sqr(Dyx);
  double const Gxy = Dxx * Dxy + Dyx * Dyy;
  double const Gyx = Dxx * Dxy + Dyy * Dyx; // = Gxy because of symmetry
  double const Gyy = Utils::sqr(Dxy) + Utils::sqr(Dyy);

  // Strain invariants, C.11 and C.12
  double const i1 = (Gxx + Gyy) - 2;
  double const i2 = ((Gxx * Gyy) - (Gxy * Gyx)) - 1;

  // Derivatives of energy density E used in chain rule below: Eq. (C.14)
  double dEdI1;
  double dEdI2;
  if (iaparams->p.ibm_triel.elasticLaw == tElasticLaw::NeoHookean) {
    // Neo-Hookean
    dEdI1 = iaparams->p.ibm_triel.k1 / 6.0;
    dEdI2 = (-1) * iaparams->p.ibm_triel.k1 / (6.0 * (i2 + 1.0) * (i2 + 1.0));
  } else {
    // Skalak
    dEdI1 = iaparams->p.ibm_triel.k1 * (i1 + 1) / 6.0;
    dEdI2 = (-1) * iaparams->p.ibm_triel.k1 / 6.0 +
            iaparams->p.ibm_triel.k2 * i2 / 6.0;
  }

  // ******** Achim's version *****************

  // Derivatives of Is (C.15)
  double const dI1dGxx = 1;
  double const dI1dGxy = 0;
  double const dI1dGyx = 0;
  double const dI1dGyy = 1;

  double const dI2dGxx = Gyy;
  double const dI2dGxy = -Gyx; // Note: Krueger has a factor 2 here, because he
                               // uses the symmetry of the G-matrix.
  double const dI2dGyx = -Gxy; // But we don't use it. So, Krueger is missing
                               // the yx term, whereas we have it.
  double const dI2dGyy = Gxx;

  // Derivatives of G (C.16)
  double const dGxxdV1x = 2 * a1 * Dxx;
  double const dGxxdV1y = 0;
  double const dGxxdV2x = 2 * a2 * Dxx;
  double const dGxxdV2y = 0;

  double const dGxydV1x = a1 * Dxy + b1 * Dxx;
  double const dGxydV1y = a1 * Dyy;
  double const dGxydV2x = a2 * Dxy + b2 * Dxx;
  double const dGxydV2y = a2 * Dyy;

  double const dGyxdV1x = a1 * Dxy + b1 * Dxx;
  double const dGyxdV1y = a1 * Dyy;
  double const dGyxdV2x = a2 * Dxy + b2 * Dxx;
  double const dGyxdV2y = a2 * Dyy;

  double const dGyydV1x = 2 * b1 * Dxy;
  double const dGyydV1y = 2 * b1 * Dyy;
  double const dGyydV2x = 2 * b2 * Dxy;
  double const dGyydV2y = 2 * b2 * Dyy;

  // Calculate forces per area in rotated system: chain rule as in appendix C of
  // Krüger (chain rule applied in eq. (C.13), but for the energy density). Only
  // two nodes are needed, third one is calculated from momentum conservation
  // Note: If you calculate the derivatives in a straightforward manner, you get
  // 8 terms (done here). Krueger exploits the symmetry of the G-matrix, which
  // results in 6 elements, but with an additional factor 2 for the xy-elements
  // (see also above at the definition of dI2dGxy).
  double f1_rot[2] = {0., 0.};
  double f2_rot[2] = {0., 0.};
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
  f1_rot[0] *= A0;
  f1_rot[1] *= A0;
  f2_rot[0] *= A0;
  f2_rot[1] *= A0;

  // ****************** Wolfgang's version ***********
  /*
   // Left here for checking, but should be identical to the version above
   double const i11 = 1.0;
   double const i12 = 1.0;
   double const i21 = Gyy;
   double const i22 = -Gyx;
   double const i23 = i22;
   double const i24 = Gxx;

   //For sake of better readability shorten the call for the triangle's
   constants: A0 = iaparams->p.stretching_force_ibm.Area0; a1 =
   iaparams->p.stretching_force_ibm.a1; a2 =
   iaparams->p.stretching_force_ibm.a2; b1 =
   iaparams->p.stretching_force_ibm.b1; b2 =
   iaparams->p.stretching_force_ibm.b2;

   f1_rot[0] = A0*((-1)*e1*((i11*2*a1*dxx)+(i12*2*b1*dxy))+
   (-1)*e2*((i21*2*a1*dxx)+(i22*(a1*dxy+b1*dxx))+(i23*(a1*dxy+b1*dxx))+(i24*2*b1*dxy)));
   f1_rot[1] = A0*((-1)*e1*((i11*0.0)+(i12*2*b1*dyy))+
   (-1)*e2*((i21*0.0)+(i22*a1*dyy)+(i23*a1*dyy)+(i24*2*b1*dyy)));

   f2_rot[0] = A0*((-1)*e1*((i11*2*a2*dxx)+(i12*2*b2*dxy))+
   (-1)*e2*((i21*2*a2*dxx)+(i22*(a2*dxy+b2*dxx))+(i23*(a2*dxy+b2*dxx))+(i24*2*b2*dxy)));
   f2_rot[1] = A0*((-1)*e1*((i11*0.0)+(i12*2*b2*dyy))+
   (-1)*e2*((i21*0.0)+(i22*a2*dyy)+(i23*a2*dyy)+(i24*2*b2*dyy)));
   */

  // Rotate forces back into original position of triangle
  double force1[3] = {0, 0, 0};
  double force2[3] = {0, 0, 0};
  RotateForces(f1_rot, f2_rot, force1, force2, vec1, vec2);

  // Calculate f3 from equilibrium and add
  for (int i = 0; i < 3; i++) {
    p1->f.f[i] += force1[i];
    p2->f.f[i] += force2[i];
    p3->f.f[i] += -force1[i] - force2[i];
  }

  return 0;
}

/****************
  IBM_Triel_ResetParams
 *****************/

int IBM_Triel_ResetParams(const int bond_type, const double k1,
                          const double l0) {

  // Check if bond exists and is of correct type
  if (bond_type >= bonded_ia_params.size()) {
    printf("bond does not exist while reading triel checkpoint\n");
    return ES_ERROR;
  }
  if (bonded_ia_params[bond_type].type != BONDED_IA_IBM_TRIEL) {
    printf("interaction type does not match while reading triel checkpoint!\n");
    return ES_ERROR;
  }

  // Check if k1 is correct
  if (fabs(bonded_ia_params[bond_type].p.ibm_triel.k1 - k1) > 1e-9) {
    printf("k1 does not match while reading triel checkpoint!\n");
    return ES_ERROR;
  }

  // Check if l0 is correct
  if (fabs(bonded_ia_params[bond_type].p.ibm_triel.l0 - l0) > 1e-9) {
    printf("l0 does not match while reading triel checkpoint!\n");
    return ES_ERROR;
  }

  // Compute cache values a1, a2, b1, b2
  double const area0 = bonded_ia_params[bond_type].p.ibm_triel.area0;
  double const lp0 = bonded_ia_params[bond_type].p.ibm_triel.lp0;
  double const sinPhi0 = bonded_ia_params[bond_type].p.ibm_triel.sinPhi0;
  double const cosPhi0 = bonded_ia_params[bond_type].p.ibm_triel.cosPhi0;
  double const area2 = 2.0 * area0;
  double const a1 = -(l0 * sinPhi0) / area2;
  double const a2 = -a1;
  double const b1 = (l0 * cosPhi0 - lp0) / area2;
  double const b2 = -(l0 * cosPhi0) / area2;

  // Hand these values over to parameter structure
  bonded_ia_params[bond_type].p.ibm_triel.a1 = a1;
  bonded_ia_params[bond_type].p.ibm_triel.a2 = a2;
  bonded_ia_params[bond_type].p.ibm_triel.b1 = b1;
  bonded_ia_params[bond_type].p.ibm_triel.b2 = b2;

  // Communicate this to whoever is interested
  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

/***********
   IBM_Triel_SetParams
************/

int IBM_Triel_SetParams(const int bond_type, const int ind1, const int ind2,
                        const int ind3, const double maxDist,
                        const tElasticLaw elasticLaw, const double k1,
                        const double k2) {
  // Create bond
  make_bond_type_exist(bond_type);

  // Get data (especially location) of three particles
  auto part1 = get_particle_data(ind1);
  auto part2 = get_particle_data(ind2);
  auto part3 = get_particle_data(ind3);

  // Calculate equilibrium lengths and angle; Note the sequence of the points!
  // lo = length between 1 and 3
  auto const templo = get_mi_vector(part3.r.p, part1.r.p);
  double const l0 = templo.norm();
  // lpo = length between 1 and 2
  auto const templpo = get_mi_vector(part2.r.p, part1.r.p);
  double const lp0 = templpo.norm();

  // cospo / sinpo angle functions between these vectors; calculated directly
  // via the products
  double const cosPhi0 = (templo * templpo) / (l0 * lp0);
  auto const vecpro = vector_product(templo, templpo);
  double const sinPhi0 = vecpro.norm() / (l0 * lp0);

  // Use the values determined above for further constants of the stretch-force
  // calculation
  double const area2 = l0 * lp0 * sinPhi0;
  double const a1 = -(l0 * sinPhi0) / area2;
  double const a2 = -a1;
  double const b1 = (l0 * cosPhi0 - lp0) / area2;
  double const b2 = -(l0 * cosPhi0) / area2;

  // General stuff
  bonded_ia_params[bond_type].type = BONDED_IA_IBM_TRIEL;
  bonded_ia_params[bond_type].num = 2;

  // Specific stuff
  bonded_ia_params[bond_type].p.ibm_triel.a1 = a1;
  bonded_ia_params[bond_type].p.ibm_triel.a2 = a2;
  bonded_ia_params[bond_type].p.ibm_triel.b1 = b1;
  bonded_ia_params[bond_type].p.ibm_triel.b2 = b2;
  bonded_ia_params[bond_type].p.ibm_triel.l0 = l0;
  bonded_ia_params[bond_type].p.ibm_triel.lp0 = lp0;
  bonded_ia_params[bond_type].p.ibm_triel.sinPhi0 = sinPhi0;
  bonded_ia_params[bond_type].p.ibm_triel.cosPhi0 = cosPhi0;
  bonded_ia_params[bond_type].p.ibm_triel.area0 = 0.5 * area2;
  bonded_ia_params[bond_type].p.ibm_triel.maxDist = maxDist;
  bonded_ia_params[bond_type].p.ibm_triel.elasticLaw = elasticLaw;
  // Always store two constants, for NeoHookean only k1 is used
  bonded_ia_params[bond_type].p.ibm_triel.k1 = k1;
  bonded_ia_params[bond_type].p.ibm_triel.k2 = k2;

  // Communicate this to whoever is interested
  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}
#endif

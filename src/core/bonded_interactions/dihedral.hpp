/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
#ifndef DIHEDRAL_H
#define DIHEDRAL_H
/** \file
 *  Routines to calculate the dihedral energy or/and
 *  force for a particle quadruple.  Note that usage of dihedrals
 *  increases the interaction range of bonded interactions to 2 times
 *  the maximal bond length!  \ref forces.cpp
 */

#include "bonded_interaction_data.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>

/** set dihedral parameters
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int dihedral_set_params(int bond_type, int mult, double bend, double phase);

/**
 * @brief Calculates the dihedral angle between particle quadruple p1, p2, p3
 * and p4.
 *
 * The dihedral angle is the angle between the planes
 * specified by the particle triples (p1,p2,p3) and (p2,p3,p4).
 * Vectors a, b and c are the bond vectors between consecutive particles.
 * If the a,b or b,c are parallel the dihedral angle is not defined in which
 * case the routine returns phi=-1. Calling functions should check for that
 *
 * @param[in]  p1 , p2 , p3 , p4 Particles forming the dihedral
 * @param[out] a Vector from @p p1 to @p p2
 * @param[out] b Vector from @p p2 to @p p3
 * @param[out] c Vector from @p p3 to @p p4
 * @param[out] aXb Vector product of a and b
 * @param[out] l_aXb |aXB|
 * @param[out] bXc Vector product of b and c
 * @param[out] l_bXc |bXc|
 * @param[out] cosphi Cosine of the dihedral angle
 * @param[out] phi Dihedral angle
 */
inline void calc_dihedral_angle(Particle const *p1, Particle const *p2,
                                Particle const *p3, Particle const *p4,
                                Utils::Vector3d &a, Utils::Vector3d &b,
                                Utils::Vector3d &c, Utils::Vector3d &aXb,
                                double *l_aXb, Utils::Vector3d &bXc,
                                double *l_bXc, double *cosphi, double *phi) {
  a = get_mi_vector(p2->r.p, p1->r.p);
  b = get_mi_vector(p3->r.p, p2->r.p);
  c = get_mi_vector(p4->r.p, p3->r.p);

  /* calculate vector product a X b and b X c */
  aXb = vector_product(a, b);
  bXc = vector_product(b, c);

  /* calculate the unit vectors */
  *l_aXb = aXb.norm();
  *l_bXc = bXc.norm();

  /* catch case of undefined dihedral angle */
  if (*l_aXb <= TINY_LENGTH_VALUE || *l_bXc <= TINY_LENGTH_VALUE) {
    *phi = -1.0;
    *cosphi = 0;
    return;
  }

  aXb /= *l_aXb;
  bXc /= *l_bXc;

  *cosphi = aXb * bXc;

  if (fabs(fabs(*cosphi) - 1) < TINY_SIN_VALUE)
    *cosphi = std::round(*cosphi);

  /* Calculate dihedral angle */
  *phi = acos(*cosphi);
  if ((aXb * c) < 0.0)
    *phi = (2.0 * Utils::pi()) - *phi;
}

/** calculate dihedral force between particles p1, p2 p3 and p4
    Written by Arijit Maitra, adapted to new force interface by Hanjo,
    more general new dihedral form by Ana.
*/
inline int calc_dihedral_force(Particle const *p2, Particle const *p1,
                               Particle const *p3, Particle const *p4,
                               Bonded_ia_parameters const *iaparams,
                               double force2[3], double force1[3],
                               double force3[3]) {
  /* vectors for dihedral angle calculation */
  Utils::Vector3d v12, v23, v34, v12Xv23, v23Xv34;
  double l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi, sinmphi_sinphi;
  /* force factors */
  double fac;

  /* dihedral angle */
  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23,
                      v23Xv34, &l_v23Xv34, &cosphi, &phi);
  /* dihedral angle not defined - force zero */
  if (phi == -1.0) {
    for (int i = 0; i < 3; i++) {
      force1[i] = 0.0;
      force2[i] = 0.0;
      force3[i] = 0.0;
    }
    return 0;
  }

  auto const f1 = (v23Xv34 - cosphi * v12Xv23) / l_v12Xv23;
  auto const f4 = (v12Xv23 - cosphi * v23Xv34) / l_v23Xv34;

  auto const v23Xf1 = vector_product(v23, f1);
  auto const v23Xf4 = vector_product(v23, f4);
  auto const v34Xf4 = vector_product(v34, f4);
  auto const v12Xf1 = vector_product(v12, f1);

  /* calculate force magnitude */
  fac = -iaparams->p.dihedral.bend * iaparams->p.dihedral.mult;

  if (fabs(sin(phi)) < TINY_SIN_VALUE) {
    /*(comes from taking the first term of the MacLaurin expansion of
      sin(n*phi - phi0) and sin(phi) and then making the division).
      The original code had a 2PI term in the cosine (cos(2PI - nPhi))
      but I removed it because it wasn't doing anything. AnaVV*/
    sinmphi_sinphi =
        iaparams->p.dihedral.mult *
        cos(iaparams->p.dihedral.mult * phi - iaparams->p.dihedral.phase) /
        cosphi;
  } else {
    sinmphi_sinphi =
        sin(iaparams->p.dihedral.mult * phi - iaparams->p.dihedral.phase) /
        sin(phi);
  }

  fac *= sinmphi_sinphi;

  /* store dihedral forces */
  for (int i = 0; i < 3; i++) {
    force1[i] = fac * v23Xf1[i];
    force2[i] = fac * (v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
    force3[i] = fac * (v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }
  return 0;
}

/** calculate dihedral energy between particles p1, p2 p3 and p4
    Written by Arijit Maitra, adapted to new force interface by Hanjo */
inline int dihedral_energy(Particle const *p1, Particle const *p2,
                           Particle const *p3, Particle const *p4,
                           Bonded_ia_parameters const *iaparams,
                           double *_energy) {
  /* vectors for dihedral calculations. */
  Utils::Vector3d v12, v23, v34, v12Xv23, v23Xv34;
  double l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23,
                      v23Xv34, &l_v23Xv34, &cosphi, &phi);

  *_energy =
      iaparams->p.dihedral.bend *
      (1. - cos(iaparams->p.dihedral.mult * phi - iaparams->p.dihedral.phase));

  return 0;
}

#endif

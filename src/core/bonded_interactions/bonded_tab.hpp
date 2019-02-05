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
#ifndef CORE_BONDED_INTERACTIONS_TABULATED_HPP
#define CORE_BONDED_INTERACTIONS_TABULATED_HPP

/** \file
 *  Routines to calculate the  energy and/or  force
 *  for a particle pair or bonds via interpolating from lookup tables.
 *  \ref forces.cpp
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 */

#include "config.hpp"

#ifdef TABULATED
#include "angle_common.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/dihedral.hpp"
#include "debug.hpp"
#include "particle_data.hpp"
#include "utils.hpp"

/** Non-Bonded tabulated potentials:
 *  Reads tabulated parameters and force and energy tables from a file.
 *  ia_params and force/energy tables are then communicated to each node
 *
 *  @param part_type_a  %Particle type for which the interaction is defined
 *  @param part_type_b  %Particle type for which the interaction is defined
 *  @param min          @copybrief TabulatedPotential::minval
 *  @param max          @copybrief TabulatedPotential::maxval
 *  @param energy       @copybrief TabulatedPotential::energy_tab
 *  @param force        @copybrief TabulatedPotential::force_tab
 *
 *  @retval 0 on success
 *  @retval 1 on particle type mismatches
 *  @retval 2 file name too long
 *  @retval 3 cannot open the file
 *  @retval 4 file too short
 *  @retval 5 file broken, cannot parse numbers
 *  @retval 6 number of points of existing potential changed
 */
int tabulated_set_params(int part_type_a, int part_type_b, double min,
                         double max, std::vector<double> const &energy,
                         std::vector<double> const &force);

/** Bonded tabulated potentials: Reads tabulated parameters and force
 *  and energy tables from a file.  ia_params and force/energy tables
 *  are then communicated to each node.
 *
 *  @param bond_type    Bond type for which the interaction is defined
 *  @param tab_type     Table type
 *  @param min          @copybrief TabulatedPotential::minval
 *  @param max          @copybrief TabulatedPotential::maxval
 *  @param energy       @copybrief TabulatedPotential::energy_tab
 *  @param force        @copybrief TabulatedPotential::force_tab
 *
 *  @retval 0 on success
 *  @retval 1 if wrong bond type
 *  @retval 2 currently unused
 *  @retval 3 cannot open the file
 *  @retval 4 file too short
 *  @retval 5 file broken, cannot parse numbers
 *  @retval 6 parameter out of bounds
 */
int tabulated_bonded_set_params(int bond_type,
                                TabulatedBondedInteraction tab_type, double min,
                                double max, std::vector<double> const &energy,
                                std::vector<double> const &force);

/* BONDED INTERACTIONS */

/** Compute a tabulated bond length force.
 *
 *  The force acts in the direction of the connecting vector between the
 *  particles. For distances smaller than the tabulated range it uses a linear
 *  extrapolation based on the first two tabulated force values.
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 *
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] force     Force.
 *  @retval 1 if the bond is broken
 *  @retval 0 otherwise
 */
inline int calc_tab_bond_force(Particle const *p1, Particle const *p2,
                               Bonded_ia_parameters const *iaparams,
                               double const dx[3], double force[3]) {
  auto const *tab_pot = iaparams->p.tab.pot;
  auto const dist = sqrt(sqrlen(dx));

  if (dist < tab_pot->cutoff()) {
    auto const fac = tab_pot->force(dist) / dist;

    for (int j = 0; j < 3; j++)
      force[j] -= fac * dx[j];

    return 0;
  } else {
    return 1;
  }
}

/** Compute a tabulated bond length energy.
 *
 *  For distances smaller than the tabulated range it uses a quadratic
 *  extrapolation based on the first two tabulated force values and the first
 *  tabulated energy value.
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 *
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @param[out] _energy   Energy.
 *  @retval 1 if the bond is broken
 *  @retval 0 otherwise
 */
inline int tab_bond_energy(Particle const *p1, Particle const *p2,
                           Bonded_ia_parameters const *iaparams,
                           double const dx[3], double *_energy) {
  auto const *tab_pot = iaparams->p.tab.pot;
  double dist = sqrt(sqrlen(dx));

  if (dist < tab_pot->cutoff()) {
    *_energy = tab_pot->energy(dist);
    return 0;
  } else {
    return 1;
  }
}

/** Compute the three-body angle interaction force.
 *
 *  The force on @p p_left and @p p_right acts perpendicular to the connecting
 *  vector between the particle and @p p_mid and in the plane defined by the
 *  three particles. The force on the middle particle balances the other two
 *  forces. The forces are scaled with the inverse length of the
 *  connecting vectors. It is assumed that the potential is tabulated
 *  for all angles between 0 and Pi.
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 *
 *  @param[in]  p_mid     Second/middle particle.
 *  @param[in]  p_left    First/left particle.
 *  @param[in]  p_right   Third/right particle.
 *  @param[in]  iaparams  Bonded parameters for the angle interaction.
 *  @param[out] force1    Force on particle 1.
 *  @param[out] force2    Force on particle 2.
 *  @return 0
 */
inline int calc_tab_angle_force(Particle const *p_mid, Particle const *p_left,
                                Particle const *p_right,
                                Bonded_ia_parameters const *iaparams,
                                double force1[3], double force2[3]) {

  auto forceFactor = [&iaparams](double &cosine) {
    double phi;
#ifdef TABANGLEMINUS
    phi = acos(-cosine);
#else
    phi = acos(cosine);
#endif
    double invsinphi = sin(phi);
    if (invsinphi < TINY_SIN_VALUE)
      invsinphi = TINY_SIN_VALUE;
    invsinphi = 1.0 / invsinphi;
    /* look up force factor */
    auto const *tab_pot = iaparams->p.tab.pot;
    auto fac = tab_pot->force(phi) * invsinphi;
    return fac;
  };

  calc_angle_generic_force(p_mid, p_left, p_right, forceFactor, force1, force2);

  return 0;
}

/* The force on each particle due to a three-body bonded tabulated
   potential is computed. */
inline void calc_angle_3body_tabulated_forces(
    Particle const *p_mid, Particle const *p_left, Particle const *p_right,
    Bonded_ia_parameters const *iaparams, Vector3d &force1, Vector3d &force2,
    Vector3d &force3) {

  auto forceFactor = [&iaparams](double &cos_phi, double &sin_phi) {
    if (cos_phi < -1.0)
      cos_phi = -TINY_COS_VALUE;
    if (cos_phi > 1.0)
      cos_phi = TINY_COS_VALUE;
#ifdef TABANGLEMINUS
    auto phi = acos(-cos_phi);
#else
    auto phi = acos(cos_phi);
#endif
    auto const *tab_pot = iaparams->p.tab.pot;
    auto dU = tab_pot->force(phi); // d/dphi of U(phi)
    // potential dependent term (dU/dphi * 1 / sin(phi))
    auto pot_dep = dU / sin_phi;
    return pot_dep;
  };

  calc_angle_generic_3body_forces(p_mid, p_left, p_right, forceFactor, force1,
                                  force2, force3);
}

/** Compute the three-body angle interaction energy.
 *  It is assumed that the potential is tabulated
 *  for all angles between 0 and Pi.
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 *
 *  @param[in]  p_mid     Second/middle particle.
 *  @param[in]  p_left    First/left particle.
 *  @param[in]  p_right   Third/right particle.
 *  @param[in]  iaparams  Bonded parameters for the angle interaction.
 *  @param[out] _energy   Energy.
 *  @return 0
 */
inline int tab_angle_energy(Particle const *p_mid, Particle const *p_left,
                            Particle const *p_right,
                            Bonded_ia_parameters const *iaparams,
                            double *_energy) {
  auto const *tab_pot = iaparams->p.tab.pot;
  /* vector from p_left to p_mid */
  auto vec1 = get_mi_vector(p_mid->r.p, p_left->r.p);
  double d1i = 1.0 / vec1.norm();
  vec1 *= d1i;
  /* vector from p_mid to p_right */
  auto vec2 = get_mi_vector(p_right->r.p, p_mid->r.p);
  double d2i = 1.0 / vec2.norm();
  vec2 *= d2i;
  /* calculate phi */
#ifdef TABANGLEMINUS
  auto phi = acos(-scalar(vec1, vec2));
#else
  auto phi = acos(scalar(vec1, vec2));
#endif

  *_energy = tab_pot->energy(phi);

  return 0;
}

/** Compute the four-body dihedral interaction force.
 *  This function is not tested yet.
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 *
 *  @param[in]  p2        Second particle.
 *  @param[in]  p1        First particle.
 *  @param[in]  p3        Third particle.
 *  @param[in]  p4        Fourth particle.
 *  @param[in]  iaparams  Bonded parameters for the dihedral interaction.
 *  @param[out] force2    Force on particle 2.
 *  @param[out] force1    Force on particle 1.
 *  @param[out] force3    Force on particle 3.
 *  @return 0
 */
inline int calc_tab_dihedral_force(Particle const *p2, Particle const *p1,
                                   Particle const *p3, Particle const *p4,
                                   Bonded_ia_parameters const *iaparams,
                                   double force2[3], double force1[3],
                                   double force3[3]) {
  int i;
  /* vectors for dihedral angle calculation */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  double v23Xf1[3], v23Xf4[3], v34Xf4[3], v12Xf1[3];
  /* dihedral angle, cosine of the dihedral angle, cosine of the bond angles */
  double phi, cosphi;
  /* force factors */
  double fac, f1[3], f4[3];
  auto const *tab_pot = iaparams->p.tab.pot;

  /* dihedral angle */
  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23,
                      v23Xv34, &l_v23Xv34, &cosphi, &phi);
  /* dihedral angle not defined - force zero */
  if (phi == -1.0) {
    for (i = 0; i < 3; i++) {
      force1[i] = 0.0;
      force2[i] = 0.0;
      force3[i] = 0.0;
    }
    return 0;
  }

  /* calculate force components (directions) */
  for (i = 0; i < 3; i++) {
    f1[i] = (v23Xv34[i] - cosphi * v12Xv23[i]) / l_v12Xv23;
    ;
    f4[i] = (v12Xv23[i] - cosphi * v23Xv34[i]) / l_v23Xv34;
  }
  vector_product(v23, f1, v23Xf1);
  vector_product(v23, f4, v23Xf4);
  vector_product(v34, f4, v34Xf4);
  vector_product(v12, f1, v12Xf1);

  /* table lookup */
  fac = tab_pot->force(phi);

  /* store dihedral forces */
  for (i = 0; i < 3; i++) {
    force1[i] = fac * v23Xf1[i];
    force2[i] = fac * (v34Xf4[i] - v12Xf1[i] - v23Xf1[i]);
    force3[i] = fac * (v12Xf1[i] - v23Xf4[i] - v34Xf4[i]);
  }

  return 0;
}

/** Compute the four-body dihedral interaction energy.
 *  This function is not tested yet.
 *  Needs feature TABULATED compiled in (see \ref config.hpp).
 *
 *  @param[in]  p2        Second particle.
 *  @param[in]  p1        First particle.
 *  @param[in]  p3        Third particle.
 *  @param[in]  p4        Fourth particle.
 *  @param[in]  iaparams  Bonded parameters for the dihedral interaction.
 *  @param[out] _energy   Energy.
 *  @return 0
 */
inline int tab_dihedral_energy(Particle const *p2, Particle const *p1,
                               Particle const *p3, Particle const *p4,
                               Bonded_ia_parameters const *iaparams,
                               double *_energy) {
  /* vectors for dihedral calculations. */
  double v12[3], v23[3], v34[3], v12Xv23[3], v23Xv34[3], l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cosphi;
  auto const *tab_pot = iaparams->p.tab.pot;

  calc_dihedral_angle(p1, p2, p3, p4, v12, v23, v34, v12Xv23, &l_v12Xv23,
                      v23Xv34, &l_v23Xv34, &cosphi, &phi);

  *_energy = tab_pot->energy(phi);

  return 0;
}

#endif
#endif

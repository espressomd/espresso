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
#ifndef CORE_BONDED_INTERACTIONS_TABULATED_HPP
#define CORE_BONDED_INTERACTIONS_TABULATED_HPP

/** \file
 *  Routines to calculate the energy and/or force for particle bonds, angles
 *  and dihedrals via interpolation of lookup tables.
 *
 *  Implementation in \ref bonded_tab.cpp.
 */

#include "config.hpp"

#include "angle_common.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "bonded_interactions/dihedral.hpp"
#include <tuple>

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

/** Set the parameters of a bonded tabulated potential.
 *  ia_params and force/energy tables are communicated to each node.
 *
 *  @param bond_type    Bond type for which the interaction is defined
 *  @param tab_type     Table type
 *  @param min          @copybrief TabulatedPotential::minval
 *  @param max          @copybrief TabulatedPotential::maxval
 *  @param energy       @copybrief TabulatedPotential::energy_tab
 *  @param force        @copybrief TabulatedPotential::force_tab
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int tabulated_bonded_set_params(int bond_type,
                                TabulatedBondedInteraction tab_type, double min,
                                double max, std::vector<double> const &energy,
                                std::vector<double> const &force);

/** Compute a tabulated bond length force.
 *
 *  The force acts in the direction of the connecting vector between the
 *  particles. For distances smaller than the tabulated range it uses a linear
 *  extrapolation based on the first two tabulated force values.
 *
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<Utils::Vector3d>
tab_bond_force(Bonded_ia_parameters const &iaparams,
               Utils::Vector3d const &dx) {
  auto const *tab_pot = iaparams.p.tab.pot;
  auto const dist = dx.norm();

  if (dist < tab_pot->cutoff()) {
    auto const fac = tab_pot->force(dist) / dist;
    return fac * dx;
  }
  return {};
}

/** Compute a tabulated bond length energy.
 *
 *  For distances smaller than the tabulated range it uses a quadratic
 *  extrapolation based on the first two tabulated force values and the first
 *  tabulated energy value.
 *
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<double>
tab_bond_energy(Bonded_ia_parameters const &iaparams,
                Utils::Vector3d const &dx) {
  auto const *tab_pot = iaparams.p.tab.pot;
  auto const dist = dx.norm();

  if (dist < tab_pot->cutoff()) {
    return tab_pot->energy(dist);
  }
  return {};
}

/** Compute the three-body angle interaction force.
 *  @param[in]  r_mid     Position of second/middle particle.
 *  @param[in]  r_left    Position of first/left particle.
 *  @param[in]  r_right   Position of third/right particle.
 *  @param[in]  iaparams  Bonded parameters for the angle interaction.
 *  @return Forces on the second, first and third particles, in that order.
 */
inline std::tuple<Utils::Vector3d, Utils::Vector3d, Utils::Vector3d>
tab_angle_force(Utils::Vector3d const &r_mid, Utils::Vector3d const &r_left,
                Utils::Vector3d const &r_right,
                Bonded_ia_parameters const &iaparams) {

  auto forceFactor = [&iaparams](double const cos_phi) {
    auto const sin_phi = sqrt(1 - Utils::sqr(cos_phi));
#ifdef TABANGLEMINUS
    auto const phi = acos(-cos_phi);
#else
    auto const phi = acos(cos_phi);
#endif
    auto const *tab_pot = iaparams.p.tab.pot;
    auto const gradient = tab_pot->force(phi);
    return -gradient / sin_phi;
  };

  return angle_generic_force(r_mid, r_left, r_right, forceFactor, true);
}

/** Compute the three-body angle interaction energy.
 *  It is assumed that the potential is tabulated
 *  for all angles between 0 and Pi.
 *
 *  @param[in]  r_mid     Position of second/middle particle.
 *  @param[in]  r_left    Position of first/left particle.
 *  @param[in]  r_right   Position of third/right particle.
 *  @param[in]  iaparams  Bonded parameters for the angle interaction.
 */
inline double tab_angle_energy(Utils::Vector3d const &r_mid,
                               Utils::Vector3d const &r_left,
                               Utils::Vector3d const &r_right,
                               Bonded_ia_parameters const &iaparams) {
  auto const vectors = calc_vectors_and_cosine(r_mid, r_left, r_right, true);
  auto const cos_phi = std::get<4>(vectors);
  /* calculate phi */
#ifdef TABANGLEMINUS
  auto const phi = acos(-cos_phi);
#else
  auto const phi = acos(cos_phi);
#endif
  return iaparams.p.tab.pot->energy(phi);
}

/** Compute the four-body dihedral interaction force.
 *  This function is not tested yet.
 *
 *  @param[in]  r1        Position of the first particle.
 *  @param[in]  r2        Position of the second particle.
 *  @param[in]  r3        Position of the third particle.
 *  @param[in]  r4        Position of the fourth particle.
 *  @param[in]  iaparams  Bonded parameters for the dihedral interaction.
 *  @return the forces on @p p2, @p p1, @p p3
 */
inline boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d,
                                  Utils::Vector3d, Utils::Vector3d>>
tab_dihedral_force(Utils::Vector3d const &r1, Utils::Vector3d const &r2,
                   Utils::Vector3d const &r3, Utils::Vector3d const &r4,
                   Bonded_ia_parameters const &iaparams) {
  /* vectors for dihedral angle calculation */
  Utils::Vector3d v12, v23, v34, v12Xv23, v23Xv34;
  double l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle, cosine of the bond angles */
  double phi, cos_phi;
  /* force factors */
  auto const *tab_pot = iaparams.p.tab.pot;

  /* dihedral angle */
  calc_dihedral_angle(r1, r2, r3, r4, v12, v23, v34, v12Xv23, &l_v12Xv23,
                      v23Xv34, &l_v23Xv34, &cos_phi, &phi);
  /* dihedral angle not defined - force zero */
  if (phi == -1.0) {
    return {};
  }

  /* calculate force components (directions) */
  auto const f1 = (v23Xv34 - cos_phi * v12Xv23) / l_v12Xv23;
  auto const f4 = (v12Xv23 - cos_phi * v23Xv34) / l_v23Xv34;

  auto const v23Xf1 = vector_product(v23, f1);
  auto const v23Xf4 = vector_product(v23, f4);
  auto const v34Xf4 = vector_product(v34, f4);
  auto const v12Xf1 = vector_product(v12, f1);

  /* table lookup */
  auto const fac = tab_pot->force(phi);

  /* store dihedral forces */
  auto const force1 = fac * v23Xf1;
  auto const force2 = fac * (v34Xf4 - v12Xf1 - v23Xf1);
  auto const force3 = fac * (v12Xf1 - v23Xf4 - v34Xf4);

  return std::make_tuple(force2, force1, force3, -(force2 + force1 + force3));
}

/** Compute the four-body dihedral interaction energy.
 *  This function is not tested yet.
 *
 *  @param[in]  r1        Position of the first particle.
 *  @param[in]  r2        Position of the second particle.
 *  @param[in]  r3        Position of the third particle.
 *  @param[in]  r4        Position of the fourth particle.
 *  @param[in]  iaparams  Bonded parameters for the dihedral interaction.
 */
inline boost::optional<double>
tab_dihedral_energy(Utils::Vector3d const &r1, Utils::Vector3d const &r2,
                    Utils::Vector3d const &r3, Utils::Vector3d const &r4,
                    Bonded_ia_parameters const &iaparams) {
  /* vectors for dihedral calculations. */
  Utils::Vector3d v12, v23, v34, v12Xv23, v23Xv34;
  double l_v12Xv23, l_v23Xv34;
  /* dihedral angle, cosine of the dihedral angle */
  double phi, cos_phi;
  auto const *tab_pot = iaparams.p.tab.pot;
  calc_dihedral_angle(r1, r2, r3, r4, v12, v23, v34, v12Xv23, &l_v12Xv23,
                      v23Xv34, &l_v23Xv34, &cos_phi, &phi);
  return tab_pot->energy(phi);
}

#endif

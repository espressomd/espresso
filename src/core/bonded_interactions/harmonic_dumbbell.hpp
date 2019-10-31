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
#ifndef _HARMONIC_DUMBBELL_HPP
#define _HARMONIC_DUMBBELL_HPP
/** \file
 *  Routines to calculate the harmonic dumbbell bond potential between particle
 *  pairs.
 *
 *  Implementation in \ref harmonic_dumbbell.cpp.
 */

#include "config.hpp"

#ifdef ROTATION
#include "bonded_interaction_data.hpp"

#include <utils/math/sqr.hpp>

/** Set the parameters for the harmonic dumbbell bond potential
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int harmonic_dumbbell_set_params(int bond_type, double k1, double k2, double r,
                                 double r_cut);

/** Compute the harmonic dumbbell bond force and torque.
 *  @param[in]  director  Director of the particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @return the force and torque
 */
inline boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d>>
harmonic_dumbbell_pair_force(Utils::Vector3d const &director,
                             Bonded_ia_parameters const &iaparams,
                             Utils::Vector3d const &dx) {
  auto const dist = dx.norm();

  if ((iaparams.p.harmonic_dumbbell.r_cut > 0.0) &&
      (dist > iaparams.p.harmonic_dumbbell.r_cut)) {
    return {};
  }

  auto const dr = dist - iaparams.p.harmonic_dumbbell.r;
  auto const normalizer = (dist > ROUND_ERROR_PREC) ? 1. / dist : 0.0;
  auto const fac = -iaparams.p.harmonic_dumbbell.k1 * dr * normalizer;
  auto const force = fac * dx;

  auto const dhat = dx * normalizer;
  auto const da = vector_product(dhat, director);
  auto const torque = iaparams.p.harmonic_dumbbell.k2 * da;

  return std::make_tuple(force, torque);
}

/** Compute the harmonic dumbbell bond energy.
 *  @param[in]  director  Director of the particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 */
inline boost::optional<double>
harmonic_dumbbell_pair_energy(Utils::Vector3d const &director,
                              Bonded_ia_parameters const &iaparams,
                              Utils::Vector3d const &dx) {
  auto const dist = dx.norm();

  if ((iaparams.p.harmonic_dumbbell.r_cut > 0.0) &&
      (dist > iaparams.p.harmonic_dumbbell.r_cut)) {
    return {};
  }

  auto const dhat = dx / dist;
  auto const da = vector_product(dhat, director);
  auto const torque = iaparams.p.harmonic_dumbbell.k2 * da;
  auto const diff = dhat - director;

  return 0.5 * iaparams.p.harmonic_dumbbell.k1 *
             Utils::sqr(dist - iaparams.p.harmonic.r) +
         0.5 * iaparams.p.harmonic_dumbbell.k2 * (torque * diff);
}

#endif // ROTATION

#endif


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

#ifndef THERMALIZED_DIST_H
#define THERMALIZED_DIST_H
/** \file
 *  Routines to thermalize the center of mass and distance of a particle pair.
 *
 *  Implementation in \ref thermalized_bond.cpp.
 */

/** number of thermalized bonds */
extern int n_thermalized_bonds;

#include "Particle.hpp"
#include "bonded_interaction_data.hpp"
#include "integrate.hpp"
#include "random.hpp"
#include "thermostat.hpp"

/** Set the parameters of a thermalized bond
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int thermalized_bond_set_params(int bond_type, double temp_com,
                                double gamma_com, double temp_distance,
                                double gamma_distance, double r_cut);

void thermalized_bond_update_params(double pref_scale);
void thermalized_bond_init();

/** Separately thermalizes the com and distance of a particle pair.
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        %Distance between the particles.
 *  @return the forces on @p p1 and @p p2
 */
inline boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d>>
thermalized_bond_forces(Particle const &p1, Particle const &p2,
                        Bonded_ia_parameters const &iaparams,
                        Utils::Vector3d const &dx) {
  // Bond broke?
  if (iaparams.p.thermalized_bond.r_cut > 0.0 &&
      dx.norm() > iaparams.p.thermalized_bond.r_cut) {
    return {};
  }

  auto const mass_tot = p1.p.mass + p2.p.mass;
  auto const mass_tot_inv = 1.0 / mass_tot;
  auto const sqrt_mass_tot = sqrt(mass_tot);
  auto const sqrt_mass_red = sqrt(p1.p.mass * p2.p.mass / mass_tot);
  auto const com_vel = mass_tot_inv * (p1.p.mass * p1.m.v + p2.p.mass * p2.m.v);
  auto const dist_vel = p2.m.v - p1.m.v;

  extern ThermalizedBondThermostat thermalized_bond;
  Utils::Vector3d force1{};
  Utils::Vector3d force2{};
  auto const noise = Random::noise_uniform<RNGSalt::THERMALIZED_BOND>(
      thermalized_bond.rng_get(), p1.p.identity, p2.p.identity);

  for (int i = 0; i < 3; i++) {
    double force_lv_com, force_lv_dist;

    // Langevin thermostat for center of mass
    if (iaparams.p.thermalized_bond.pref2_com > 0.0) {
      force_lv_com =
          -iaparams.p.thermalized_bond.pref1_com * com_vel[i] +
          sqrt_mass_tot * iaparams.p.thermalized_bond.pref2_com * noise[i];
    } else {
      force_lv_com = -iaparams.p.thermalized_bond.pref1_com * com_vel[i];
    }

    // Langevin thermostat for distance p1->p2
    if (iaparams.p.thermalized_bond.pref2_dist > 0.0) {
      force_lv_dist =
          -iaparams.p.thermalized_bond.pref1_dist * dist_vel[i] +
          sqrt_mass_red * iaparams.p.thermalized_bond.pref2_dist * noise[i];
    } else {
      force_lv_dist = -iaparams.p.thermalized_bond.pref1_dist * dist_vel[i];
    }
    // Add forces
    force1[i] = p1.p.mass * mass_tot_inv * force_lv_com - force_lv_dist;
    force2[i] = p2.p.mass * mass_tot_inv * force_lv_com + force_lv_dist;
  }

  return std::make_tuple(force1, force2);
}

#endif

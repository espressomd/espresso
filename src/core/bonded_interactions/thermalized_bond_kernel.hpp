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

#ifndef THERMALIZED_DIST_KERNEL_H
#define THERMALIZED_DIST_KERNEL_H

#include "thermalized_bond.hpp"

#include "Particle.hpp"
#include "random.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>

#include <cmath>
#include <tuple>

/** Separately thermalizes the com and distance of a particle pair.
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  dx        %Distance between the particles.
 *  @return the forces on @p p1 and @p p2
 */
inline boost::optional<std::tuple<Utils::Vector3d, Utils::Vector3d>>
ThermalizedBond::forces(Particle const &p1, Particle const &p2,
                        Utils::Vector3d const &dx) const {
  // Bond broke?
  if (r_cut > 0.0 && dx.norm() > r_cut) {
    return {};
  }

  auto const mass_tot = p1.mass() + p2.mass();
  auto const mass_tot_inv = 1.0 / mass_tot;
  auto const sqrt_mass_tot = sqrt(mass_tot);
  auto const sqrt_mass_red = sqrt(p1.mass() * p2.mass() / mass_tot);
  auto const com_vel = mass_tot_inv * (p1.mass() * p1.v() + p2.mass() * p2.v());
  auto const dist_vel = p2.v() - p1.v();

  extern ThermalizedBondThermostat thermalized_bond;
  Utils::Vector3d force1{};
  Utils::Vector3d force2{};
  auto const noise = Random::noise_uniform<RNGSalt::THERMALIZED_BOND>(
      thermalized_bond.rng_counter(), thermalized_bond.rng_seed(), p1.id(),
      p2.id());

  for (unsigned int i = 0u; i < 3u; ++i) {
    double force_lv_com, force_lv_dist;

    // Langevin thermostat for center of mass
    if (pref2_com > 0.0) {
      force_lv_com =
          -pref1_com * com_vel[i] + sqrt_mass_tot * pref2_com * noise[i];
    } else {
      force_lv_com = -pref1_com * com_vel[i];
    }

    // Langevin thermostat for distance p1->p2
    if (pref2_dist > 0.0) {
      force_lv_dist =
          -pref1_dist * dist_vel[i] + sqrt_mass_red * pref2_dist * noise[i];
    } else {
      force_lv_dist = -pref1_dist * dist_vel[i];
    }
    // Add forces
    force1[i] = p1.mass() * mass_tot_inv * force_lv_com - force_lv_dist;
    force2[i] = p2.mass() * mass_tot_inv * force_lv_com + force_lv_dist;
  }

  return std::make_tuple(force1, force2);
}

#endif

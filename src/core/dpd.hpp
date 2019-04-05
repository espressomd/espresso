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
#ifndef CORE_DPD_HPP
#define CORE_DPD_HPP
/** \file
 *  Routines to use dpd as thermostat or pair force
 *  T. Soddemann, B. Duenweg and K. Kremer, Phys. Rev. E 68, 046702 (2003)
 *
 *  Implementation in forces.cpp.
 */

#include "config.hpp"

#ifdef DPD

#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_data.hpp"
#include "random.hpp"

#include "utils/uniform.hpp"
#include <Random123/philox.h>

/** philox functiontality: increment, get/set */
void dpd_rng_counter_increment();

void dpd_heat_up();
void dpd_cool_down();
void dpd_switch_off();
int dpd_set_params(int part_type_a, int part_type_b, double gamma, double r_c,
                   int wf, double tgamma, double tr_c, int twf, uint64_t seed);
void dpd_init();
void dpd_update_params(double pref2_scale);

Vector3d dpd_pair_force(Particle const *p1, Particle const *p2,
                        IA_parameters *ia_params, double const *d, double dist,
                        double dist2);


/** Return a random 3d vector with the philox thermostat.
    Random numbers depend on
    1. dpd_rng_counter (initialized by seed) which is increased on
   integration
    2. Salt (decorrelates different counter)
    3. Particle ID (decorrelates particles, gets rid of seed-per-node)
*/
inline Vector3d v_noise(int particle_id) {

  using rng_type = r123::Philox4x64;
  using ctr_type = rng_type::ctr_type;
  using key_type = rng_type::key_type;

  ctr_type c{
      {dpd_rng_counter.value(), static_cast<uint64_t>(RNGSalt::DPD)}};

  auto f_random = [&c](int id) -> Vector3d {
    key_type k{{static_cast<uint32_t>(id)}};

    auto const noise = rng_type{}(c, k);

    using Utils::uniform;
    return Vector3d{uniform(noise[0]), uniform(noise[1]), uniform(noise[2])} -
           Vector3d::broadcast(0.5);
  };

  return f_random(particle_id);
}

#endif

#endif

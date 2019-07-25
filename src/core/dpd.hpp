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
#include "particle_data.hpp"

#include <utils/Vector.hpp>

struct IA_parameters;

struct DPDParameters {
  double gamma = 0.;
  double cutoff = -1.;
  int wf = 0;
  double pref = 0.0;
};

#include "utils/uniform.hpp"
#include "utils/u32_to_u64.hpp"
#include <Random123/philox.h>

/** philox functiontality: increment, get/set */
extern std::unique_ptr<Utils::Counter<uint64_t>> dpd_rng_counter;
bool dpd_is_seed_required();
void dpd_rng_counter_increment();
void dpd_set_rng_state(uint64_t counter);
uint64_t dpd_get_rng_state();


void dpd_heat_up();
void dpd_cool_down();
int dpd_set_params(int part_type_a, int part_type_b, double gamma, double r_c,
                   int wf, double tgamma, double tr_c, int twf);
void dpd_init();
void dpd_update_params(double pref2_scale);

Utils::Vector3d dpd_pair_force(Particle const *p1, Particle const *p2,
                               const IA_parameters *ia_params, double const *d,
                               double dist, double dist2);
Utils::Vector9d dpd_stress();

/** Return a random 4d vector with the philox thermostat.
    Random numbers depend on
    1. dpd_rng_counter (initialized by seed) which is increased on
   integration
    2. Salt (decorrelates different counter)
    3. Two particle IDs (order-independent, decorrelates particles, gets rid of seed-per-node)
*/
inline Vector4d dpd_noise(int pid1, int pid2) {

  using rng_type = r123::Philox4x64;
  using ctr_type = rng_type::ctr_type;
  using key_type = rng_type::key_type;

  ctr_type c{
      {dpd_rng_counter->value(), static_cast<uint64_t>(RNGSalt::SALT_DPD)}};

  uint64_t merged_ids;
  if (pid1 > pid2)  {
    merged_ids = Utils::u32_to_u64(pid1, pid2);
  } else {
    merged_ids = Utils::u32_to_u64(pid2, pid1);
  }
  key_type k{merged_ids};

  auto const noise = rng_type{}(c, k);

  using Utils::uniform;
  return Vector4d{uniform(noise[0]), uniform(noise[1]), uniform(noise[2]), uniform(noise[3])} -
         Vector4d::broadcast(0.5);
}

#endif
#endif

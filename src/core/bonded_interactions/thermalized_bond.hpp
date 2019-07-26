
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

#ifndef THERMALIZED_DIST_H
#define THERMALIZED_DIST_H
/** \file
 *  Routines to thermalize the com and distance of a particle pair.
 *  \ref forces.cpp
 */

/** number of thermalized bonds */
extern int n_thermalized_bonds;

#include "bonded_interaction_data.hpp"
#include "debug.hpp"
#include "integrate.hpp"
#include "random.hpp"

#include <Random123/philox.h>
#include <utils/Counter.hpp>
#include <utils/u32_to_u64.hpp>
#include <utils/uniform.hpp>
#include <utils/Vector.hpp>

extern std::unique_ptr<Utils::Counter<uint64_t>> thermalized_bond_rng_counter;

/** Philox RNG python interface */
bool thermalized_bond_is_seed_required();
void thermalized_bond_set_rng_state(uint64_t counter);
uint64_t thermalized_bond_get_rng_state();

/** Called in integration loop */
void thermalized_bond_rng_counter_increment();

/** Set the parameters of a thermalized bond
 *
 *  @retval ES_OK on success
 *  @retval ES_ERROR on error
 */
int thermalized_bond_set_params(int bond_type, double temp_com,
                                double gamma_com, double temp_distance,
                                double gamma_distance, double r_cut);

void thermalized_bond_heat_up();
void thermalized_bond_cool_down();
void thermalized_bond_update_params(double pref_scale);
void thermalized_bond_init();

/** Return a random 3d vector with the philox thermostat.
    Random numbers depend on
    1. rng_counter (initialized by seed) which is increased on
   integration
    2. Salt (decorrelates different counter)
    3. Particle ID, Partner ID 
*/

inline Utils::Vector3d v_noise(int particle_id, int partner_id) {

  using rng_type = r123::Philox4x64;
  using ctr_type = rng_type::ctr_type;
  using key_type = rng_type::key_type;

  ctr_type c{{thermalized_bond_rng_counter->value(),
              static_cast<uint64_t>(RNGSalt::THERMALIZED_BOND)}};

  /** Bond is stored on particle_id, so concatenation with the
      partner_id is unique. This will give the same RN for 
	  multiple thermalized bonds on the same particle pair, which
	  should not be allowed.
   */
  uint64_t merged_ids;
  auto const id1 = static_cast<uint32_t>(particle_id);
  auto const id2 = static_cast<uint32_t>(partner_id);
  merged_ids = Utils::u32_to_u64(id1, id2);
  key_type k{merged_ids};

  auto const noise = rng_type{}(c, k);

  using Utils::uniform;
  return Utils::Vector3d{uniform(noise[0]), uniform(noise[1]),
                         uniform(noise[2])} -
         Utils::Vector3d::broadcast(0.5);
}

/** Separately thermalizes the com and distance of a particle pair.
 *  @param[in]  p1        First particle.
 *  @param[in]  p2        Second particle.
 *  @param[in]  iaparams  Bonded parameters for the pair interaction.
 *  @param[in]  dx        Distance between the particles.
 *  @param[out] force1    Force on particle @p p1
 *  @param[out] force2    Force on particle @p p2
 *  @retval 1 if the bond is broken
 *  @retval 0 otherwise
 */
inline int calc_thermalized_bond_forces(const Particle *p1, const Particle *p2,
                                        const Bonded_ia_parameters *iaparams,
                                        const Utils::Vector3d &dx,
                                        double *force1, double *force2) {
  // Bond broke?
  if (iaparams->p.thermalized_bond.r_cut > 0.0 &&
      dx.norm() > iaparams->p.thermalized_bond.r_cut) {
    return 1;
  }

  double force_lv_com, force_lv_dist, com_vel, dist_vel;
  double mass_tot = p1->p.mass + p2->p.mass;
  double mass_tot_inv = 1.0 / mass_tot;
  double sqrt_mass_tot = sqrt(mass_tot);
  double sqrt_mass_red = sqrt(p1->p.mass * p2->p.mass / mass_tot);

  Utils::Vector3d noise = v_noise(p1->p.identity, p2->p.identity); 

  for (int i = 0; i < 3; i++) {

    // Langevin thermostat for center of mass
    com_vel =
        mass_tot_inv * (p1->p.mass * p1->m.v[i] + p2->p.mass * p2->m.v[i]);
    if (iaparams->p.thermalized_bond.pref2_com > 0.0) {
      force_lv_com = -iaparams->p.thermalized_bond.pref1_com * com_vel +
                     sqrt_mass_tot * iaparams->p.thermalized_bond.pref2_com *
                         noise[i];
    } else {
      force_lv_com = -iaparams->p.thermalized_bond.pref1_com * com_vel;
    }

    // Langevin thermostat for distance p1->p2
    dist_vel = p2->m.v[i] - p1->m.v[i];
    if (iaparams->p.thermalized_bond.pref2_dist > 0.0) {
      force_lv_dist = -iaparams->p.thermalized_bond.pref1_dist * dist_vel +
                      sqrt_mass_red * iaparams->p.thermalized_bond.pref2_dist *
                          noise[i];
    } else {
      force_lv_dist = -iaparams->p.thermalized_bond.pref1_dist * dist_vel;
    }
    // Add forces
    force1[i] = p1->p.mass * mass_tot_inv * force_lv_com - force_lv_dist;
    force2[i] = p2->p.mass * mass_tot_inv * force_lv_com + force_lv_dist;
  }

  ONEPART_TRACE(if (p1->p.identity == check_id) fprintf(
      stderr, "%d: OPT: THERMALIZED BOND f = (%.3e,%.3e,%.3e)\n", this_node,
      p1->f.f[0] + force1[0], p1->f.f[1] + force1[1], p1->f.f[2] + force1[2]));
  ONEPART_TRACE(if (p2->p.identity == check_id) fprintf(
      stderr, "%d: OPT: THERMALIZED BOND f = (%.3e,%.3e,%.3e)\n", this_node,
      p2->f.f[0] + force2[0], p2->f.f[1] + force2[1], p2->f.f[2] + force2[2]));
  return 0;
}

#endif

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
/** \file
 *  Implementation of dpd.hpp.
 */
#include "dpd.hpp"

#ifdef DPD
#include "communication.hpp"
#include "event.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "random.hpp"
#include "short_range_loop.hpp"
#include "thermostat.hpp"

#include <Random123/philox.h>
#include <boost/mpi/collectives/reduce.hpp>
#include <utils/NoOp.hpp>
#include <utils/constants.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/u32_to_u64.hpp>
#include <utils/uniform.hpp>

using Utils::Vector3d;

std::unique_ptr<Utils::Counter<uint64_t>> dpd_rng_counter;

/** Return a random 3d vector with the philox thermostat.
    Random numbers depend on
    1. dpd_rng_counter (initialized by seed) which is increased on
   integration
    2. Salt (decorrelates different counter)
    3. Two particle IDs (order-independent, decorrelates particles, gets rid of
   seed-per-node)
*/
Vector3d dpd_noise(uint32_t pid1, uint32_t pid2) {

  using rng_type = r123::Philox4x64;
  using ctr_type = rng_type::ctr_type;
  using key_type = rng_type::key_type;

  ctr_type c{
      {dpd_rng_counter->value(), static_cast<uint64_t>(RNGSalt::SALT_DPD)}};

  uint64_t merged_ids;
  auto const id1 = static_cast<uint32_t>(pid1);
  auto const id2 = static_cast<uint32_t>(pid2);

  if (id1 > id2) {
    merged_ids = Utils::u32_to_u64(id1, id2);
  } else {
    merged_ids = Utils::u32_to_u64(id2, id1);
  }
  key_type k{merged_ids};

  auto const noise = rng_type{}(c, k);

  using Utils::uniform;
  return Vector3d{uniform(noise[0]), uniform(noise[1]), uniform(noise[2])} -
         Vector3d::broadcast(0.5);
}

void mpi_bcast_dpd_rng_counter_slave(const uint64_t counter) {
  dpd_rng_counter = std::make_unique<Utils::Counter<uint64_t>>(counter);
}

REGISTER_CALLBACK(mpi_bcast_dpd_rng_counter_slave)

void mpi_bcast_dpd_rng_counter(const uint64_t counter) {
  mpi_call(mpi_bcast_dpd_rng_counter_slave, counter);
}

void dpd_rng_counter_increment() { dpd_rng_counter->increment(); }

bool dpd_is_seed_required() {
  /* Seed is required if rng is not initialized */
  return dpd_rng_counter == nullptr;
}

void dpd_set_rng_state(const uint64_t counter) {
  mpi_bcast_dpd_rng_counter(counter);
  dpd_rng_counter = std::make_unique<Utils::Counter<uint64_t>>(counter);
}

uint64_t dpd_get_rng_state() { return dpd_rng_counter->value(); }

void dpd_heat_up() {
  double pref_scale = sqrt(3);
  dpd_update_params(pref_scale);
}

void dpd_cool_down() {
  double pref_scale = 1.0 / sqrt(3);
  dpd_update_params(pref_scale);
}

int dpd_set_params(int part_type_a, int part_type_b, double gamma, double r_c,
                   int wf, double tgamma, double tr_c, int twf) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  data->dpd_radial = DPDParameters{
      gamma, r_c, wf, sqrt(24.0 * temperature * gamma / time_step)};
  data->dpd_trans = DPDParameters{
      tgamma, tr_c, twf, sqrt(24.0 * temperature * tgamma / time_step)};

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

void dpd_init() {
  for (int type_a = 0; type_a < max_seen_particle_type; type_a++) {
    for (int type_b = 0; type_b < max_seen_particle_type; type_b++) {
      auto data = get_ia_param(type_a, type_b);

      data->dpd_radial.pref =
          sqrt(24.0 * temperature * data->dpd_radial.gamma / time_step);
      data->dpd_trans.pref =
          sqrt(24.0 * temperature * data->dpd_trans.gamma / time_step);
    }
  }
}

void dpd_update_params(double pref_scale) {
  for (int type_a = 0; type_a < max_seen_particle_type; type_a++) {
    for (int type_b = 0; type_b < max_seen_particle_type; type_b++) {
      auto data = get_ia_param(type_a, type_b);

      data->dpd_radial.pref *= pref_scale;
      data->dpd_trans.pref *= pref_scale;
    }
  }
}

static double weight(int type, double r_cut, double r) {
  if (type == 0) {
    return 1.;
  }
  return 1. - r / r_cut;
}

Vector3d dpd_pair_force(DPDParameters const &params, Vector3d const &v,
                        double dist, Vector3d const &noise) {
  if (dist < params.cutoff) {
    auto const omega = weight(params.wf, params.cutoff, dist);
    auto const omega2 = Utils::sqr(omega);

    auto const f_d = params.gamma * omega2 * v;
    auto const f_r = params.pref * omega * noise;

    return f_r - f_d;
  }

  return {};
}

Utils::Vector3d dpd_pair_force(Particle const *const p1,
                               Particle const *const p2,
                               IA_parameters const *const ia_params,
                               Utils::Vector3d const &d, double dist,
                               double dist2) {
  if (ia_params->dpd_radial.cutoff <= 0.0 &&
      ia_params->dpd_trans.cutoff <= 0.0) {
    return {};
  }

  auto const v21 = p1->m.v - p2->m.v;
  auto const noise_vec =
      (ia_params->dpd_radial.pref > 0.0 || ia_params->dpd_trans.pref > 0.0)
          ? dpd_noise(p1->p.identity, p2->p.identity)
          : Vector3d{};

  auto const f_r = dpd_pair_force(ia_params->dpd_radial, v21, dist, noise_vec);
  auto const f_t = dpd_pair_force(ia_params->dpd_trans, v21, dist, noise_vec);

  /* Projection operator to radial direction */
  auto const P = tensor_product(d / dist2, d);
  /* This is equivalent to P * f_r + (1 - P) * f_t, but with
   * doing only one matrix-vector multiplication */
  auto const force = P * (f_r - f_t) + f_t;
  return force;
}

static auto dpd_viscous_stress_local() {
  on_observable_calc();

  Utils::Vector<Vector3d, 3> stress{};
  short_range_loop(
      Utils::NoOp{},
      [&stress](const Particle &p1, const Particle &p2, Distance const &d) {
        auto const v21 = p1.m.v - p2.m.v;

        auto ia_params = get_ia_param(p1.p.type, p2.p.type);
        auto const dist = std::sqrt(d.dist2);

        auto const f_r = dpd_pair_force(ia_params->dpd_radial, v21, dist, {});
        auto const f_t = dpd_pair_force(ia_params->dpd_trans, v21, dist, {});

        /* Projection operator to radial direction */
        auto const P = tensor_product(d.vec21 / d.dist2, d.vec21);
        /* This is equivalent to P * f_r + (1 - P) * f_t, but with
         * doing only one matrix-vector multiplication */
        auto const f = P * (f_r - f_t) + f_t;

        stress += tensor_product(d.vec21, f);
      });

  return stress;
}
REGISTER_CALLBACK_REDUCTION(dpd_viscous_stress_local, std::plus<>())

/**
 * @brief Viscous stress tensor of the DPD interaction.
 *
 * This calculates the total viscous stress contribution of the
 * DPD interaction. It contains only the dissipative contributions
 * of the interaction without noise. It's calculated as the
 * sum over all pair virials as
 * \f[
 *    \sigma^{\nu\mu} = V^{-1}\sum_i \sum_{j < i} r_{i,j}^{\nu} (- \gamma_{i,j}
 * v_{i,j})^{\mu} \f] where \f$\gamma_{i,j}\f$ is the (in general tensor valued)
 * DPD friction coefficient for particles i and j, \f$v_{i,j}\f$, \f$r_{i,j}\f$
 * are their relative velocity and distance and \f$V\f$ is the box volume.
 *
 * @return Stress tensor contribution.
 */
Utils::Vector9d dpd_stress() {
  auto const stress = mpi_call(Communication::Result::reduction, std::plus<>(),
                               dpd_viscous_stress_local);
  auto const box_l = box_geo.length();

  return Utils::Vector9d{stress[0][0], stress[0][1], stress[0][2],
                         stress[1][0], stress[1][1], stress[1][2],
                         stress[2][0], stress[2][1], stress[2][2]} /
         (box_l[0] * box_l[1] * box_l[2]);
}
#endif

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
/** \file
 *  Implementation of dpd.hpp.
 */
#include "config.hpp"

#ifdef DPD

#include "dpd.hpp"

#include "MpiCallbacks.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "interactions.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "random.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>
#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/matrix.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <functional>

using Utils::Vector3d;

/** Return a random uniform 3D vector with the Philox thermostat.
 *  Random numbers depend on
 *  1. dpd_rng_counter (initialized by seed) which is increased on integration
 *  2. Salt (decorrelates different counters)
 *  3. Two particle IDs (order-independent, decorrelates particles, gets rid of
 *     seed-per-node)
 */
Vector3d dpd_noise(uint32_t pid1, uint32_t pid2) {
  return Random::noise_uniform<RNGSalt::SALT_DPD>(
      dpd.rng_counter(), dpd.rng_seed(), (pid1 < pid2) ? pid2 : pid1,
      (pid1 < pid2) ? pid1 : pid2);
}

int dpd_set_params(int part_type_a, int part_type_b, double gamma, double k,
                   double r_c, int wf, double tgamma, double tr_c, int twf) {
  IA_parameters &ia_params = *get_ia_param_safe(part_type_a, part_type_b);

  ia_params.dpd_radial = DPDParameters{
      gamma, k, r_c, wf, sqrt(24.0 * temperature * gamma / time_step)};
  ia_params.dpd_trans = DPDParameters{
      tgamma, k, tr_c, twf, sqrt(24.0 * temperature * tgamma / time_step)};

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

void dpd_init() {
  for (int type_a = 0; type_a < max_seen_particle_type; type_a++) {
    for (int type_b = 0; type_b < max_seen_particle_type; type_b++) {
      IA_parameters &ia_params = *get_ia_param(type_a, type_b);

      ia_params.dpd_radial.pref =
          sqrt(24.0 * temperature * ia_params.dpd_radial.gamma / time_step);
      ia_params.dpd_trans.pref =
          sqrt(24.0 * temperature * ia_params.dpd_trans.gamma / time_step);
    }
  }
}

static double weight(int type, double r_cut, double k, double r) {
  if (type == 0) {
    return 1.;
  }
  return 1. - pow((r / r_cut), k);
}

Vector3d dpd_pair_force(DPDParameters const &params, Vector3d const &v,
                        double dist, Vector3d const &noise) {
  if (dist < params.cutoff) {
    auto const omega = weight(params.wf, params.cutoff, params.k, dist);
    auto const omega2 = Utils::sqr(omega);

    auto const f_d = params.gamma * omega2 * v;
    auto const f_r = params.pref * omega * noise;

    return f_r - f_d;
  }

  return {};
}

Utils::Vector3d dpd_pair_force(Particle const &p1, Particle const &p2,
                               IA_parameters const &ia_params,
                               Utils::Vector3d const &d, double dist,
                               double dist2) {
  if (ia_params.dpd_radial.cutoff <= 0.0 && ia_params.dpd_trans.cutoff <= 0.0) {
    return {};
  }

  auto const v21 = p1.m.v - p2.m.v;
  auto const noise_vec =
      (ia_params.dpd_radial.pref > 0.0 || ia_params.dpd_trans.pref > 0.0)
          ? dpd_noise(p1.p.identity, p2.p.identity)
          : Vector3d{};

  auto const f_r = dpd_pair_force(ia_params.dpd_radial, v21, dist, noise_vec);
  auto const f_t = dpd_pair_force(ia_params.dpd_trans, v21, dist, noise_vec);

  /* Projection operator to radial direction */
  auto const P = tensor_product(d / dist2, d);
  /* This is equivalent to P * f_r + (1 - P) * f_t, but with
   * doing only one matrix-vector multiplication */
  auto const force = P * (f_r - f_t) + f_t;
  return force;
}

static auto dpd_viscous_stress_local() {
  on_observable_calc();

  Utils::Matrix<double, 3, 3> stress{};
  cell_structure.non_bonded_loop(
      [&stress](const Particle &p1, const Particle &p2, Distance const &d) {
        auto const v21 = p1.m.v - p2.m.v;

        IA_parameters const &ia_params = *get_ia_param(p1.p.type, p2.p.type);
        auto const dist = std::sqrt(d.dist2);

        auto const f_r = dpd_pair_force(ia_params.dpd_radial, v21, dist, {});
        auto const f_t = dpd_pair_force(ia_params.dpd_trans, v21, dist, {});

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
  auto const volume = box_geo.volume();

  return Utils::flatten(stress) / volume;
}
#endif

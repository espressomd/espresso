/*
 * Copyright (C) 2010-2026 The ESPResSo project
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
#include <config/config.hpp>

#ifdef ESPRESSO_DPD

#include "dpd.hpp"

#include "BoxGeometry.hpp"
#include "cell_system/CellStructure.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "random.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/Vector.hpp>
#include <utils/math/tensor_product.hpp>
#include <utils/matrix.hpp>

#include <boost/mpi/collectives/reduce.hpp>

#include <cmath>
#include <functional>
#include <type_traits>

/** Return a random uniform 3D vector with the Philox thermostat.
 *  Random numbers depend on
 *  1. dpd_rng_counter (initialized by seed) which is increased on integration
 *  2. Salt (decorrelates different counters)
 *  3. Two particle IDs (order-independent, decorrelates particles, gets rid of
 *     seed-per-node)
 */
Utils::Vector3d dpd_noise(DPDThermostat const &dpd, int pid1, int pid2) {
  auto const pref = (pid1 < pid2) ? 1.0 : -1.0;
  return pref * Random::noise_uniform<RNGSalt::SALT_DPD>(
                    dpd.rng_counter(), dpd.rng_seed(),
                    (pid1 < pid2) ? pid2 : pid1, (pid1 < pid2) ? pid1 : pid2);
}

void InteractionsNonBonded::dpd_init(double kT, double time_step) {
  auto const max_type = get_max_seen_particle_type();
  for (int type_a = 0; type_a <= max_type; type_a++) {
    for (int type_b = type_a; type_b <= max_type; type_b++) {
      auto &ia_params = get_ia_param(type_a, type_b);

      ia_params.dpd.radial.pref =
          sqrt(24.0 * kT * ia_params.dpd.radial.gamma / time_step);
      ia_params.dpd.trans.pref =
          sqrt(24.0 * kT * ia_params.dpd.trans.gamma / time_step);
    }
  }
}

Utils::Vector3d dpd_pair_force(
    Utils::Vector3d const &p1_position, Utils::Vector3d const &p1_velocity,
    int p1_id, Utils::Vector3d const &p2_position,
    Utils::Vector3d const &p2_velocity, int p2_id, DPDThermostat const &dpd,
    BoxGeometry const &box_geo, IA_parameters const &ia_params,
    Utils::Vector3d const &d, double dist, double dist2) {
  if (ia_params.dpd.radial.cutoff <= 0.0 && ia_params.dpd.trans.cutoff <= 0.0) {
    return {};
  }

  auto const v21 = box_geo.velocity_difference(p1_position, p2_position,
                                               p1_velocity, p2_velocity);
  auto const noise_vec =
      (ia_params.dpd.radial.pref > 0.0 || ia_params.dpd.trans.pref > 0.0)
          ? dpd_noise(dpd, p1_id, p2_id)
          : Utils::Vector3d{};

  auto const f_r = dpd_pair_force(ia_params.dpd.radial, v21, dist, noise_vec);
  auto const f_t = dpd_pair_force(ia_params.dpd.trans, v21, dist, noise_vec);

  /* Projection operator to radial direction */
  auto const P = tensor_product(d / dist2, d);
  /* This is equivalent to P * f_r + (1 - P) * f_t, but with
   * doing only one matrix-vector multiplication */
  auto const force = P * (f_r - f_t) + f_t;
  return force;
}

static auto dpd_pressure_local(System::System &system) {
  auto const &box_geo = *system.box_geo;
  auto const &nonbonded_ias = *system.nonbonded_ias;
  auto &cell_structure = *system.cell_structure;
  auto &dpd = *system.thermostat->dpd;
  system.on_observable_calc();

  Utils::Matrix<double, 3, 3> pressure{};
  cell_structure.non_bonded_loop([&pressure, &box_geo, &nonbonded_ias,
                                  &dpd](Particle const &p1, Particle const &p2,
                                        Distance const &d) {
    auto const v21 =
        box_geo.velocity_difference(p1.pos(), p2.pos(), p1.v(), p2.v());

    auto const &ia_params = nonbonded_ias.get_ia_param(p1.type(), p2.type());
    auto const noise_vec =
        (ia_params.dpd.radial.pref > 0.0 || ia_params.dpd.trans.pref > 0.0)
            ? dpd_noise(dpd, p1.id(), p2.id())
            : Utils::Vector3d{};
    auto const dist = std::sqrt(d.dist2);

    auto const f_r = dpd_pair_force(ia_params.dpd.radial, v21, dist, noise_vec);
    auto const f_t = dpd_pair_force(ia_params.dpd.trans, v21, dist, noise_vec);

    /* Projection operator to radial direction */
    auto const P = tensor_product(d.vec21 / d.dist2, d.vec21);
    /* This is equivalent to P * f_r + (1 - P) * f_t, but with
     * doing only one matrix-vector multiplication */
    auto const f = P * (f_r - f_t) + f_t;

    pressure += tensor_product(d.vec21, f);
  });

  return pressure;
}

/**
 * @brief Pressure tensor contribution of the DPD interaction.
 *
 * This calculates the total contribution of the DPD interaction to
 * the pressure tensor, i.e. the dissipative and random (noise) forces.
 * It's calculated as the sum over all pair virials as
 * \f[
 *    P^{\nu\mu} = V^{-1}\sum_i \sum_{j < i} r_{i,j}^{\nu} F_{i,j}^{\mu}
 * \f]
 * where \f$F_{i,j}\f$ is the DPD force (dissipative plus noise) exerted
 * by particle j on particle i, \f$r_{i,j}\f$ is their distance vector
 * and \f$V\f$ is the box volume.
 *
 * @return Pressure tensor contribution.
 */
Utils::Vector9d dpd_pressure(System::System &system,
                             boost::mpi::communicator const &comm) {
  auto const local_pressure = dpd_pressure_local(system);
  std::remove_const_t<decltype(local_pressure)> global_pressure{};

  boost::mpi::reduce(comm, local_pressure, global_pressure, std::plus<>(), 0);

  return Utils::flatten(global_pressure) / system.box_geo->volume();
}

#endif // ESPRESSO_DPD

/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_NPT
#include "velocity_verlet_npt.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"
#include "thermostats/npt_inline.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives.hpp>

#include <cmath>
#include <functional>

static constexpr Utils::Vector3i nptgeom_dir{{1, 2, 4}};

/**
 * @brief Scale and communicate instantaneous NpT pressure and
 * propagate the conjugate momentum for volume.
 * @f$ p_{\epsilon}(t+dt) = p_{\epsilon}(t)
 *                         + 3 V (P_{\text{inst}} - P_{\text{ext}}) dt / 2 @f$
 */
static void
velocity_verlet_npt_propagate_p_eps(NptIsoParameters &nptiso,
                                    InstantaneousPressure &npt_inst_pressure,
                                    double time_step) {
  /* finalize derivation of p_inst */
  npt_inst_pressure.p_inst = {0., 0.};

  for (auto i = 0u; i < 3u; ++i) {
    if (nptiso.geometry & NptIsoParameters::nptgeom_dir[i]) {
      npt_inst_pressure.p_inst[0] +=
          npt_inst_pressure.p_vir[i] + npt_inst_pressure.p_vel[i];
      npt_inst_pressure.p_inst[1] += npt_inst_pressure.p_vir[i];
    }
  }

  Utils::Vector2d p_sum = {0., 0.};
  boost::mpi::reduce(::comm_cart, npt_inst_pressure.p_inst, p_sum,
                     std::plus<Utils::Vector2d>(), 0);

  /* propagate p_epsilon */
  if (::this_node == 0) {
    npt_inst_pressure.p_inst = p_sum / (nptiso.dimension * nptiso.volume);
    nptiso.p_epsilon += nptiso.volume *
                        (npt_inst_pressure.p_inst[0] - nptiso.p_ext) * 1.5 *
                        time_step;
    if (nptiso.particle_number > 1) {
      nptiso.p_epsilon +=
          (p_sum[0] - p_sum[1]) *
          (0.5 / static_cast<double>(nptiso.particle_number - 1)) * time_step;
    }
  }
  boost::mpi::broadcast(::comm_cart, nptiso.p_epsilon, 0);
}

/* propagate positions with dt/2;
 * @f$ pos[t+0.5*dt] = pos[t] + vel[t+0.5*dt] * 0.5 * dt @f$
 */
static void velocity_verlet_npt_propagate_pos(ParticleRangeNPT const &particles,
                                              double time_step) {
  for (auto &p : particles) {
    for (auto j = 0u; j < 3u; ++j) {
      if (!p.is_fixed_along(j)) {
        p.pos()[j] = p.pos()[j] + p.v()[j] * 0.5 * time_step;
      }
    }
  }
}

/* rescaling positions based on MTK equation;
 * @f$ pos[t] = \exp(0.5 * dt * p_{\epsilon} / W) * pos[t] @f$
 */
static void
velocity_verlet_npt_propagate_pos_MTK(NptIsoParameters &nptiso,
                                      ParticleRangeNPT const &particles) {
  auto const propagator =
      std::exp(nptiso.half_dt_inv_piston * nptiso.p_epsilon);

  for (auto &p : particles) {
    for (auto j = 0u; j < 3u; ++j) {
      if (!p.is_fixed_along(j)) {
        if (nptiso.geometry & NptIsoParameters::nptgeom_dir[j]) {
          p.pos()[j] *= propagator;
        }
      }
    }
  }
}

/**
 * @brief propagate positions and the volume and add thermal fluctuation.
 * A and V are the position and volume propagators for half-time step.
 * O is the propagator corresponding to Ornstein-Uhlenbeck process
 * representing the stochastic thermostat.
 * The time evolution follows the sequence A-V-O-V-A in this function,
 * with propagators applied right to left.
 */
static void velocity_verlet_npt_propagate_AVOVA_MTK(
    ParticleRangeNPT const &particles, IsotropicNptThermostat const &npt_iso,
    double time_step, System::System &system) {

  auto &box_geo = *system.box_geo;
  auto &cell_structure = *system.cell_structure;
  auto &nptiso = *system.nptiso;
  double L_new = 0.0;

  /* 1st propagation pos_MTK and pos*/
  velocity_verlet_npt_propagate_pos_MTK(nptiso, particles);
  velocity_verlet_npt_propagate_pos(particles, time_step);

  /* stochastic reservoirs for conjugate momentum for particles
   *  @f$ p(t+0.5*dt) = p(t+0.5*dt) \exp(- \gamma_0 dt / m)
   *                      + \sqrt{k_B T (1 - \exp(-2 \gamma_0 dt)}N(0,1) @f$
   */
  for (auto &p : particles) {
    auto const v_therm =
        propagate_therm0_nptiso(npt_iso, p.v(), p.mass(), p.id());
    for (unsigned int j = 0; j < 3; j++) {
      if (!p.is_fixed_along(j)) {
        if (nptiso.geometry & NptIsoParameters::nptgeom_dir[j]) {
          p.v()[j] = v_therm[j];
        }
      }
    }
  }

  /* 1st propagation: volume with dt/2
   * @f$ V(t+0.5*dt) = \exp(0.5 * dt * 3 * p_{\epsilon} / W) * V(t) @f$,
   * stochastic reservoirs for conjugate momentum for V
   * @f$ p_{\epsilon} = p_{epsilon}(t) \exp(- \gamma_V dt / W)
   *                    + \sqrt{k_B T (1 - \exp(-2 \gamma_V dt / W)}N(0,1) @f$,
   * 2nd propagation: volume with dt/2
   * @f$ V(t+dt) = \exp(0.5 * dt * 3 * p_{\epsilon} / W) * V(t+0.5*dt) @f$,
   */
  if (::this_node == 0) {
    nptiso.volume *=
        std::exp(1.5 * nptiso.inv_piston * nptiso.p_epsilon * time_step);
    nptiso.p_epsilon = propagate_thermV_nptiso(npt_iso, nptiso.p_epsilon);
    nptiso.volume *=
        std::exp(1.5 * nptiso.inv_piston * nptiso.p_epsilon * time_step);
    L_new = std::pow(nptiso.volume, 1.0 / nptiso.dimension);
  }

  /* 2nd propagation pos and pos_MTK*/
  velocity_verlet_npt_propagate_pos(particles, time_step);
  velocity_verlet_npt_propagate_pos_MTK(nptiso, particles);

  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);

  /* Apply new volume to the box-length, communicate it, and account for
   * necessary adjustments to the cell geometry */
  Utils::Vector3d new_box;

  if (::this_node == 0) {
    new_box = box_geo.length();
    for (auto i = 0u; i < 3u; ++i) {
      if (nptiso.cubic_box ||
          nptiso.geometry & NptIsoParameters::nptgeom_dir[i]) {
        new_box[i] = L_new;
      }
    }
  }

  boost::mpi::broadcast(::comm_cart, new_box, 0);

  box_geo.set_length(new_box);
  // fast box length update
  system.on_boxl_change(true);
}

/* rescaling momentum based on MTK equation;
 * @f$ vel[t] = \exp(-0.5 * dt * p_{\epsilon} * (1 + 1 / (N - 1)) / W) * vel[t]
 * @f$
 */
static void
velocity_verlet_npt_propagate_vel_MTK(NptIsoParameters const &nptiso,
                                      InstantaneousPressure &npt_inst_pressure,
                                      ParticleRangeNPT const &particles) {
  npt_inst_pressure.p_vel = {};
  auto const propagater =
      std::exp(nptiso.half_dt_inv_piston_and_Nf * nptiso.p_epsilon);

  for (auto &p : particles) {
    for (auto j = 0u; j < 3u; ++j) {
      if (!p.is_fixed_along(j)) {
        if (nptiso.geometry & ::nptgeom_dir[j]) {
          p.v()[j] *= propagater;
          npt_inst_pressure.p_vel[j] += Utils::sqr(p.v()[j]) * p.mass();
        }
      }
    }
  }
}

void velocity_verlet_npt_MTK_step_1(ParticleRangeNPT const &particles,
                                    IsotropicNptThermostat const &npt_iso,
                                    double time_step, System::System &system) {
  auto &nptiso = *system.nptiso;
  auto &npt_inst_pressure = *system.npt_inst_pressure;
  velocity_verlet_npt_propagate_vel_MTK(nptiso, npt_inst_pressure, particles);
  velocity_verlet_npt_propagate_p_eps(nptiso, npt_inst_pressure, time_step);
  velocity_verlet_npt_propagate_vel(nptiso, npt_inst_pressure, particles,
                                    time_step);
  velocity_verlet_npt_propagate_AVOVA_MTK(particles, npt_iso, time_step,
                                          system);
}

void velocity_verlet_npt_MTK_step_2(ParticleRangeNPT const &particles,
                                    double time_step, System::System &system) {
  auto &nptiso = *system.nptiso;
  auto &npt_inst_pressure = *system.npt_inst_pressure;
  velocity_verlet_npt_propagate_vel(nptiso, npt_inst_pressure, particles,
                                    time_step);
  velocity_verlet_npt_propagate_p_eps(nptiso, npt_inst_pressure, time_step);
  velocity_verlet_npt_propagate_vel_MTK(nptiso, npt_inst_pressure, particles);
}

#endif // ESPRESSO_NPT

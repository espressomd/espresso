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

#include "config/config.hpp"

#ifdef NPT
#include "velocity_verlet_npt.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "thermostats/npt_inline.hpp"

#include <utils/Vector.hpp>

#include <boost/mpi/collectives.hpp>

#include <cmath>
#include <functional>

/**
 * @brief Scale and communicate instantaneous NpT pressure and
 * propagate the conjugate momentum for volume.
 * @f$ p_{\epsilon}(t+dt) = p_{\epsilon}(t)
 *                          + (P_{\text{inst}} - P_{\text{ext}})*0.5*dt @f$
 */
static void
velocity_verlet_npt_finalize_p_inst(NptIsoParameters &nptiso,
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
    nptiso.p_epsilon +=
        (npt_inst_pressure.p_inst[0] - nptiso.p_ext) * 0.5 * time_step;
  }
}

/**
 * @brief propagete positions and the volume and add thermal fluctuation.
 * A and V are the position and volume propagators for half-time step.
 * O is the propagator corresponding to Ornstein-Uhlenbeck process
 * representing the stochastic thermostat.
 * The time evolution follows the sequence A-V-O-V-A in this function,
 * with propagators applied right to left.
 */
static void velocity_verlet_npt_propagate_AVOVA_And(
    ParticleRangeNPT const &particles, IsotropicNptThermostat const &npt_iso,
    double time_step, System::System &system) {

  auto &box_geo = *system.box_geo;
  auto &cell_structure = *system.cell_structure;
  auto &npt_inst_pressure = *system.npt_inst_pressure;
  auto &nptiso = *system.nptiso;
  Utils::Vector3d scal{};
  double L_halfdt = 0.0;
  double L_dt = 0.0;

  /* finalize derivation of p_inst and propagate p_epsilon;
   * @f$ p_{\epsilon}(t+dt) = p_{\epsilon}(t)
   *                          + (P_{\text{inst}} - P_{\text{ext}})*0.5*dt @f$
   */
  velocity_verlet_npt_finalize_p_inst(nptiso, npt_inst_pressure, time_step);

  /* 1st adjust \ref NptIsoParameters::nptiso.volume with dt/2;
   * @f$ V(t+0.5*dt) = V(t) + 0.5*p_{\epsilon}*dt @f$
   * and prepare pos and vel-rescaling
   */
  if (::this_node == 0) {
    nptiso.volume += nptiso.inv_piston * nptiso.p_epsilon * 0.5 * time_step;
    // L(t)**2 / L(t+0.5*dt)**2, where the numerator follows time in position,
    // and the denominator follows time in velocity
    scal[2] = Utils::sqr(box_geo.length()[nptiso.non_const_dim]) /
              pow(nptiso.volume, 2.0 / nptiso.dimension);
    if (nptiso.volume < 0.0) {
      runtimeErrorMsg()
          << "your choice of piston= " << nptiso.piston << ", dt= " << time_step
          << ", p_epsilon= " << nptiso.p_epsilon
          << " just caused the volume to become negative, decrease dt";
      nptiso.volume = box_geo.volume();
      scal[2] = 1;
    }

    L_halfdt = pow(nptiso.volume, 1.0 / nptiso.dimension);

    // L(t+0.5*dt) / L(t)
    scal[1] = L_halfdt * box_geo.length_inv()[nptiso.non_const_dim];
    // L(t) / L(t+0.5*dt)
    scal[0] = 1. / scal[1];
  }
  boost::mpi::broadcast(::comm_cart, scal, 0);

  /* 1st propagate positions with dt/2 and rescaling pos;
   * @f$ pos[t+0.5*dt] = (L(t+0.5*dt) / L(t)) *
   *    (pos[t] + (L(t)**2 / L(t+0.5*dt)**2) * vel[t+0.5*dt] * 0.5 * dt) @f$
   */
  for (auto &p : particles) {
    for (auto j = 0u; j < 3u; ++j) {
      if (!p.is_fixed_along(j)) {
        if (nptiso.geometry & NptIsoParameters::nptgeom_dir[j]) {
          p.pos()[j] =
              scal[1] * (p.pos()[j] + scal[2] * p.v()[j] * 0.5 * time_step);
          p.pos_at_last_verlet_update()[j] *= scal[1];
          p.v()[j] *= scal[0];
        } else {
          p.pos()[j] += p.v()[j] * 0.5 * time_step;
        }
      }
    }
  }

  /* stochastic reservoirs for conjugate momentum for V;
   * @f$ p_{\epsilon} = p_{epsilon}(t) \exp(- \gamma_V dt / W)
   *                    + \sqrt{k_B T (1 - \exp(-2 \gamma_V dt / W)}N(0,1) @f$,
   * 2nd adjust \ref NptIsoParameters::nptiso.volume with dt/2;
   * @f$ V(t+dt) = V(t+0.5*dt) + 0.5*p_{\epsilon}*dt @f$,
   * and prepare pos- and vel-rescaling
   */
  if (::this_node == 0) {
    nptiso.p_epsilon =
        propagate_thermV_nptiso(npt_iso, nptiso.p_epsilon, nptiso.piston);
    nptiso.volume += nptiso.inv_piston * nptiso.p_epsilon * 0.5 * time_step;
    L_dt = pow(nptiso.volume, 1.0 / nptiso.dimension);

    scal[2] = 1.0;
    scal[1] = L_dt / L_halfdt;
    scal[0] = 1. / scal[1];
  }
  boost::mpi::broadcast(::comm_cart, scal, 0);

  /* stochastic reserviors for velocities;
   *  @f$ p(t+0.5*dt) = p(t+0.5*dt) \exp(- \gamma_0 dt / m)
   *                      + \sqrt{k_B T (1 - \exp(-2 \gamma_0 dt)}N(0,1) @f$
   * and 2nd propagate positions with dt/2 while rescaling positions velocities
   *  @f$ pos[t+dt] = (L(t+dt) / L(t+0.5*dt)) *
   *                  (pos[t+0.5*dt] + vel[t+0.5*dt]) @f$
   */
  for (auto &p : particles) {
    auto const v_therm =
        propagate_therm0_nptiso(npt_iso, p.v(), p.mass(), p.id());
    for (auto j = 0u; j < 3u; ++j) {
      if (!p.is_fixed_along(j)) {
        if (nptiso.geometry & NptIsoParameters::nptgeom_dir[j]) {
          p.v()[j] = v_therm[j];
          p.pos()[j] =
              scal[1] * (p.pos()[j] + scal[2] * p.v()[j] * 0.5 * time_step);
          p.pos_at_last_verlet_update()[j] *= scal[1];
          p.v()[j] *= scal[0];
        } else {
          p.pos()[j] += p.v()[j] * 0.5 * time_step;
        }
      }
    }
  }

  cell_structure.set_resort_particles(Cells::RESORT_LOCAL);

  /* Apply new volume to the box-length, communicate it, and account for
   * necessary adjustments to the cell geometry */
  Utils::Vector3d new_box;

  if (::this_node == 0) {
    new_box = box_geo.length();
    for (auto i = 0u; i < 3u; ++i) {
      if (nptiso.cubic_box ||
          nptiso.geometry & NptIsoParameters::nptgeom_dir[i]) {
        new_box[i] = L_dt;
      }
    }
  }

  boost::mpi::broadcast(::comm_cart, new_box, 0);

  box_geo.set_length(new_box);
  // fast box length update
  system.on_boxl_change(true);
}

void velocity_verlet_npt_Andersen_step_1(ParticleRangeNPT const &particles,
                                         IsotropicNptThermostat const &npt_iso,
                                         double time_step,
                                         System::System &system) {
  auto &nptiso = *system.nptiso;
  auto &npt_inst_pressure = *system.npt_inst_pressure;
  velocity_verlet_npt_propagate_vel(nptiso, npt_inst_pressure, particles,
                                    time_step);
  velocity_verlet_npt_propagate_AVOVA_And(particles, npt_iso, time_step,
                                          system);
}

void velocity_verlet_npt_Andersen_step_2(ParticleRangeNPT const &particles,
                                         double time_step,
                                         System::System &system) {
  auto &nptiso = *system.nptiso;
  auto &npt_inst_pressure = *system.npt_inst_pressure;
  velocity_verlet_npt_propagate_vel_final(nptiso, npt_inst_pressure, particles,
                                          time_step);
  velocity_verlet_npt_finalize_p_inst(nptiso, npt_inst_pressure, time_step);
}

#endif // NPT

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

#include "lb/particle_coupling.hpp"
#include "BoxGeometry.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "errorhandling.hpp"
#include "random.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"

#include <utils/Counter.hpp>
#include <utils/Vector.hpp>

#include <boost/mpi.hpp>

#ifdef CALIPER
#include <caliper/cali.h>
#endif

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <optional>
#include <ranges>
#include <stdexcept>
#include <vector>

static Thermostat::GammaType lb_handle_particle_anisotropy(Particle const &p,
                                                           double lb_gamma) {
#ifdef THERMOSTAT_PER_PARTICLE
  auto const &partcl_gamma = p.gamma();
#ifdef PARTICLE_ANISOTROPY
  auto const default_gamma = Thermostat::GammaType::broadcast(lb_gamma);
#else
  auto const default_gamma = lb_gamma;
#endif // PARTICLE_ANISOTROPY
  return Thermostat::handle_particle_gamma(partcl_gamma, default_gamma);
#else
  return lb_gamma;
#endif // THERMOSTAT_PER_PARTICLE
}

static Utils::Vector3d lb_drag_force(Particle const &p, double lb_gamma,
                                     Utils::Vector3d const &v_fluid) {
#ifdef LB_ELECTROHYDRODYNAMICS
  auto const v_drift = v_fluid + p.mu_E();
#else
  auto const &v_drift = v_fluid;
#endif
  auto const gamma = lb_handle_particle_anisotropy(p, lb_gamma);

  /* calculate viscous force (eq. (9) @cite ahlrichs99a) */
  return Utils::hadamard_product(gamma, v_drift - p.v());
}

Utils::Vector3d lb_drag_force(LB::Solver const &lb, double lb_gamma,
                              Particle const &p,
                              Utils::Vector3d const &shifted_pos) {
  auto const v_fluid = lb.get_coupling_interpolated_velocity(shifted_pos);
  return lb_drag_force(p, lb_gamma, v_fluid);
}

/**
 * @brief Check if a position is within the local box + halo.
 *
 * @param local_box Local geometry
 * @param pos Position to check
 * @param halo Halo
 *
 * @return True iff the point is inside of the box up to halo.
 */
static bool in_local_domain(LocalBox const &local_box,
                            Utils::Vector3d const &pos, double halo = 0.) {
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const lower_corner = local_box.my_left() - halo_vec;
  auto const upper_corner = local_box.my_right() + halo_vec;

  return pos >= lower_corner and pos < upper_corner;
}

static bool in_box(Utils::Vector3d const &pos,
                   Utils::Vector3d const &lower_corner,
                   Utils::Vector3d const &upper_corner) {
  return pos >= lower_corner and pos < upper_corner;
}

bool in_local_halo(LocalBox const &local_box, Utils::Vector3d const &pos,
                   double agrid) {
  auto const halo = 0.5 * agrid;
  return in_local_domain(local_box, pos, halo);
}

static void positions_in_halo_impl(Utils::Vector3d const &pos_folded,
                                   Utils::Vector3d const &halo_lower_corner,
                                   Utils::Vector3d const &halo_upper_corner,
                                   BoxGeometry const &box_geo,
                                   std::vector<Utils::Vector3d> &res) {

  // Lees-Edwards: pre-calc positional offset folded into the simulation box
  double folded_le_offset = 0.;
  auto const &le = box_geo.lees_edwards_bc();
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    folded_le_offset = Algorithm::periodic_fold(
        le.pos_offset, box_geo.length()[le.shear_direction]);
  }

  for (int i : {-1, 0, 1}) {
    for (int j : {-1, 0, 1}) {
      for (int k : {-1, 0, 1}) {
        Utils::Vector3d shift{{double(i), double(j), double(k)}};

        auto pos_shifted = pos_folded;
        // Lees Edwards: folded position incl. LE pos offset
        // This is needed to ensure that the position from which `pos_shifted`
        // is calculated below, is always in the primary simulation box.
        if (box_geo.type() == BoxType::LEES_EDWARDS) {
          pos_shifted[le.shear_direction] = Algorithm::periodic_fold(
              pos_folded[le.shear_direction] +
                  shift[le.shear_plane_normal] * folded_le_offset,
              box_geo.length()[le.shear_direction]);
        }
        pos_shifted += Utils::hadamard_product(box_geo.length(), shift);

        if (in_box(pos_shifted, halo_lower_corner, halo_upper_corner)) {
          res.emplace_back(pos_shifted);
        }
      }
    }
  }
}

/**
 * @brief Return a vector of positions shifted by +,- box length in each
 * coordinate
 */
std::vector<Utils::Vector3d> positions_in_halo(Utils::Vector3d const &pos,
                                               BoxGeometry const &box_geo,
                                               LocalBox const &local_box,
                                               double agrid) {
  auto const halo = 0.5 * agrid;
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const fully_inside_lower = local_box.my_left() + 2. * halo_vec;
  auto const fully_inside_upper = local_box.my_right() - 2. * halo_vec;
  auto const pos_folded = box_geo.folded_position(pos);
  if (in_box(pos_folded, fully_inside_lower, fully_inside_upper)) {
    return {pos_folded};
  }
  auto const halo_lower_corner = local_box.my_left() - halo_vec;
  auto const halo_upper_corner = local_box.my_right() + halo_vec;
  std::vector<Utils::Vector3d> res;
  positions_in_halo_impl(pos_folded, halo_lower_corner, halo_upper_corner,
                         box_geo, res);
  return res;
}

static auto lees_edwards_vel_shift(Utils::Vector3d const &pos_shifted_by_box_l,
                                   Utils::Vector3d const &orig_pos,
                                   BoxGeometry const &box_geo) {
  Utils::Vector3d vel_shift{{0., 0., 0.}};
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    auto le = box_geo.lees_edwards_bc();
    auto normal_shift =
        (pos_shifted_by_box_l - orig_pos)[le.shear_plane_normal];
    // normal_shift is +,- box_l or 0 up to floating point errors
    if (normal_shift > std::numeric_limits<double>::epsilon()) {
      vel_shift[le.shear_direction] -= le.shear_velocity;
    }
    if (normal_shift < -std::numeric_limits<double>::epsilon()) {
      vel_shift[le.shear_direction] += le.shear_velocity;
    }
  }
  return vel_shift;
}

namespace LB {

Utils::Vector3d ParticleCoupling::get_noise_term(Particle const &p) const {
  if (not m_thermalized) {
    return Utils::Vector3d{};
  }
  using std::sqrt;
  using Utils::sqrt;
  auto const gamma = lb_handle_particle_anisotropy(p, m_thermostat.gamma);
  auto const noise = Random::noise_uniform<RNGSalt::PARTICLES>(
      m_thermostat.rng_counter(), m_thermostat.rng_seed(), p.id());
  return m_noise_pref_wo_gamma * Utils::hadamard_product(sqrt(gamma), noise);
}

/** @brief Collect particles to couple to LB and initiate obtaining Lb
 * interpolated velocities */
ParticleCouplingState
ParticleCoupling::prepare_coupling(std::vector<Particle *> const &particles) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  ParticleCouplingState state;
  state.coupled_particle_data.reserve(particles.size());

  if (particles.empty()) {
    return state;
  }

  auto const halo = 0.5 * m_lb.get_agrid();
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const fully_inside_lower = m_local_box.my_left() + 2. * halo_vec;
  auto const fully_inside_upper = m_local_box.my_right() - 2. * halo_vec;
  auto const halo_lower_corner = m_local_box.my_left() - halo_vec;
  auto const halo_upper_corner = m_local_box.my_right() + halo_vec;

  // First pass: determine positions and coupling modes
  for (auto ptr : particles) {
    auto &p = *ptr;
    auto const folded_pos = m_box_geo.folded_position(p.pos());

    ParticleCouplingState::CoupledParticleData data;
    data.particle = ptr;
    data.mode = ParticleCouplingState::none;
    data.force_positions_start = state.all_force_positions.size();
    data.force_positions_count = 0;

    // Temporary vector to collect positions for this particle
    std::vector<Utils::Vector3d> temp_positions;
    temp_positions.reserve(8);

    if (in_box(folded_pos, fully_inside_lower, fully_inside_upper)) {
      // Determine force coupling positions
      // If the folded position of the particle is further than 1/2 agrid away
      // from any wall, ghosts of the particle (shifted by +- box_l in any
      // coordinate) are not within the domain of this MPi rank including halos
      // (-agrid/2 +local_box_min..local_box_max+agrid/2)
      temp_positions.emplace_back(folded_pos);
    } else {
      // the particle is close to a boundary. We might have to also couple
      // halos of the particle shifted by +- box_l in any coordinate
      positions_in_halo_impl(folded_pos, halo_lower_corner, halo_upper_corner,
                             m_box_geo, temp_positions);
    }

    // Check if particle should be coupled
#ifdef ENGINE
    if (p.swimming().is_engine_force_on_fluid) {
      // swimmers excert a force on the fluid, but there is no
      // coupling force based on the velocity difference between particle and
      // fluid
      data.mode = ParticleCouplingState::swimmer_force_on_fluid;
    }
#endif
    if (data.mode == ParticleCouplingState::none) {
      // Check if any position is within velocity coupling region
      for (auto const &pos : temp_positions) {
        if (pos >= halo_lower_corner and pos < halo_upper_corner) {
          // Coupling force on particle and fluid based on velocity difference.
          // I.e., flwo velocity has to be obtained at this position.
          // Note: Only one interpolation per particle is needed
          // since velocities in the halo layers are identical
          data.velocity_coupling_index =
              state.positions_velocity_coupling.size();
          state.positions_velocity_coupling.emplace_back(pos);
          data.mode = ParticleCouplingState::particle_force;
          break;
        }
      }
    }

    // Only keep particles that will be coupled
    if (data.mode != ParticleCouplingState::none) {
      // Add positions directly to global vector
      state.all_force_positions.insert(state.all_force_positions.end(),
                                       temp_positions.begin(),
                                       temp_positions.end());
      data.force_positions_count = temp_positions.size();
      state.coupled_particle_data.emplace_back(std::move(data));
    }
  }

  if (!state.coupled_particle_data.empty()) {
    state.interpolated_velocities = m_lb.get_coupling_interpolated_velocities(
        state.positions_velocity_coupling);
  }

  return state;
}

void ParticleCoupling::apply_forces(ParticleCouplingState &state) {
  if (state.coupled_particle_data.empty()) {
    return;
  }

  auto const &domain_lower_corner = m_local_box.my_left();
  auto const &domain_upper_corner = m_local_box.my_right();

  // Pre-allocate forces vector with same size as positions
  std::vector<Utils::Vector3d> all_forces;
  all_forces.reserve(state.all_force_positions.size());

  // Calculate and apply forces
  for (auto &data : state.coupled_particle_data) {
    auto &p = *data.particle;
    Utils::Vector3d force_on_particle = {};

    if (data.mode == ParticleCouplingState::particle_force) {
#ifndef THERMOSTAT_PER_PARTICLE
      if (m_thermostat.gamma > 0.)
#endif
      {
        auto v_fluid =
            state.interpolated_velocities[*data.velocity_coupling_index];
        auto const &vel_pos =
            state.positions_velocity_coupling[*data.velocity_coupling_index];

        if (m_box_geo.type() == BoxType::LEES_EDWARDS) {
          // Account for the case where the interpolated velocity has been read
          // from a ghost of the particle across the LE boundary (or vice versa)
          // Then the particle velocity is shifted by +,- the LE shear velocity
          auto const vel_correction =
              lees_edwards_vel_shift(vel_pos, p.pos(), m_box_geo);
          v_fluid += vel_correction;
        }
        auto const drag_force = lb_drag_force(p, m_thermostat.gamma, v_fluid);
        auto const random_force = get_noise_term(p);
        force_on_particle = drag_force + random_force;
      }
    }

    auto force_on_fluid = -force_on_particle;
#ifdef ENGINE
    if (data.mode == ParticleCouplingState::swimmer_force_on_fluid) {
      force_on_fluid = p.calc_director() * p.swimming().f_swim;
    }
#endif

    // Apply forces to particle and collect forces for LB
    bool particle_force_applied = false;
    for (size_t i = 0; i < data.force_positions_count; ++i) {
      auto const &pos =
          state.all_force_positions[data.force_positions_start + i];
      if (pos >= domain_lower_corner and pos < domain_upper_corner) {
        /* Particle is in our LB volume, so this node
         * is responsible to adding its force */
        if (!particle_force_applied) {
          p.force() += force_on_particle;
          particle_force_applied = true;
        }
      }
      all_forces.emplace_back(force_on_fluid);
    }
  }

  m_lb.add_forces_at_pos(state.all_force_positions, all_forces);
}

#if defined(THERMOSTAT_PER_PARTICLE) and defined(PARTICLE_ANISOTROPY)
static void lb_coupling_sanity_checks(Particle const &p) {
  /*
  lb does (at the moment) not support rotational particle coupling.
  Consequently, anisotropic particles are also not supported.
  */
  auto const &p_gamma = p.gamma();
  if (p_gamma[0] != p_gamma[1] or p_gamma[1] != p_gamma[2]) {
    runtimeErrorMsg() << "anisotropic particle (id " << p.id()
                      << ") coupled to LB.";
  }
}
#endif

} // namespace LB

LB::ParticleCouplingState System::System::lb_prepare_particle_coupling() {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  assert(thermostat->lb != nullptr);

  if (thermostat->lb->couple_to_md) {
    if (not lb.is_solver_set()) {
      runtimeErrorMsg() << "The LB thermostat requires a LB fluid";
      return {};
    }
    auto const real_particles = cell_structure->local_particles();
    auto const ghost_particles = cell_structure->ghost_particles();
    LB::ParticleCoupling coupling{*thermostat->lb, lb, *box_geo, *local_geo};
    LB::CouplingBookkeeping bookkeeping{*cell_structure};
    lb.ghost_communication_vel();
    std::vector<Particle *> particles{};
    for (auto const *particle_range : {&real_particles, &ghost_particles}) {
      for (auto &p : *particle_range) {
        if (not LB::is_tracer(p) and bookkeeping.should_be_coupled(p)) {
#if defined(THERMOSTAT_PER_PARTICLE) and defined(PARTICLE_ANISOTROPY)
          LB::lb_coupling_sanity_checks(p);
#endif
          particles.emplace_back(&p);
        }
      }
    }

    return coupling.prepare_coupling(particles);
  }
  return {};
}

void System::System::lb_apply_particle_forces(
    LB::ParticleCouplingState &state) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  assert(thermostat->lb != nullptr);
  if (thermostat->lb->couple_to_md && !state.coupled_particle_data.empty()) {
    if (not lb.is_solver_set()) {
      runtimeErrorMsg() << "The LB thermostat requires a LB fluid";
      return;
    }
    LB::ParticleCoupling coupling{*thermostat->lb, lb, *box_geo, *local_geo};
    coupling.apply_forces(state);
  }
}

void System::System::lb_couple_particles() {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  auto coupling_state = lb_prepare_particle_coupling();
  lb_apply_particle_forces(coupling_state);
}

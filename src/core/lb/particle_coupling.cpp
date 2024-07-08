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
  for (int i : {-1, 0, 1}) {
    for (int j : {-1, 0, 1}) {
      for (int k : {-1, 0, 1}) {
        Utils::Vector3d shift{{double(i), double(j), double(k)}};
        Utils::Vector3d pos_shifted =
            pos_folded + Utils::hadamard_product(box_geo.length(), shift);

        if (box_geo.type() == BoxType::LEES_EDWARDS) {
          // note: modulo convention: 0 <= offset < length
          auto const &le = box_geo.lees_edwards_bc();
          auto const length = box_geo.length()[le.shear_direction];
          auto folded_offset = std::fmod(le.pos_offset, length);
          if (folded_offset < 0.) {
            folded_offset += length;
          }
          pos_shifted[le.shear_direction] +=
              shift[le.shear_plane_normal] * folded_offset;
        }

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

void ParticleCoupling::kernel(std::vector<Particle *> const &particles) {
  if (particles.empty()) {
    return;
  }
  enum coupling_modes { none, particle_force, swimmer_force_on_fluid };
  auto const halo = 0.5 * m_lb.get_agrid();
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const fully_inside_lower = m_local_box.my_left() + 2. * halo_vec;
  auto const fully_inside_upper = m_local_box.my_right() - 2. * halo_vec;
  auto const halo_lower_corner = m_local_box.my_left() - halo_vec;
  auto const halo_upper_corner = m_local_box.my_right() + halo_vec;
  std::vector<Utils::Vector3d> positions_velocity_coupling;
  std::vector<Utils::Vector3d> positions_force_coupling;
  std::vector<Utils::Vector3d> force_coupling_forces;
  std::vector<uint8_t> positions_force_coupling_counter;
  std::vector<Particle *> coupled_particles;
  for (auto ptr : particles) {
    auto &p = *ptr;
    auto span_size = 1u;
    auto const folded_pos = m_box_geo.folded_position(p.pos());
    if (in_box(folded_pos, fully_inside_lower, fully_inside_upper)) {
      positions_force_coupling.emplace_back(folded_pos);
    } else {
      auto const old_size = positions_force_coupling.size();
      positions_in_halo_impl(folded_pos, halo_lower_corner, halo_upper_corner,
                             m_box_geo, positions_force_coupling);
      auto const new_size = positions_force_coupling.size();
      span_size = static_cast<uint8_t>(new_size - old_size);
    }
    auto coupling_mode = none;
#ifdef ENGINE
    if (p.swimming().is_engine_force_on_fluid) {
      coupling_mode = swimmer_force_on_fluid;
    }
#endif
    if (coupling_mode == none) {
      for (auto end = positions_force_coupling.end(), it = end - span_size;
           it != end; ++it) {
        auto const &pos = *it;
        if (pos >= halo_lower_corner and pos < halo_upper_corner) {
          positions_velocity_coupling.emplace_back(pos);
          coupling_mode = particle_force;
          break;
        }
      }
    }
    if (coupling_mode == none) {
      positions_force_coupling.erase(positions_force_coupling.end() - span_size,
                                     positions_force_coupling.end());
    } else {
      coupled_particles.emplace_back(ptr);
      positions_force_coupling_counter.emplace_back(span_size);
    }
  }

  if (coupled_particles.empty()) {
    return;
  }
  auto interpolated_velocities =
      m_lb.get_coupling_interpolated_velocities(positions_velocity_coupling);

  auto const &domain_lower_corner = m_local_box.my_left();
  auto const &domain_upper_corner = m_local_box.my_right();
  auto it_interpolated_velocities = interpolated_velocities.begin();
  auto it_positions_force_coupling = positions_force_coupling.begin();
  auto it_positions_velocity_coupling = positions_velocity_coupling.begin();
  auto it_positions_force_coupling_counter =
      positions_force_coupling_counter.begin();
  for (auto ptr : coupled_particles) {
    auto &p = *ptr;
    auto coupling_mode = particle_force;
#ifdef ENGINE
    if (p.swimming().is_engine_force_on_fluid) {
      coupling_mode = swimmer_force_on_fluid;
    }
#endif
    Utils::Vector3d force_on_particle = {};
    if (coupling_mode == particle_force) {
      auto &v_fluid = *it_interpolated_velocities;
      if (m_box_geo.type() == BoxType::LEES_EDWARDS) {
        // Account for the case where the interpolated velocity has been read
        // from a ghost of the particle across the LE boundary (or vice verssa)
        // Then the particle velocity is shifted by +,- the LE shear velocity
        auto const vel_correction = lees_edwards_vel_shift(
            *it_positions_velocity_coupling, p.pos(), m_box_geo);
        v_fluid += vel_correction;
      }
      auto const drag_force = lb_drag_force(p, m_thermostat.gamma, v_fluid);
      auto const random_force = get_noise_term(p);
      force_on_particle = drag_force + random_force;
      ++it_interpolated_velocities;
      ++it_positions_velocity_coupling;
    }

    auto force_on_fluid = -force_on_particle;
#ifdef ENGINE
    if (coupling_mode == swimmer_force_on_fluid) {
      force_on_fluid = p.calc_director() * p.swimming().f_swim;
    }
#endif

    auto const span_size = *it_positions_force_coupling_counter;
    ++it_positions_force_coupling_counter;
    for (uint8_t i{0u}; i < span_size; ++i) {
      auto &pos = *it_positions_force_coupling;
      if (pos >= domain_lower_corner and pos < domain_upper_corner) {
        /* Particle is in our LB volume, so this node
         * is responsible to adding its force */
        p.force() += force_on_particle;
      }
      force_coupling_forces.emplace_back(force_on_fluid);
      ++it_positions_force_coupling;
    }
  }
  m_lb.add_forces_at_pos(positions_force_coupling, force_coupling_forces);
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

void System::System::lb_couple_particles() {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  assert(thermostat->lb != nullptr);
  if (thermostat->lb->couple_to_md) {
    if (not lb.is_solver_set()) {
      runtimeErrorMsg() << "The LB thermostat requires a LB fluid";
      return;
    }
    auto const real_particles = cell_structure->local_particles();
    auto const ghost_particles = cell_structure->ghost_particles();
    LB::ParticleCoupling coupling{*thermostat->lb, lb, *box_geo, *local_geo};
    LB::CouplingBookkeeping bookkeeping{*cell_structure};
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
    coupling.kernel(particles);
  }
}

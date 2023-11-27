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
#include <boost/serialization/optional.hpp>

#ifdef CALIPER
#include <caliper/cali.h>
#endif

#include <cmath>
#include <cstdint>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <vector>

LB::ParticleCouplingConfig lb_particle_coupling;

static auto is_lb_active() { return System::get_system().lb.is_solver_set(); }

void mpi_bcast_lb_particle_coupling_local() {
  boost::mpi::broadcast(comm_cart, lb_particle_coupling, 0);
}

REGISTER_CALLBACK(mpi_bcast_lb_particle_coupling_local)

void mpi_bcast_lb_particle_coupling() {
  mpi_call_all(mpi_bcast_lb_particle_coupling_local);
}

void lb_lbcoupling_activate() { lb_particle_coupling.couple_to_md = true; }

void lb_lbcoupling_deactivate() {
  if (is_lb_active() and this_node == 0 and lb_particle_coupling.gamma > 0.) {
    runtimeWarningMsg()
        << "Recalculating forces, so the LB coupling forces are not "
           "included in the particle force the first time step. This "
           "only matters if it happens frequently during sampling.";
  }

  lb_particle_coupling.couple_to_md = false;
}

void lb_lbcoupling_set_gamma(double gamma) {
  lb_particle_coupling.gamma = gamma;
}

double lb_lbcoupling_get_gamma() { return lb_particle_coupling.gamma; }

bool lb_lbcoupling_is_seed_required() {
  if (is_lb_active()) {
    return not lb_particle_coupling.rng_counter_coupling.is_initialized();
  }
  return false;
}

uint64_t lb_coupling_get_rng_state_cpu() {
  return lb_particle_coupling.rng_counter_coupling->value();
}

uint64_t lb_lbcoupling_get_rng_state() {
  if (is_lb_active()) {
    return lb_coupling_get_rng_state_cpu();
  }
  throw std::runtime_error("No LB active");
}

void lb_lbcoupling_set_rng_state(uint64_t counter) {
  if (is_lb_active()) {
    lb_particle_coupling.rng_counter_coupling =
        Utils::Counter<uint64_t>(counter);
  } else
    throw std::runtime_error("No LB active");
}

void add_md_force(LB::Solver &lb, Utils::Vector3d const &pos,
                  Utils::Vector3d const &force, double time_step) {
  /* transform momentum transfer to lattice units
     (eq. (12) @cite ahlrichs99a) */
  auto const delta_j = (time_step / lb.get_lattice_speed()) * force;
  lb.add_force_density(pos, delta_j);
}

static Thermostat::GammaType lb_handle_particle_anisotropy(Particle const &p) {
  auto const lb_gamma = lb_lbcoupling_get_gamma();
#ifdef THERMOSTAT_PER_PARTICLE
  auto const &partcl_gamma = p.gamma();
#ifdef PARTICLE_ANISOTROPY
  auto const default_gamma = Thermostat::GammaType::broadcast(lb_gamma);
#else
  auto const default_gamma = lb_gamma;
#endif // PARTICLE_ANISOTROPY
  return handle_particle_gamma(partcl_gamma, default_gamma);
#else
  return lb_gamma;
#endif // THERMOSTAT_PER_PARTICLE
}

Utils::Vector3d lb_drag_force(LB::Solver const &lb, Particle const &p,
                              Utils::Vector3d const &shifted_pos,
                              Utils::Vector3d const &vel_offset) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation (eq. (11) @cite ahlrichs99a) */
  auto const v_fluid = lb.get_coupling_interpolated_velocity(shifted_pos);
  auto const v_drift = v_fluid + vel_offset;
  auto const gamma = lb_handle_particle_anisotropy(p);

  /* calculate viscous force (eq. (9) @cite ahlrichs99a) */
  return Utils::hadamard_product(gamma, v_drift - p.v());
}

/**
 * @brief Check if a position is within the local box + halo.
 *
 * @param pos Position to check
 * @param halo Halo
 *
 * @return True iff the point is inside of the box up to halo.
 */
static bool in_local_domain(Utils::Vector3d const &pos, double halo = 0.) {
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const &local_geo = *System::get_system().local_geo;
  auto const lower_corner = local_geo.my_left() - halo_vec;
  auto const upper_corner = local_geo.my_right() + halo_vec;

  return pos >= lower_corner and pos < upper_corner;
}

static bool in_box(Utils::Vector3d const &pos,
                   Utils::Vector3d const &lower_corner,
                   Utils::Vector3d const &upper_corner) {
  return pos >= lower_corner and pos < upper_corner;
}

static bool in_local_halo(Utils::Vector3d const &pos, double agrid) {
  auto const halo = 0.5 * agrid;
  return in_local_domain(pos, halo);
}

bool in_local_halo(Utils::Vector3d const &pos) {
  return in_local_domain(pos, System::get_system().lb.get_agrid());
}

/**
 * @brief Return a vector of positions shifted by +,- box length in each
 * coordinate
 */
std::vector<Utils::Vector3d> positions_in_halo(Utils::Vector3d const &pos,
                                               BoxGeometry const &box,
                                               double agrid) {
  auto const &system = System::get_system();
  auto const &box_geo = *system.box_geo;
  auto const &local_geo = *system.local_geo;
  auto const halo = 0.5 * agrid;
  auto const halo_vec = Utils::Vector3d::broadcast(halo);
  auto const fully_inside_lower = local_geo.my_left() + 2. * halo_vec;

  // If the particle is at least one agrid away from the node boundary
  // any ghosts shifte dby +- box_length cannot be in the lb volume
  // accessible by this node (-agrid/2 to box_length +agrid/2)
  auto const fully_inside_upper = local_geo.my_right() - 2. * halo_vec;
  if (in_box(pos, fully_inside_lower, fully_inside_upper)) {
    return {pos};
  }

  // For remaingin particles, positoins shifted by +- box_length
  // can potentially be in the locally accessible LB volume
  auto const halo_lower_corner = local_geo.my_left() - halo_vec;
  auto const halo_upper_corner = local_geo.my_right() + halo_vec;

  std::vector<Utils::Vector3d> res;
  for (int i : {-1, 0, 1}) {
    for (int j : {-1, 0, 1}) {
      for (int k : {-1, 0, 1}) {
        Utils::Vector3d shift{{double(i), double(j), double(k)}};
        Utils::Vector3d pos_shifted =
            pos + Utils::hadamard_product(box.length(), shift);
        // Apply additional shift from Lees-Edwards boundary conditions
        if (box_geo.type() == BoxType::LEES_EDWARDS) {
          auto le = box_geo.lees_edwards_bc();
          auto normal_shift = (pos_shifted - pos)[le.shear_plane_normal];
          if (normal_shift > std::numeric_limits<double>::epsilon())
            pos_shifted[le.shear_direction] += le.pos_offset;
          if (normal_shift < -std::numeric_limits<double>::epsilon())
            pos_shifted[le.shear_direction] -= le.pos_offset;
        }

        if (in_box(pos_shifted, halo_lower_corner, halo_upper_corner)) {
          res.push_back(pos_shifted);
        }
      }
    }
  }
  return res;
}

//* Determine Lees-Edwards velocity offset for positions shifted
/*  by +- box_length in one or more coordinates,
/*  i.e., those obtained form positoins_in_halo()
*/
Utils::Vector3d lees_edwards_vel_shift(const Utils::Vector3d &shifted_pos,
                                       const Utils::Vector3d &orig_pos,
                                       const BoxGeometry &box_geo) {
  Utils::Vector3d vel_shift{{0., 0., 0.}};
  if (box_geo.type() == BoxType::LEES_EDWARDS) {
    auto le = box_geo.lees_edwards_bc();
    auto normal_shift = (shifted_pos - orig_pos)[le.shear_plane_normal];
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
  auto const &rng_counter = lb_particle_coupling.rng_counter_coupling;
  if (not rng_counter) {
    throw std::runtime_error(
        "Access to uninitialized LB particle coupling RNG counter");
  }
  using std::sqrt;
  using Utils::sqrt;
  auto const counter = rng_counter->value();
  auto const gamma = lb_handle_particle_anisotropy(p);
  return m_noise_pref_wo_gamma *
         Utils::hadamard_product(
             sqrt(gamma),
             Random::noise_uniform<RNGSalt::PARTICLES>(counter, 0, p.id()));
}

void ParticleCoupling::kernel(Particle &p) {
  if (p.is_virtual() and not m_couple_virtual)
    return;

  auto const agrid = m_lb.get_agrid();
  auto const &box_geo = *System::get_system().box_geo;

  // Calculate coupling force
  Utils::Vector3d force_on_particle = {};
  auto folded_pos = box_geo.folded_position(p.pos());

#ifdef ENGINE
  if (not p.swimming().is_engine_force_on_fluid)
#endif
    for (auto pos : positions_in_halo(folded_pos, box_geo, agrid)) {
      if (in_local_halo(pos, agrid)) {
        auto const vel_offset = lb_drift_velocity_offset(p) +
                                lees_edwards_vel_shift(pos, p.pos(), box_geo);
        auto const drag_force = lb_drag_force(m_lb, p, pos, vel_offset);
        auto const random_force = get_noise_term(p);
        force_on_particle = drag_force + random_force;
        break;
      }
    }

  auto force_on_fluid = -force_on_particle;
#ifdef ENGINE
  if (p.swimming().is_engine_force_on_fluid) {
    force_on_fluid = p.calc_director() * p.swimming().f_swim;
  }
#endif

  // couple positions including shifts by one box length to add
  // forces to ghost layers
  for (auto pos : positions_in_halo(folded_pos, box_geo, agrid)) {
    if (in_local_domain(pos)) {
      /* Particle is in our LB volume, so this node
       * is responsible to adding its force */
      p.force() += force_on_particle;
    }
    add_md_force(m_lb, pos, force_on_fluid, m_time_step);
  }
}

bool CouplingBookkeeping::is_ghost_for_local_particle(Particle const &p) const {
  auto const &system = System::get_system();
  auto const &cell_structure = *system.cell_structure;
  return not cell_structure.get_local_particle(p.id())->is_ghost();
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

void couple_particles(bool couple_virtual, ParticleRange const &real_particles,
                      ParticleRange const &ghost_particles, double time_step) {
#ifdef CALIPER
  CALI_CXX_MARK_FUNCTION;
#endif
  if (lb_particle_coupling.couple_to_md) {
    auto &lb = System::get_system().lb;
    if (lb.is_solver_set()) {
      ParticleCoupling coupling{lb, couple_virtual, time_step};
      CouplingBookkeeping bookkeeping{};
      for (auto const &particle_range : {real_particles, ghost_particles}) {
        for (auto &p : particle_range) {
          if (bookkeeping.should_be_coupled(p)) {
#if defined(THERMOSTAT_PER_PARTICLE) and defined(PARTICLE_ANISOTROPY)
            lb_coupling_sanity_checks(p);
#endif
            coupling.kernel(p);
          }
        }
      }
    }
  }
}

} // namespace LB

void lb_lbcoupling_propagate() {
  auto const &lb = System::get_system().lb;
  if (lb.is_solver_set() and lb.get_kT() > 0.0) {
    lb_particle_coupling.rng_counter_coupling->increment();
  }
}

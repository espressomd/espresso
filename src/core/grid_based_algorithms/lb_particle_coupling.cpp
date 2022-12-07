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
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "random.hpp"

#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_interpolation.hpp"
#include "grid_based_algorithms/lb_particle_coupling.hpp"

#include <profiler/profiler.hpp>
#include <utils/Counter.hpp>
#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/mpi.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <utility>

LB_Particle_Coupling lb_particle_coupling;

void mpi_bcast_lb_particle_coupling_local() {
  boost::mpi::broadcast(comm_cart, lb_particle_coupling, 0);
}

REGISTER_CALLBACK(mpi_bcast_lb_particle_coupling_local)

void mpi_bcast_lb_particle_coupling() {
  mpi_call_all(mpi_bcast_lb_particle_coupling_local);
}

void lb_lbcoupling_activate() { lb_particle_coupling.couple_to_md = true; }

void lb_lbcoupling_deactivate() {
  if (lattice_switch != ActiveLB::NONE && this_node == 0 &&
      lb_particle_coupling.gamma > 0.) {
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
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
    return not lb_particle_coupling.rng_counter_coupling.is_initialized();
  }
  return false;
}

uint64_t lb_coupling_get_rng_state_cpu() {
  return lb_particle_coupling.rng_counter_coupling->value();
}

uint64_t lb_lbcoupling_get_rng_state() {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
    return lb_coupling_get_rng_state_cpu();
  }
  throw std::runtime_error("No LB active");
}

void lb_lbcoupling_set_rng_state(uint64_t counter) {
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
    lb_particle_coupling.rng_counter_coupling =
        Utils::Counter<uint64_t>(counter);
  } else
    throw std::runtime_error("No LB active");
}

void add_md_force(Utils::Vector3d const &pos, Utils::Vector3d const &force,
                  double time_step) {
  /* transform momentum transfer to lattice units
     (eq. (12) @cite ahlrichs99a) */
  auto const delta_j = -(time_step / LB::get_lattice_speed()) * force;
  lb_lbinterpolation_add_force_density(pos, delta_j);
}

Utils::Vector3d lb_particle_coupling_drift_vel_offset(const Particle &p) {
  Utils::Vector3d vel_offset{};
#ifdef ENGINE
  if (p.swimming().swimming) {
    vel_offset += p.swimming().v_swim * p.calc_director();
  }
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  vel_offset += p.mu_E();
#endif
  return vel_offset;
}

Utils::Vector3d lb_drag_force(Particle const &p,
                              Utils::Vector3d const &shifted_pos,
                              Utils::Vector3d const &vel_offset) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation (eq. (11) @cite ahlrichs99a) */
  auto const interpolated_u =
      lb_lbinterpolation_get_interpolated_velocity(shifted_pos) *
      LB::get_lattice_speed();

  Utils::Vector3d v_drift = interpolated_u + vel_offset;
  /* calculate viscous force (eq. (9) @cite ahlrichs99a) */
  return -lb_lbcoupling_get_gamma() * (p.v() - v_drift);
}

template <class T, std::size_t N>
using Box = std::pair<Utils::Vector<T, N>, Utils::Vector<T, N>>;

/**
 * @brief Check if a position is in a box.
 *
 * The left boundary belong to the box, the
 * right one does not. Periodic boundaries are
 * not considered.
 *
 * @param pos Position to check
 * @param box Box to check
 *
 * @return True iff the point is inside of the box.
 */
template <class T, std::size_t N>
bool in_box(Utils::Vector<T, N> const &pos, Box<T, N> const &box) {
  return (pos >= box.first) and (pos < box.second);
}

/**
 * @brief Check if a position is within the local box + halo.
 *
 * @param pos Position to check
 * @param halo Halo
 *
 * @return True iff the point is inside of the box up to halo.
 */
inline bool in_local_domain(Utils::Vector3d const &pos, double halo = 0.) {
  auto const halo_vec = Utils::Vector3d::broadcast(halo);

  return in_box(
      pos, {local_geo.my_left() - halo_vec, local_geo.my_right() + halo_vec});
}

bool in_local_halo(Utils::Vector3d const &pos) {
  auto const halo = 0.5 * LB::get_agrid();

  return in_local_domain(pos, halo);
}

/** @brief Return a vector of positions shifted by +,- box length in each
 ** coordinate */
std::vector<Utils::Vector3d> positions_in_halo(Utils::Vector3d pos,
                                               const BoxGeometry &box) {
  std::vector<Utils::Vector3d> res;
  for (int i : {-1, 0, 1}) {
    for (int j : {-1, 0, 1}) {
      for (int k : {-1, 0, 1}) {
        Utils::Vector3d shift{{double(i), double(j), double(k)}};
        Utils::Vector3d pos_shifted =
            pos + Utils::hadamard_product(box.length(), shift);

        if (box_geo.type() == BoxType::LEES_EDWARDS) {
          auto le = box_geo.lees_edwards_bc();
          auto normal_shift = (pos_shifted - pos)[le.shear_plane_normal];
          if (normal_shift > std::numeric_limits<double>::epsilon())
            pos_shifted[le.shear_direction] += le.pos_offset;
          if (normal_shift < -std::numeric_limits<double>::epsilon())
            pos_shifted[le.shear_direction] -= le.pos_offset;
        }

        if (in_local_halo(pos_shifted)) {
          res.push_back(pos_shifted);
        }
      }
    }
  }
  return res;
}

/** @brief Return if locally there exists a physical particle
 ** for a given (ghost) particle */
bool is_ghost_for_local_particle(const Particle &p) {
  return !cell_structure.get_local_particle(p.id())->is_ghost();
}

/** @brief Determine if a given particle should be coupled.
 ** In certain cases, there may be more than one ghost for the same particle.
 ** To make sure, that these are only coupled once, ghosts' ids are stored
 ** in an unordered_set. */
bool should_be_coupled(const Particle &p,
                       std::unordered_set<int> &coupled_ghost_particles) {
  // always couple physical particles
  if (not p.is_ghost()) {
    return true;
  }
  // for ghosts check that we don't have the physical particle and
  // that a ghost for the same particle has not yet been coupled
  if ((not is_ghost_for_local_particle(p)) and
      (coupled_ghost_particles.find(p.id()) == coupled_ghost_particles.end())) {
    coupled_ghost_particles.insert(p.id());
    return true;
  }
  return false;
}

#ifdef ENGINE
void add_swimmer_force(Particle const &p, double time_step) {
  if (p.swimming().swimming) {
    // calculate source position
    auto const magnitude = p.swimming().dipole_length;
    auto const direction = static_cast<double>(p.swimming().push_pull);
    auto const director = p.calc_director();
    auto const source_position = p.pos() + direction * magnitude * director;
    auto const force = p.swimming().f_swim * director;

    // couple positions including shifts by one box length to add forces
    // to ghost layers
    for (auto pos : positions_in_halo(source_position, box_geo)) {
      add_md_force(pos, force, time_step);
    }
  }
}
#endif

Utils::Vector3d lb_particle_coupling_noise(bool enabled, int part_id,
                                           const OptionalCounter &rng_counter) {
  if (enabled) {
    if (rng_counter) {
      return Random::noise_uniform<RNGSalt::PARTICLES>(rng_counter->value(), 0,
                                                       part_id);
    }
    throw std::runtime_error(
        "Access to uninitialized LB particle coupling RNG counter");
  }
  return {};
}

void couple_particle(Particle &p, bool couple_virtual, double noise_amplitude,
                     const OptionalCounter &rng_counter, double time_step) {

  if (p.is_virtual() and not couple_virtual)
    return;

  // Calculate coupling force
  Utils::Vector3d coupling_force = {};
  for (auto pos : positions_in_halo(p.pos(), box_geo)) {
    if (in_local_halo(pos)) {
      auto const drag_force =
          lb_drag_force(p, pos, lb_particle_coupling_drift_vel_offset(p));
      auto const random_force =
          noise_amplitude * lb_particle_coupling_noise(noise_amplitude > 0.0,
                                                       p.id(), rng_counter);
      coupling_force = drag_force + random_force;
      break;
    }
  }

  // couple positions including shifts by one box length to add
  // forces to ghost layers
  for (auto pos : positions_in_halo(p.pos(), box_geo)) {
    if (in_local_domain(pos)) {
      /* Particle is in our LB volume, so this node
       * is responsible to adding its force */
      p.force() += coupling_force;
    }
    add_md_force(pos, coupling_force, time_step);
  }

#ifdef ENGINE
  add_swimmer_force(p, time_step);
#endif
}

void lb_lbcoupling_calc_particle_lattice_ia(bool couple_virtual,
                                            const ParticleRange &particles,
                                            const ParticleRange &more_particles,
                                            double time_step) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
  if (lattice_switch == ActiveLB::WALBERLA_LB) {
    if (lb_particle_coupling.couple_to_md) {
      bool using_regular_cell_structure =
          (local_geo.cell_structure_type() ==
           CellStructureType::CELL_STRUCTURE_REGULAR);
      if (not using_regular_cell_structure) {
        if (n_nodes > 1) {
          throw std::runtime_error("LB only works with regular decomposition, "
                                   "if using more than one MPI node");
        }
      }
      using Utils::sqr;
      auto const kT = LB::get_kT() * sqr(LB::get_lattice_speed());
      /* Eq. (16) @cite ahlrichs99a.
       * The factor 12 comes from the fact that we use random numbers
       * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
       * time_step comes from the discretization.
       */
      auto const noise_amplitude =
          (kT > 0.)
              ? std::sqrt(12. * 2. * lb_lbcoupling_get_gamma() * kT / time_step)
              : 0.0;

      std::unordered_set<int> coupled_ghost_particles;

      /* Couple particles ranges */
      for (auto &p : particles) {
        if (should_be_coupled(p, coupled_ghost_particles)) {
          couple_particle(p, couple_virtual, noise_amplitude,
                          lb_particle_coupling.rng_counter_coupling, time_step);
        }
      }

      for (auto &p : more_particles) {
        if (should_be_coupled(p, coupled_ghost_particles)) {
          couple_particle(p, couple_virtual, noise_amplitude,
                          lb_particle_coupling.rng_counter_coupling, time_step);
        }
      }
    }
  }
}

void lb_lbcoupling_propagate() {
  if (lattice_switch != ActiveLB::NONE) {
    if (LB::get_kT() > 0.0) {
      if (lattice_switch == ActiveLB::WALBERLA_LB) {
        lb_particle_coupling.rng_counter_coupling->increment();
      }
    }
  }
}

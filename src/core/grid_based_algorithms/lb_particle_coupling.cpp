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
#include "lb_particle_coupling.hpp"
#include "LocalBox.hpp"
#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "config/config.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/OptionalCounter.hpp"
#include "integrate.hpp"
#include "lb_interface.hpp"
#include "lb_interpolation.hpp"
#include "lbgpu.hpp"
#include "random.hpp"

#include <profiler/profiler.hpp>
#include <utils/Counter.hpp>
#include <utils/Vector.hpp>

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
  mpi_bcast_lb_particle_coupling();
}

double lb_lbcoupling_get_gamma() { return lb_particle_coupling.gamma; }

bool lb_lbcoupling_is_seed_required() {
  if (lattice_switch == ActiveLB::CPU) {
    return not lb_particle_coupling.rng_counter_coupling.is_initialized();
  }
#ifdef CUDA
  if (lattice_switch == ActiveLB::GPU) {
    return not rng_counter_coupling_gpu.is_initialized();
  }
#endif
  return false;
}

uint64_t lb_coupling_get_rng_state_cpu() {
  return lb_particle_coupling.rng_counter_coupling->value();
}

uint64_t lb_lbcoupling_get_rng_state() {
  if (lattice_switch == ActiveLB::CPU) {
    return lb_coupling_get_rng_state_cpu();
  }
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    return lb_coupling_get_rng_state_gpu();
#endif
  }
  return {};
}

void lb_lbcoupling_set_rng_state(uint64_t counter) {
  if (lattice_switch == ActiveLB::CPU) {
    lb_particle_coupling.rng_counter_coupling =
        Utils::Counter<uint64_t>(counter);
    mpi_bcast_lb_particle_coupling();
  } else if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    lb_coupling_set_rng_state_gpu(counter);
#endif
  }
}

namespace {
/**
 * @brief Add a force to the lattice force density.
 * @param pos Position of the force
 * @param force Force in MD units.
 * @param time_step MD time step.
 */
void add_md_force(Utils::Vector3d const &pos, Utils::Vector3d const &force,
                  double time_step) {
  /* transform momentum transfer to lattice units
     (eq. (12) @cite ahlrichs99a) */
  auto const delta_j = -(time_step / lb_lbfluid_get_lattice_speed()) * force;
  lb_lbinterpolation_add_force_density(pos, delta_j);
}
} // namespace

/** Coupling of a single particle to viscous fluid with Stokesian friction.
 *
 *  Section II.C. @cite ahlrichs99a
 *
 *  @param[in] p             The coupled particle.
 *  @param[in] pos           Local position of particle or its ghost.
 *  @param[in] f_random      Additional force to be included.
 *
 *  @return The viscous coupling force plus @p f_random.
 */
Utils::Vector3d lb_viscous_coupling(Particle const &p,
                                    Utils::Vector3d const &pos,
                                    Utils::Vector3d const &f_random) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation (eq. (11) @cite ahlrichs99a) */
  auto const interpolated_u =
      lb_lbinterpolation_get_interpolated_velocity(pos) *
      lb_lbfluid_get_lattice_speed();

  Utils::Vector3d v_drift = interpolated_u;
#ifdef ENGINE
  if (p.swimming().swimming) {
    v_drift += p.swimming().v_swim * p.calc_director();
  }
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  v_drift += p.mu_E();
#endif

  /* calculate viscous force (eq. (9) @cite ahlrichs99a) */
  auto const force = -lb_lbcoupling_get_gamma() * (p.v() - v_drift) + f_random;

  return force;
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

/**
 * @brief Check if a position is within the local LB domain
 *       plus halo.
 *
 * @param pos Position to check
 *
 * @return True iff the point is inside of the domain.
 */
bool in_local_halo(Utils::Vector3d const &pos) {
  auto const halo = 0.5 * lb_lbfluid_get_agrid();

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
    const double direction =
        double(p.swimming().push_pull) * p.swimming().dipole_length;
    auto const director = p.calc_director();
    auto const source_position = p.pos() + direction * director;
    auto const force = p.swimming().f_swim * director;

    // couple positions including shifts by one box length to add forces
    // to ghost layers
    for (auto pos : positions_in_halo(source_position, box_geo)) {
      add_md_force(pos, force, time_step);
    }
  }
}
#endif

void lb_lbcoupling_calc_particle_lattice_ia(bool couple_virtual,
                                            const ParticleRange &particles,
                                            const ParticleRange &more_particles,
                                            double time_step) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    if (lb_particle_coupling.couple_to_md && this_node == 0) {
      switch (lb_lbinterpolation_get_interpolation_order()) {
      case (InterpolationOrder::linear):
        lb_calc_particle_lattice_ia_gpu<8>(
            couple_virtual, lb_lbcoupling_get_gamma(), time_step);
        break;
      case (InterpolationOrder::quadratic):
        lb_calc_particle_lattice_ia_gpu<27>(
            couple_virtual, lb_lbcoupling_get_gamma(), time_step);
        break;
      }
    }
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
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
      switch (lb_lbinterpolation_get_interpolation_order()) {
      case (InterpolationOrder::quadratic):
        throw std::runtime_error("The non-linear interpolation scheme is not "
                                 "implemented for the CPU LB.");
      case (InterpolationOrder::linear): {
        auto const kT = lb_lbfluid_get_kT();
        /* Eq. (16) @cite ahlrichs99a.
         * The factor 12 comes from the fact that we use random numbers
         * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
         * time_step comes from the discretization.
         */
        auto const noise_amplitude =
            (kT > 0.) ? std::sqrt(12. * 2. * lb_lbcoupling_get_gamma() * kT /
                                  time_step)
                      : 0.0;

        auto f_random = [noise_amplitude](int id) -> Utils::Vector3d {
          if (noise_amplitude > 0.0) {
            return Random::noise_uniform<RNGSalt::PARTICLES>(
                lb_particle_coupling.rng_counter_coupling->value(), 0, id);
          }
          return {};
        };

        auto couple_particle = [&](Particle &p) -> void {
          if (p.is_virtual() and !couple_virtual)
            return;

          // Calculate coupling force
          Utils::Vector3d force = {};
          for (auto pos : positions_in_halo(p.pos(), box_geo)) {
            if (in_local_halo(pos)) {
              force = lb_viscous_coupling(p, pos,
                                          noise_amplitude * f_random(p.id()));
              break;
            }
          }

          // couple positions including shifts by one box length to add
          // forces to ghost layers
          for (auto pos : positions_in_halo(p.pos(), box_geo)) {
            if (in_local_domain(pos)) {
              /* if the particle is in our LB volume, this node
               * is responsible to adding its force */
              p.force() += force;
            }
            add_md_force(pos, force, time_step);
          }

#ifdef ENGINE
          add_swimmer_force(p, time_step);
#endif
        };

        std::unordered_set<int> coupled_ghost_particles;

        /* Couple particles ranges */
        for (auto &p : particles) {
          if (should_be_coupled(p, coupled_ghost_particles)) {
            couple_particle(p);
          }
        }

        for (auto &p : more_particles) {
          if (should_be_coupled(p, coupled_ghost_particles)) {
            couple_particle(p);
          }
        }

        break;
      }
      }
    }
  }
}

void lb_lbcoupling_propagate() {
  if (lattice_switch != ActiveLB::NONE) {
    if (lb_lbfluid_get_kT() > 0.0) {
      if (lattice_switch == ActiveLB::CPU) {
        lb_particle_coupling.rng_counter_coupling->increment();
      } else if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
        rng_counter_coupling_gpu->increment();
#endif
      }
    }
  }
}

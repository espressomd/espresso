/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "cells.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "errorhandling.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lattice.hpp"
#include "integrate.hpp"
#include "lb_interface.hpp"
#include "lb_interpolation.hpp"
#include "lbgpu.hpp"
#include "particle_data.hpp"
#include "random.hpp"

#include <profiler/profiler.hpp>
#include <utils/Counter.hpp>
#include <utils/u32_to_u64.hpp>
#include <utils/uniform.hpp>

#include <Random123/philox.h>
#include <boost/mpi.hpp>

LB_Particle_Coupling lb_particle_coupling;

void mpi_bcast_lb_particle_coupling_slave() {
  boost::mpi::broadcast(comm_cart, lb_particle_coupling, 0);
}

REGISTER_CALLBACK(mpi_bcast_lb_particle_coupling_slave)

void mpi_bcast_lb_particle_coupling() {
  mpi_call(mpi_bcast_lb_particle_coupling_slave);
  boost::mpi::broadcast(comm_cart, lb_particle_coupling, 0);
}

void lb_lbcoupling_activate() { lb_particle_coupling.couple_to_md = true; }

void lb_lbcoupling_deactivate() {
  if (lattice_switch != ActiveLB::NONE && this_node == 0 &&
      lb_particle_coupling.gamma > 0.) {
    runtimeWarningMsg()
        << "Recalculating forces, so the LB coupling forces are not "
           "included in the particle force the first time step. This "
           "only matters if it happens frequently during "
           "sampling.";
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
 */
void add_md_force(Utils::Vector3d const &pos, Utils::Vector3d const &force) {
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
 *  @param[in]     f_random  Additional force to be included.
 *
 *  @return The viscous coupling force plus f_random.
 */
Utils::Vector3d lb_viscous_coupling(Particle const &p,
                                    Utils::Vector3d const &f_random) {
  /* calculate fluid velocity at particle's position
     this is done by linear interpolation (eq. (11) @cite ahlrichs99a) */
  auto const interpolated_u =
      lb_lbinterpolation_get_interpolated_velocity(p.r.p) *
      lb_lbfluid_get_lattice_speed();

  Utils::Vector3d v_drift = interpolated_u;
#ifdef ENGINE
  if (p.p.swim.swimming) {
    v_drift += p.p.swim.v_swim * p.r.calc_director();
  }
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  v_drift += p.p.mu_E;
#endif

  /* calculate viscous force (eq. (9) @cite ahlrichs99a) */
  auto const force = -lb_lbcoupling_get_gamma() * (p.m.v - v_drift) + f_random;

  add_md_force(p.r.p, force);

  return force;
}

namespace {
using Utils::Vector;
using Utils::Vector3d;

template <class T, size_t N> using Box = std::pair<Vector<T, N>, Vector<T, N>>;

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
template <class T, size_t N>
bool in_box(Vector<T, N> const &pos, Box<T, N> const &box) {
  return (pos >= box.first) and (pos < box.second);
}

/**
 * @brief Check if a position is within the local box + halo.
 *
 * @param pos Position to check
 * @param local_box Geometry to check
 * @param halo Halo
 *
 * @return True iff the point is inside of the box up to halo.
 */
template <class T>
bool in_local_domain(Vector<T, 3> const &pos, LocalBox<T> const &local_box,
                     T const &halo = {}) {
  auto const halo_vec = Vector<T, 3>::broadcast(halo);

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
bool in_local_halo(Vector3d const &pos) {
  auto const halo = 0.5 * lb_lbfluid_get_agrid();

  return in_local_domain(pos, local_geo, halo);
}

#ifdef ENGINE
void add_swimmer_force(Particle &p) {
  if (p.p.swim.swimming) {
    // calculate source position
    const double direction =
        double(p.p.swim.push_pull) * p.p.swim.dipole_length;
    auto const director = p.r.calc_director();
    auto const source_position = p.r.p + direction * director;

    if (not in_local_halo(source_position)) {
      return;
    }

    add_md_force(source_position, p.p.swim.f_swim * director);
  }
}
#endif

} // namespace

void lb_lbcoupling_calc_particle_lattice_ia(
    bool couple_virtual, const ParticleRange &particles,
    const ParticleRange &more_particles) {
  ESPRESSO_PROFILER_CXX_MARK_FUNCTION;
  if (lattice_switch == ActiveLB::GPU) {
#ifdef CUDA
    if (lb_particle_coupling.couple_to_md && this_node == 0) {
      switch (lb_lbinterpolation_get_interpolation_order()) {
      case (InterpolationOrder::linear):
        lb_calc_particle_lattice_ia_gpu<8>(couple_virtual,
                                           lb_lbcoupling_get_gamma());
        break;
      case (InterpolationOrder::quadratic):
        lb_calc_particle_lattice_ia_gpu<27>(couple_virtual,
                                            lb_lbcoupling_get_gamma());
        break;
      }
    }
#endif
  } else if (lattice_switch == ActiveLB::CPU) {
    if (lb_particle_coupling.couple_to_md) {
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
                lb_particle_coupling.rng_counter_coupling->value(), id);
          }
          return {};
        };

        auto couple_particle = [&](Particle &p) -> void {
          if (p.p.is_virtual and !couple_virtual)
            return;

          /* Particle is in our LB volume, so this node
           * is responsible to adding its force */
          if (in_local_domain(p.r.p, local_geo)) {
            auto const force = lb_viscous_coupling(
                p, noise_amplitude * f_random(p.identity()));
            /* add force to the particle */
            p.f.f += force;
            /* Particle is not in our domain, but adds to the force
             * density in our domain, only calculate contribution to
             * the LB force density. */
          } else if (in_local_halo(p.r.p)) {
            lb_viscous_coupling(p, noise_amplitude * f_random(p.identity()));
          }

#ifdef ENGINE
          add_swimmer_force(p);
#endif
        };

        /* Couple particles ranges */
        for (auto &p : particles) {
          couple_particle(p);
        }

        for (auto &p : more_particles) {
          couple_particle(p);
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

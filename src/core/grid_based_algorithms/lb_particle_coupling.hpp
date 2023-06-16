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
#ifndef LB_PARTICLE_COUPLING_HPP
#define LB_PARTICLE_COUPLING_HPP

#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "grid.hpp"

#include <utils/Counter.hpp>
#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <boost/optional.hpp>
#include <boost/serialization/access.hpp>
#include <boost/serialization/optional.hpp>

#include <cassert>
#include <cmath>
#include <cstdint>
#include <unordered_set>
#include <vector>

using OptionalCounter = boost::optional<Utils::Counter<uint64_t>>;

void lb_lbcoupling_propagate();
uint64_t lb_lbcoupling_get_rng_state();
void lb_lbcoupling_set_rng_state(uint64_t counter);
void lb_lbcoupling_set_gamma(double friction);
double lb_lbcoupling_get_gamma();
bool lb_lbcoupling_is_seed_required();

/**
 * @brief Activate the coupling between LB and MD particles.
 * @note This is a collective function and needs to be called from all
 * processes.
 */
void lb_lbcoupling_activate();

/**
 * @brief Deactivate the coupling between LB and MD particles.
 * @note This is a collective function and needs to be called from all
 * processes.
 */
void lb_lbcoupling_deactivate();

/**
 * @brief Check if a position is within the local LB domain plus halo.
 *
 * @param pos Position to check
 *
 * @return True iff the point is inside of the domain.
 */
bool in_local_halo(Utils::Vector3d const &pos);

/**
 * @brief Add a force to the lattice force density.
 * @param pos Position of the force
 * @param force Force in MD units.
 * @param time_step MD time step.
 */
void add_md_force(Utils::Vector3d const &pos, Utils::Vector3d const &force,
                  double time_step);

// internal function exposed for unit testing
std::vector<Utils::Vector3d> positions_in_halo(Utils::Vector3d pos,
                                               const BoxGeometry &box);

// internal function exposed for unit testing
void add_swimmer_force(Particle const &p, double time_step);

/**
 * @brief Calculate particle drift velocity offset due to ENGINE and
 * ELECTROHYDRODYNAMICS.
 */
Utils::Vector3d lb_particle_coupling_drift_vel_offset(const Particle &p);

void mpi_bcast_lb_particle_coupling();

/** @brief Calculate drag force on a single particle.
 *
 *  See section II.C. @cite ahlrichs99a
 *
 *  @param[in] p           The coupled particle
 *  @param[in] shifted_pos The particle position with optional shift
 *  @param[in] vel_offset  Velocity offset to be added to interpolated LB
 *                         velocity before calculating the force
 *
 *  @return The viscous coupling force
 */
Utils::Vector3d lb_drag_force(Particle const &p,
                              Utils::Vector3d const &shifted_pos,
                              Utils::Vector3d const &vel_offset);

namespace LB {
struct ParticleCouplingConfig {
  OptionalCounter rng_counter_coupling = {};
  /** @brief Friction coefficient for the particle coupling. */
  double gamma = 0.0;
  bool couple_to_md = false;

private:
  friend class boost::serialization::access;

  template <class Archive> void serialize(Archive &ar, const unsigned int) {
    ar &rng_counter_coupling;
    ar &gamma;
    ar &couple_to_md;
  }
};
} // namespace LB

// internal global exposed for unit testing
extern LB::ParticleCouplingConfig lb_particle_coupling;

namespace LB {

/** @brief Calculate particle-lattice interactions. */
void couple_particles(bool couple_virtual, ParticleRange const &real_particles,
                      ParticleRange const &ghost_particles, double time_step);

class ParticleCoupling {
  bool m_couple_virtual;
  bool m_thermalized;
  double m_time_step;
  double m_noise_pref_wo_gamma;

public:
  ParticleCoupling(bool couple_virtual, double time_step, double kT)
      : m_couple_virtual{couple_virtual}, m_thermalized{kT != 0.},
        m_time_step{time_step} {
    assert(kT >= 0.);
    /* Eq. (16) @cite ahlrichs99a, without the gamma term.
     * The factor 12 comes from the fact that we use random numbers
     * from -0.5 to 0.5 (equally distributed) which have variance 1/12.
     * The time step comes from the discretization.
     */
    auto constexpr variance_inv = 12.;
    m_noise_pref_wo_gamma = std::sqrt(variance_inv * 2. * kT / time_step);
  }

  ParticleCoupling(bool couple_virtual, double time_step)
      : ParticleCoupling(couple_virtual, time_step,
                         LB::get_kT() * Utils::sqr(LB::get_lattice_speed())) {}

  Utils::Vector3d get_noise_term(Particle const &p) const;
  void kernel(Particle &p);
};

/**
 * @brief Keep track of ghost particles that have already been coupled once.
 * In certain cases, there may be more than one ghost for the same particle.
 * To make sure these are only coupled once, ghosts' ids are recorded.
 */
class CouplingBookkeeping {
  std::unordered_set<int> m_coupled_ghosts;

  /** @brief Check if there is locally a real particle for the given ghost. */
  bool is_ghost_for_local_particle(Particle const &p) const;

public:
  /** @brief Determine if a given particle should be coupled. */
  bool should_be_coupled(Particle const &p) {
    // real particles: always couple
    if (not p.is_ghost()) {
      return true;
    }
    // ghosts: check we don't have the corresponding real particle on the same
    // node, and that a ghost for the same particle hasn't been coupled already
    if (m_coupled_ghosts.count(p.id()) != 0 or is_ghost_for_local_particle(p)) {
      return false;
    }
    m_coupled_ghosts.insert(p.id());
    return true;
  }
};

} // namespace LB

#endif

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
#ifndef LB_PARTICLE_COUPLING_HPP
#define LB_PARTICLE_COUPLING_HPP

#include "BoxGeometry.hpp"
#include "OptionalCounter.hpp"
#include "ParticleRange.hpp"

#include <boost/serialization/access.hpp>

#include <cstdint>
#include <unordered_set>

/** Calculate particle lattice interactions.
 *  So far, only viscous coupling with Stokesian friction is implemented.
 *  Include all particle-lattice forces in this function.
 *  The function is called from \ref force_calc.
 */
void lb_lbcoupling_calc_particle_lattice_ia(bool couple_virtual,
                                            const ParticleRange &particles,
                                            const ParticleRange &more_particles,
                                            double time_step);
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

struct LB_Particle_Coupling {
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

// expose functions that are also used to couple lb_inertialess_tracers
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

bool in_local_halo(Utils::Vector3d const &pos);
std::vector<Utils::Vector3d> positions_in_halo(Utils::Vector3d pos,
                                               const BoxGeometry &box);
bool is_ghost_for_local_particle(const Particle &p);
bool should_be_coupled(const Particle &p,
                       std::unordered_set<int> &coupled_ghost_particles);

#endif

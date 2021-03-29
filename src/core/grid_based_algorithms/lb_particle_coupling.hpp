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

#include "OptionalCounter.hpp"
#include "ParticleRange.hpp"

#include <utils/Counter.hpp>
#include <utils/Vector.hpp>

#include <boost/serialization/access.hpp>

#include <cstdint>

/** Calculate particle lattice interactions.
 *  So far, only viscous coupling with Stokesian friction is implemented.
 *  Include all particle-lattice forces in this function.
 *  The function is called from \ref force_calc.
 */
void lb_lbcoupling_calc_particle_lattice_ia(
    bool couple_virtual, const ParticleRange &particles,
    const ParticleRange &more_particles);
void lb_lbcoupling_propagate();
uint64_t lb_lbcoupling_get_rng_state();
void lb_lbcoupling_set_rng_state(uint64_t counter);
void lb_lbcoupling_set_gamma(double friction);
double lb_lbcoupling_get_gamma();
bool lb_lbcoupling_is_seed_required();

/*
 * @brief Activate the coupling between LB and MD particles.
 * @note This is a collective function and needs to be called from all
 * processes.
 */
void lb_lbcoupling_activate();

/*
 * @brief Deactivate the coupling between LB and MD particles.
 * @note This is a collective function and needs to be called from all
 * processes.
 */
void lb_lbcoupling_deactivate();

/**
 * @brief Check if a position is within the local LB domain
 *       plus halo.
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
 */
void add_md_force(Utils::Vector3d const &pos, Utils::Vector3d const &force);

Utils::Vector3d lb_particle_coupling_noise(bool enabled, int part_id,
                                           const OptionalCounter &rng_counter);

Utils::Vector3d lb_particle_coupling_drift_vel_offset(const Particle &p);

void mpi_bcast_lb_particle_coupling();

/** calculate drag force on a single particle
 *
 *  Section II.C. @cite ahlrichs99a
 *
 *  @param[in] p             The coupled particle.
 *  @param vel_offset        Velocity offset to be added to interpolated LB
 * velocity before calculating the force
 *
 *  @return The viscous coupling force
 */
Utils::Vector3d lb_drag_force(Particle const &p,
                              const Utils::Vector3d &vel_offset);

struct LB_Particle_Coupling {
  OptionalCounter rng_counter_coupling;
  /*
   * @brief Friction constant for the particle coupling.
   */
  double gamma = 0.0;
  bool couple_to_md = false;

private:
  friend class boost::serialization::access;

  template <class Archive>
  void serialize(Archive &ar, const unsigned int version) {
    ar &rng_counter_coupling;
    ar &gamma;
    ar &couple_to_md;
  }
};

extern LB_Particle_Coupling lb_particle_coupling;

#endif

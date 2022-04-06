/*
 * Copyright (C) 2010-2019 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
#ifndef CORE_PARTICLE_NODE_HPP
#define CORE_PARTICLE_NODE_HPP
/** \file
 *  Particles creation and deletion.
 *
 *  This file contains everything related to particle storage and tracking.
 *
 *  Implementation in particle_node.cpp.
 */

#include "Particle.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <cstddef>
#include <vector>

/**
 * @brief Get particle data.
 *
 *  @param part the identity of the particle to fetch
 */
const Particle &get_particle_data(int p_id);

/**
 * @brief Fetch a range of particle into the fetch cache.
 *
 *
 * If the range is larger than the cache size, only
 * the particle that fit into the cache are fetched.
 *
 * The particles have to exist, an exception it throw
 * if one of the the particles can not be found.
 *
 * @param ids Ids of the particles that should be fetched.
 */
void prefetch_particle_data(Utils::Span<const int> ids);

/** @brief Invalidate the fetch cache for get_particle_data. */
void invalidate_fetch_cache();

/** @brief Return the maximal number of particles that are
 *         kept in the fetch cache.
 */
std::size_t fetch_cache_max_size();

/** Invalidate \ref particle_node. This has to be done
 *  at the beginning of the integration.
 */
void clear_particle_node();

/** Call only on the head node.
 *  Move a particle to a new position.
 *  If it does not exist, it is created.
 *  @param p_id     identity of the particle to move (or create)
 *  @param pos      position
 */
void place_particle(int p_id, Utils::Vector3d const &pos);

/** Remove particle with a given identity. Also removes all bonds to the
 *  particle.
 *  @param p_id     identity of the particle to remove
 */
void remove_particle(int p_id);

/** Remove all particles. */
void remove_all_particles();

void init_type_map(int type);
void on_particle_type_change(int p_id, int type);

/** Find a particle of given type and return its id */
int get_random_p_id(int type, int random_index_in_type_map);
int number_of_particles_with_type(int type);

/**
 * @brief Check if particle exists.
 *
 * @param p_id     identity of the particle
 * @return True iff the particle exists.
 */
bool particle_exists(int p_id);

/**
 *  @brief Get the MPI rank which owns the a specific particle.
 *
 *  @param p_id     identity of the particle
 *  @return The MPI rank the particle is on.
 */
int get_particle_node(int p_id);

/**
 * @brief Get all particle ids.
 *
 * @return Sorted ids of all existing particles.
 */
std::vector<int> get_particle_ids();

/**
 * @brief Get maximal particle id.
 */
int get_maximal_particle_id();

/**
 * @brief Get number of particles.
 */
int get_n_part();

#endif

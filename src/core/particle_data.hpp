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
#ifndef _PARTICLE_DATA_H
#define _PARTICLE_DATA_H
/** \file
 *  Particles and particle lists.
 *
 *  This file contains everything related to particle storage. If you want to
 *  add a new property to the particles, it is probably a good idea to modify
 *  Particle to give scripts access to that property. You always have to modify
 *  two positions: first the print section, where you should add your new
 *  data at the end, and second the read section where you have to find a nice
 *  and short name for your property to appear in the Python code. Then you
 *  just parse your part out of argc and argv.
 *
 *  Implementation in particle_data.cpp.
 */

#include "config.hpp"

#include "Particle.hpp"
#include "ParticleList.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>

#include <memory>

/************************************************
 * defines
 ************************************************/

enum {
  /// ok code for \ref place_particle
  ES_PART_OK = 0,
  /// error code for \ref place_particle
  ES_PART_ERROR = -1,
  /// ok code for \ref place_particle, particle is new
  ES_PART_CREATED = 1
};

#ifdef EXTERNAL_FORCES
/**
 *  \ref ParticleProperties::ext_flag "ext_flag" value for fixed coordinate
 *  coord.
 */
#define COORD_FIXED(coord) (2u << (coord))
/** \ref ParticleProperties::ext_flag "ext_flag" mask to check whether any of
 *  the coordinates is fixed. */
#define COORDS_FIX_MASK (COORD_FIXED(0) | COORD_FIXED(1) | COORD_FIXED(2))
#else
#define COORD_FIXED(coord) (0)
#endif

/************************************************
 * Functions
 ************************************************/

/*    Other Functions                           */
/************************************************/

/** Invalidate \ref particle_node. This has to be done
 *  at the beginning of the integration.
 */
void clear_particle_node();

/**
 * @brief Get particle data.
 *
 *   @param part the identity of the particle to fetch
 *   @return Pointer to copy of particle if it exists,
 *           nullptr otherwise;
 */
const Particle &get_particle_data(int part);

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
size_t fetch_cache_max_size();

/** Call only on the master node.
 *  Move a particle to a new position.
 *  If it does not exist, it is created.
 *  @param part the identity of the particle to move
 *  @param p    its new position
 *  @retval ES_PART_OK if particle existed
 *  @retval ES_PART_CREATED if created
 *  @retval ES_PART_ERROR if id is illegal
 */
int place_particle(int part, const double *p);

/** Call only on the master node: set particle velocity.
 *  @param part the particle.
 *  @param v its new velocity.
 */
void set_particle_v(int part, double *v);

#ifdef ENGINE
/** Call only on the master node: set particle velocity.
 *  @param part the particle.
 *  @param swim struct containing swimming parameters
 */
void set_particle_swimming(int part, ParticleParametersSwimming swim);
#endif

/** Call only on the master node: set particle force.
 *  @param part the particle.
 *  @param F its new force.
 */
void set_particle_f(int part, const Utils::Vector3d &F);

/** Call only on the master node: set particle mass.
 *  @param part the particle.
 *  @param mass its new mass.
 */
void set_particle_mass(int part, double mass);

#ifdef ROTATIONAL_INERTIA
/** Call only on the master node: set particle rotational inertia.
 *  @param part the particle.
 *  @param rinertia its new inertia.
 */
void set_particle_rotational_inertia(int part, double *rinertia);
#endif

/** Call only on the master node: Specifies whether a particle's rotational
 *  degrees of freedom are integrated or not. If set to zero, the content of
 *  the torque and omega variables are meaningless
 *  @param part the particle.
 *  @param rot the degrees of freedom flag.
 */
void set_particle_rotation(int part, int rot);

/** @brief rotate a particle around an axis
 *
 *  @param part particle id
 *  @param axis rotation axis
 *  @param angle rotation angle
 */
void rotate_particle(int part, const Utils::Vector3d &axis, double angle);

/** Call only on the master node: set particle charge.
 *  @param part the particle.
 *  @param q its new charge.
 */
void set_particle_q(int part, double q);

#ifdef LB_ELECTROHYDRODYNAMICS
/** Call only on the master node: set particle electrophoretic mobility.
 *  @param part the particle.
 *  @param mu_E its new mobility.
 */
void set_particle_mu_E(int part, Utils::Vector3d const &mu_E);
void get_particle_mu_E(int part, Utils::Vector3d &mu_E);
#endif

/** Call only on the master node: set particle type.
 *  @param p_id the particle.
 *  @param type its new type.
 */
void set_particle_type(int p_id, int type);

/** Call only on the master node: set particle's molecule id.
 *  @param part the particle.
 *  @param mid  its new mol id.
 */
void set_particle_mol_id(int part, int mid);

#ifdef ROTATION
/** Call only on the master node: set particle orientation using quaternions.
 *  @param part the particle.
 *  @param quat its new value for quaternions.
 */
void set_particle_quat(int part, double *quat);

/** Call only on the master node: set particle angular velocity from lab frame.
 *  @param part the particle.
 *  @param omega_lab its new angular velocity.
 */
void set_particle_omega_lab(int part, const Utils::Vector3d &omega_lab);

/** Call only on the master node: set particle angular velocity in body frame.
 *  @param part the particle.
 *  @param omega its new angular velocity.
 */
void set_particle_omega_body(int part, const Utils::Vector3d &omega);

/** Call only on the master node: set particle torque from lab frame.
 *  @param part the particle.
 *  @param torque_lab its new torque.
 */
void set_particle_torque_lab(int part, const Utils::Vector3d &torque_lab);

#endif

#ifdef DIPOLES
/** Call only on the master node: set particle dipole orientation.
 *  @param part the particle.
 *  @param dip its new dipole orientation.
 */
void set_particle_dip(int part, double const *dip);

/** Call only on the master node: set particle dipole moment (absolute value).
 *  @param part the particle.
 *  @param dipm its new dipole moment.
 */
void set_particle_dipm(int part, double dipm);
#endif

#ifdef VIRTUAL_SITES
/** Call only on the master node: set particle virtual flag.
 *  @param part the particle.
 *  @param is_virtual new @ref ParticleProperties::is_virtual "is_virtual" flag.
 */
void set_particle_virtual(int part, bool is_virtual);
#endif
#ifdef VIRTUAL_SITES_RELATIVE
void set_particle_vs_quat(int part, Utils::Vector4d const &vs_relative_quat);
void set_particle_vs_relative(int part, int vs_relative_to, double vs_distance,
                              Utils::Vector4d const &rel_ori);
#endif

#ifdef LANGEVIN_PER_PARTICLE
/** Call only on the master node: set particle temperature.
 *  @param part the particle.
 *  @param T its new temperature.
 */
void set_particle_temperature(int part, double T);

/** Call only on the master node: set particle frictional coefficient.
 *  @param part the particle.
 *  @param gamma its new frictional coefficient.
 */
#ifndef PARTICLE_ANISOTROPY
void set_particle_gamma(int part, double gamma);
#else
void set_particle_gamma(int part, Utils::Vector3d gamma);
#endif
#ifdef ROTATION
#ifndef PARTICLE_ANISOTROPY
void set_particle_gamma_rot(int part, double gamma);
#else
void set_particle_gamma_rot(int part, Utils::Vector3d gamma_rot);
#endif
#endif
#endif // LANGEVIN_PER_PARTICLE

#ifdef EXTERNAL_FORCES
#ifdef ROTATION
/** Call only on the master node: set particle external torque.
 *  @param part  the particle.
 *  @param torque new value for ext_torque.
 */
void set_particle_ext_torque(int part, const Utils::Vector3d &torque);
#endif
/** Call only on the master node: set particle external force.
 *  @param part  the particle.
 *  @param force new value for ext_force.
 */
void set_particle_ext_force(int part, const Utils::Vector3d &force);
/** Call only on the master node: set coordinate axes for which the particles
 *  motion is fixed.
 *  @param part  the particle.
 *  @param flag new value for flagged coordinate axes to be fixed
 */
void set_particle_fix(int part, uint8_t flag);
#endif

/** Call only on the master node: remove bond from particle.
 *  @param part     identity of principal atom of the bond.
 *  @param bond     field containing the bond type number and the identity
 *                  of all bond partners (secondary atoms of the bond).
 */
void delete_particle_bond(int part, Utils::Span<const int> bond);

/** Call only on the master node: remove all bonds from particle.
 *  @param part     identity of principal atom of the bond.
 */
void delete_particle_bonds(int part);

/** Call only on the master node: Add bond to particle.
 *  @param part     identity of principal atom of the bond.
 *  @param bond     field containing the bond type number and the
 *  identity of all bond partners (secondary atoms of the bond).
 */
void add_particle_bond(int part, Utils::Span<const int> bond);

const std::vector<BondView> &get_particle_bonds(int part);

#ifdef EXCLUSIONS
/** Call only on the master node: change particle constraints.
 *  @param part     identity of particle for which the exclusion is set.
 *  @param part2    identity of particle for which the exclusion is set.
 *                  If -1, delete all exclusions.
 *  @param _delete  if true, do not add the exclusion, rather delete it if
 *                  found
 *  @retval ES_OK on success
 *  @retval ES_ERROR on failure (e.g. particles do not exist / did not have
 *          exclusion set)
 */
int change_exclusion(int part, int part2, int _delete);

/** remove all exclusions. */
void remove_all_exclusions();
#endif

/** Remove particle with a given identity. Also removes all bonds to the
 *  particle.
 *  @param part     identity of the particle to remove
 *  @retval ES_OK on success
 *  @retval ES_ERROR if particle does not exist
 */
int remove_particle(int part);

/** Remove all particles. */
void remove_all_particles();

/** Used by \ref mpi_place_particle, should not be used elsewhere.
 *  Move a particle to a new position. If it does not exist, it is created.
 *  The position must be on the local node!
 *
 *  @param id    the identity of the particle to move
 *  @param pos   its new position
 *  @param _new  if true, the particle is allocated, else has to exists already
 *
 *  @return Pointer to the particle.
 */
Particle *local_place_particle(int id, const Utils::Vector3d &pos, int _new);

/** Used for example by \ref mpi_send_exclusion.
 *  Locally add an exclusion to a particle.
 *  @param part1 the identity of the first exclusion partner
 *  @param part2 the identity of the second exclusion partner
 *  @param _delete if true, delete the exclusion instead of add
 */
void local_change_exclusion(int part1, int part2, int _delete);

/** Used by \ref mpi_rescale_particles, should not be used elsewhere.
 *  Locally rescale all particles on current node.
 *  @param dir   direction to scale (0/1/2 = x/y/z, 3 = x+y+z isotropically)
 *  @param scale factor by which to rescale (>1: stretch, <1: contract)
 */
void local_rescale_particles(int dir, double scale);

/** Automatically add the next \<distance\> neighbors in each molecule to the
 *  exclusion list.
 *  This uses the bond topology obtained directly from the particles.
 *  To easily setup the bonds, all data should be on a single node,
 *  therefore the \ref partCfg array is used. With large amounts of
 *  particles, you should avoid this function and setup exclusions manually.
 */
void auto_exclusions(int distance);

void init_type_map(int type);

/* find a particle of given type and return its id */
int get_random_p_id(int type, int random_index_in_type_map);
int number_of_particles_with_type(int type);

// The following functions are used by the python interface to obtain
// properties of a particle, which are only compiled in in some configurations
// This is needed, because cython does not support conditional compilation
// within a ctypedef definition

#ifdef ROTATION
void pointer_to_omega_body(Particle const *p, double const *&res);

inline Utils::Vector3d get_torque_body(const Particle &p) { return p.f.torque; }

void pointer_to_quat(Particle const *p, double const *&res);

#endif

void pointer_to_q(Particle const *p, double const *&res);

#ifdef VIRTUAL_SITES
void pointer_to_virtual(Particle const *p, bool const *&res);
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void pointer_to_vs_quat(Particle const *p, double const *&res);
void pointer_to_vs_relative(Particle const *p, int const *&res1,
                            double const *&res2, double const *&res3);
#endif

void pointer_to_dipm(Particle const *P, double const *&res);

#ifdef EXTERNAL_FORCES
void pointer_to_ext_force(Particle const *p, double const *&res2);
#ifdef ROTATION
void pointer_to_ext_torque(Particle const *p, double const *&res2);
#endif
void pointer_to_fix(Particle const *p, const uint8_t *&res);
#endif

#ifdef LANGEVIN_PER_PARTICLE
void pointer_to_gamma(Particle const *p, double const *&res);
void pointer_to_temperature(Particle const *p, double const *&res);
#ifdef ROTATION
void pointer_to_gamma_rot(Particle const *p, double const *&res);
#endif
#endif // LANGEVIN_PER_PARTICLE
#ifdef ROTATION
#endif

#ifdef ENGINE
void pointer_to_swimming(Particle const *p,
                         ParticleParametersSwimming const *&swim);
#endif

#ifdef ROTATIONAL_INERTIA
void pointer_to_rotational_inertia(Particle const *p, double const *&res);
#endif

/**
 * @brief Check if particle exists.
 *
 * @param part Id of the particle
 * @return True iff the particle exists.
 */
bool particle_exists(int part);

/**
 *  @brief Get the mpi rank which owns the particle with id.
 *
 *  @param id Id of the particle
 *  @return The MPI rank the particle is on.
 */
int get_particle_node(int id);

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

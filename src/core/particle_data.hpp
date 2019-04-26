/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
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

#include <utils/List.hpp>
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
/** \ref ParticleProperties::ext_flag "ext_flag" value for particle subject to
 *  an external force
 */
#define PARTICLE_EXT_FORCE 1
/** \ref ParticleProperties::ext_flag "ext_flag" value for fixed coordinate
 *  coord. */
#define COORD_FIXED(coord) (2L << (coord))
/** \ref ParticleProperties::ext_flag "ext_flag" mask to check whether any of
 *  the coordinates is fixed. */
#define COORDS_FIX_MASK (COORD_FIXED(0) | COORD_FIXED(1) | COORD_FIXED(2))
/** \ref ParticleProperties::ext_flag "ext_flag" mask to check whether all of
 *  the coordinates are fixed. */
#define COORDS_ALL_FIXED (COORD_FIXED(0) & COORD_FIXED(1) & COORD_FIXED(2))

#ifdef ROTATION
/** \ref ParticleProperties::ext_flag "ext_flag" value for particle subject to
 *  an external torque. */
#define PARTICLE_EXT_TORQUE 16
#endif

#endif

/************************************************
 * data types
 ************************************************/

/** Properties of a particle which are not supposed to
 *  change during the integration, but have to be known
 *  for all ghosts. Ghosts are particles which are
 *  needed in the interaction calculation, but are just copies of
 *  particles stored on different nodes.
 */
struct ParticleProperties {
  /** unique identifier for the particle. */
  int identity = -1;
  /** Molecule identifier. */
  int mol_id = 0;
  /** particle type, used for non bonded interactions. */
  int type = 0;

#ifdef MASS
  /** particle mass */
  double mass = 1.0;
#else
  constexpr static double mass{1.0};
#endif /* MASS */

#ifdef ROTATIONAL_INERTIA
  /** rotational inertia */
  Utils::Vector3d rinertia = {1., 1., 1.};
#else
  static constexpr Utils::Vector3d rinertia = {1., 1., 1.};
#endif

#ifdef AFFINITY
  /** parameters for affinity mechanisms */
  Utils::Vector3d bond_site = {-1., -1., -1.};
#endif

#ifdef MEMBRANE_COLLISION
  /** parameters for membrane collision mechanisms */
  Utils::Vector3d out_direction = {0., 0., 0.};
#endif

  // Determines, whether a particle's rotational degrees of freedom are
  // integrated
  int rotation = 0;

  /** charge. */
#ifdef ELECTROSTATICS
  double q = 0.0;
#else
  constexpr static double q{0.0};
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  /** electrophoretic mobility times E-field: mu_0 * E */
  Utils::Vector3d mu_E = {0., 0., 0.};
#endif

#ifdef DIPOLES
  /** dipole moment (absolute value)*/
  double dipm = 0.;
#endif

#ifdef VIRTUAL_SITES
  /** is particle virtual
      0 = real particle
      else = virtual particle */
  int is_virtual = 0;
#ifdef VIRTUAL_SITES_RELATIVE
  /** In case, the "relative" implementation of virtual sites is enabled, the
  following properties define, with respect to which real particle a virtual
  site is placed and in what distance. The relative orientation of the vector
  pointing from real particle to virtual site with respect to the orientation
  of the real particle is stored in the virtual site's quaternion attribute.
  */
  struct VirtualSitesRelativeParameteres {
    int to_particle_id = 0;
    double distance = 0;
    // Store relative position of the virtual site.
    Utils::Vector4d rel_orientation = {0., 0., 0., 0.};
    // Store the orientation of the virtual particle in the body fixed frame.
    Utils::Vector4d quat = {0., 0., 0., 0.};

    template <class Archive> void serialize(Archive &ar, long int) {
      ar &to_particle_id;
      ar &distance;
      ar &rel_orientation;
      ar &quat;
    }
  } vs_relative;

#endif
#else  /* VIRTUAL_SITES */
  static constexpr const int is_virtual = 0;
#endif /* VIRTUAL_SITES */

#ifdef LANGEVIN_PER_PARTICLE
  double T = -1.;
#ifndef PARTICLE_ANISOTROPY
  double gamma = -1.;
#else
  Utils::Vector3d gamma = {-1., -1., -1.};
#endif // PARTICLE_ANISOTROPY
/* Friction coefficient gamma for rotation */
#ifdef ROTATION
#ifndef PARTICLE_ANISOTROPY
  double gamma_rot = -1.;
#else
  Utils::Vector3d gamma_rot = {-1., -1., -1.};
#endif // ROTATIONAL_INERTIA
#endif // ROTATION
#endif // LANGEVIN_PER_PARTICLE

#ifdef SWIMMER_REACTIONS
  int catalyzer_count = 0;
#endif

#ifdef EXTERNAL_FORCES
  /** flag whether to fix a particle in space.
      Values:
      <ul> <li> 0 no external influence
           <li> 1 apply external force \ref ParticleProperties::ext_force
           <li> 2,3,4 fix particle coordinate 0,1,2
           <li> 5 apply external torque \ref ParticleProperties::ext_torque
      </ul>
  */
  int ext_flag = 0;
  /** External force, apply if \ref ParticleProperties::ext_flag == 1. */
  Utils::Vector3d ext_force = {0, 0, 0};

#ifdef ROTATION
  /** External torque, apply if \ref ParticleProperties::ext_flag == 16. */
  Utils::Vector3d ext_torque = {0, 0, 0};
#endif
#endif
};

/** Positional information on a particle. Information that is
    communicated to calculate interactions with ghost particles. */
struct ParticlePosition {
  /** periodically folded position. */
  Utils::Vector3d p = {0, 0, 0};

#ifdef ROTATION
  /** quaternions to define particle orientation */
  Utils::Vector4d quat = {1., 0., 0., 0.};
  /** unit director calculated from the quaternions */
  inline const Utils::Vector3d calc_director() const {
    return {2 * (quat[1] * quat[3] + quat[0] * quat[2]),
            2 * (quat[2] * quat[3] - quat[0] * quat[1]),
            quat[0] * quat[0] - quat[1] * quat[1] - quat[2] * quat[2] +
                quat[3] * quat[3]};
  };
#endif

#ifdef BOND_CONSTRAINT
  /**stores the particle position at the previous time step*/
  Utils::Vector3d p_old = {0., 0., 0.};
#endif
};

/** Force information on a particle. Forces of ghost particles are
    collected and added up to the force of the original particle. */
struct ParticleForce {
  ParticleForce() = default;
  ParticleForce(ParticleForce const &) = default;
  ParticleForce(const Utils::Vector3d &f) : f(f) {}
#ifdef ROTATION
  ParticleForce(const Utils::Vector3d &f, const Utils::Vector3d &torque)
      : f(f), torque(torque) {}
#endif

  ParticleForce &operator+=(ParticleForce const &rhs) {
    f += rhs.f;
#ifdef ROTATION
    torque += rhs.torque;
#endif

    return *this;
  }

  /** force. */
  Utils::Vector3d f = {0., 0., 0.};

#ifdef ROTATION
  /** torque */
  Utils::Vector3d torque = {0., 0., 0.};
#endif
};

/** Momentum information on a particle. Information not contained in
    communication of ghost particles so far, but a communication would
    be necessary for velocity dependent potentials. */
struct ParticleMomentum {
  /** velocity. */
  Utils::Vector3d v = {0., 0., 0.};

#ifdef ROTATION
  /** angular velocity
      ALWAYS IN PARTICLE FIXED, I.E., CO-ROTATING COORDINATE SYSTEM */
  Utils::Vector3d omega = {0., 0., 0.};
#endif
};

/** Information on a particle that is needed only on the
 *  node the particle belongs to
 */
struct ParticleLocal {
  /** position in the last time step before last Verlet list update. */
  Utils::Vector3d p_old = {0, 0, 0};
  /** index of the simulation box image where the particle really sits. */
  Utils::Vector3i i = {0, 0, 0};

  /** check whether a particle is a ghost or not */
  int ghost = 0;
};

struct ParticleParametersSwimming {
// ifdef inside because we need this type for some MPI prototypes
#ifdef ENGINE
  bool swimming = false;
  double f_swim = 0.;
  double v_swim = 0.;
#if defined(LB) || defined(LB_GPU)
  int push_pull = 0;
  double dipole_length = 0.;
  Utils::Vector3d v_center;
  Utils::Vector3d v_source;
  double rotational_friction = 0.;
#endif
#endif

  template <typename Archive> void serialize(Archive &ar, long int) {
#ifdef ENGINE
    ar &swimming &f_swim &v_swim
#if defined(LB) || defined(LB_GPU)
        &push_pull &dipole_length &v_center &v_source &rotational_friction
#endif
        ;
#endif
  }
};

/** Struct holding all information for one particle. */
struct Particle {
  int &identity() { return p.identity; }
  int const &identity() const { return p.identity; }

  bool operator==(Particle const &rhs) const {
    return identity() == rhs.identity();
  }

  bool operator!=(Particle const &rhs) const {
    return identity() != rhs.identity();
  }

  /**
   * @brief Return a copy of the particle with
   *        only the fixed size parts.
   *
   * This creates a copy of the particle with
   * only the parts than can be copied w/o heap
   * allocation, e.g. w/o bonds and exclusions.
   * This is more efficient if these parts are
   * not actually needed.
   */
  Particle flat_copy() const {
    Particle ret;

    ret.p = p;
    ret.r = r;
    ret.m = m;
    ret.f = f;
    ret.l = l;
#ifdef ENGINE
    ret.swim = swim;
#endif

    return ret;
  }

  ///
  ParticleProperties p;
  ///
  ParticlePosition r;
#ifdef DIPOLES
  inline const Utils::Vector3d calc_dip() const {
    return r.calc_director() * p.dipm;
  }
#endif
  ///
  ParticleMomentum m;
  ///
  ParticleForce f;
  ///
  ParticleLocal l;

  /** Bonded interactions list
   *
   *  The format is pretty simple: just the bond type, and then the particle
   *  ids. The number of particle ids can be determined easily from the
   *  bonded_ia_params entry for the type.
   */
  IntList bl;

  IntList &bonds() { return bl; }
  IntList const &bonds() const { return bl; }

  IntList &exclusions() {
#ifdef EXCLUSIONS
    return el;
#else
    throw std::runtime_error{"Exclusions not enabled."};
#endif
  }

  IntList const &exclusions() const {
#ifdef EXCLUSIONS
    return el;
#else
    throw std::runtime_error{"Exclusions not enabled."};
#endif
  }

#ifdef EXCLUSIONS
  /** list of particles, with which this particle has no nonbonded
   *  interactions
   */
  IntList el;
#endif

#ifdef ENGINE
  ParticleParametersSwimming swim;
#endif
};

/**
 * These functions cause a compile time error if
 * Particles are copied by memmove or memcpy,
 * which do not keep class invariants.
 *
 * These are templates so that the error is caused
 * at the place they are used.
 */
template <typename Size> void memmove(Particle *, Particle *, Size) {
  static_assert(sizeof(Size) == 0, "Particles can not be copied like this.");
}
template <typename Size> void memmove(Particle *, Particle const *, Size) {
  static_assert(sizeof(Size) == 0, "Particles can not be copied like this.");
}

template <typename Size> void memcpy(Particle *, Particle *, Size) {
  static_assert(sizeof(Size) == 0, "Particles can not be copied like this.");
}
template <typename Size> void memcpy(Particle *, Particle const *, Size) {
  static_assert(sizeof(Size) == 0, "Particles can not be copied like this.");
}

template <typename Size, typename... Ts>
void MPI_Send(Particle *, Size, Ts...) {
  static_assert(sizeof(Size) == 0, "Particles can not be copied like this.");
}

template <typename Size, typename... Ts>
void MPI_Send(Particle const *, Size, Ts...) {
  static_assert(sizeof(Size) == 0, "Particles can not be copied like this.");
}

/** List of particles. The particle array is resized using a sophisticated
 *  (we hope) algorithm to avoid unnecessary resizes.
 *  Access using \ref realloc_particlelist, ...
 */
struct ParticleList {
  ParticleList() : part{nullptr}, n{0}, max{0} {}
  /** The particles payload */
  Particle *part;
  /** Number of particles contained */
  int n;
  /** Number of particles that fit in until a resize is needed */
  int max;
};

/************************************************
 * exported variables
 ************************************************/

/** Highest particle number seen so far. If you leave out some
    particle numbers, this number might be higher than the
    true number of particles. On the other hand, if you start
    your particle numbers at 0, the total number of particles
    is larger by 1.
*/
extern int max_seen_particle;
/** total number of particles on all nodes. */
extern int n_part;
/** flag that active swimming particles exist */
extern bool swimming_particles_exist;

/** id->particle mapping on all nodes. This is used to find partners
    of bonded interactions. */
extern Particle **local_particles;
extern int max_local_particles;

/************************************************
 * Functions
 ************************************************/

/*       Functions acting on Particles          */
/************************************************/

/** Deallocate the dynamic storage of a particle. */
void free_particle(Particle *part);

/*    Functions acting on Particle Lists        */
/************************************************/

/** Initialize a particle list.
 *  Use with care and ONLY for initialization! */
void init_particlelist(ParticleList *pList);

/** Allocate storage for local particles and ghosts. This version
    does \em not care for the bond information to be freed if necessary.
    \param plist the list on which to operate
    \param size the size to provide at least. It is rounded
    up to multiples of \ref PART_INCREMENT.
    \return true iff particle addresses have changed */
int realloc_particlelist(ParticleList *plist, int size);

/** Append a particle at the end of a particle List.
    reallocates particles if necessary!
    This procedure does not care for \ref local_particles.
    \param l List to append the particle to.
    \param part  Particle to append. */
void append_unindexed_particle(ParticleList *l, Particle &&part);

/** Append a particle at the end of a particle List.
    reallocates particles if necessary!
    This procedure cares for \ref local_particles.
    \param plist List to append the particle to.
    \param part  Particle to append.
    \return Pointer to new location of the particle. */
Particle *append_indexed_particle(ParticleList *plist, Particle &&part);

/** Remove a particle from one particle List and append it to another.
    Refill the sourceList with last particle and update its entry in
    local_particles. reallocates particles if necessary.  This
    procedure does not care for \ref local_particles.
    NOT IN USE AT THE MOMENT.
    \param destList   List where the particle is appended.
    \param sourceList List where the particle will be removed.
    \param ind        Index of the particle in the sourceList.
    \return Pointer to new location of the particle.
 */
Particle *move_unindexed_particle(ParticleList *destList,
                                  ParticleList *sourceList, int ind);

/** Remove a particle from one particle List and append it to another.
    Refill the sourceList with last particle and update its entry in
    local_particles. Reallocates particles if necessary.  This
    procedure cares for \ref local_particles.
    \param destList   List where the particle is appended.
    \param sourceList List where the particle will be removed.
    \param ind        Index of the particle in the sourceList.
    \return Pointer to new location of the particle.
 */
Particle *move_indexed_particle(ParticleList *destList,
                                ParticleList *sourceList, int ind);

Particle extract_indexed_particle(ParticleList *sl, int i);

/*    Other Functions                           */
/************************************************/

/** Update the entries in \ref local_particles for all particles in the list pl.
 *  @param pl the list to put in.
 */
void update_local_particles(ParticleList *pl);

/** Invalidate \ref particle_node. This has to be done
 *  at the beginning of the integration.
 */
void clear_particle_node();

/** Realloc \ref local_particles. */
void realloc_local_particles(int part);

/**
 * @brief Get particle data.
 *
 *   @param part the identity of the particle to fetch
 *   @return Pointer to copy of particle if it exists,
 *          nullptr otherwise;
 */
const Particle &get_particle_data(int part);

/**
 * @brief Fetch a range of particle into the fetch cache.
 *
 * If the range is larger than the cache size, only
 * the particle that fit into the cache are fetched.
 *
 * @param ids Ids of the particles that should be
 *        fetched.
 */
void prefetch_particle_data(std::vector<int> ids);

/** @brief Invalidate the fetch cache for get_particle_data. */
void invalidate_fetch_cache();

/** Call only on the master node.
 *  Move a particle to a new position.
 *  If it does not exist, it is created.
 *  @param part the identity of the particle to move
 *  @param p    its new position
 *  @retval ES_PART_OK if particle existed
 *  @retval ES_PART_CREATED if created
 *  @retval ES_PART_ERROR if id is illegal
 */
int place_particle(int part, double p[3]);

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

/** Call only on the master node: set particle solvation free energy.
 *  @param part the particle.
 *  @param solvation its new solvation free energy.
 */
void set_particle_solvation(int part, double *solvation);

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

#ifdef AFFINITY
/** Call only on the master node: set particle affinity.
 *  @param part the particle.
 *  @param bond_site its new site of the affinity bond.
 */
void set_particle_affinity(int part, double *bond_site);
#endif

#ifdef MEMBRANE_COLLISION
/** Call only on the master node: set particle out_direction.
 *  @param part the particle.
 *  @param out_direction its new outward direction with respect to membrane.
 */
void set_particle_out_direction(int part, double *out_direction);
#endif

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
void set_particle_mu_E(int part, double *mu_E);
void get_particle_mu_E(int part, double (&mu_E)[3]);
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
/** Call only on the master node: set particle dipole moment (absolute value).
 *  @param part the particle.
 *  @param is_virtual its new is_virtual.
 */
void set_particle_virtual(int part, int is_virtual);
#endif
#ifdef VIRTUAL_SITES_RELATIVE
void set_particle_vs_quat(int part, double *vs_relative_quat);
void set_particle_vs_relative(int part, int vs_relative_to, double vs_distance,
                              double *rel_ori);
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
void set_particle_fix(int part, int flag);
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

/** For all local particles, remove bonds incorporating the specified particle.
 *  @param part     identity of the particle to free from bonds
 */
void remove_all_bonds_to(int part);

/** Used by \ref mpi_place_particle, should not be used elsewhere.
 *  Move a particle to a new position. If it does not exist, it is created.
 *  The position must be on the local node!
 *
 *  @param part the identity of the particle to move
 *  @param p    its new position
 *  @param _new  if true, the particle is allocated, else has to exists already
 *
 *  @return Pointer to the particle.
 */
Particle *local_place_particle(int part, const double p[3], int _new);

/** Used by \ref mpi_place_particle, should not be used elsewhere.
 *  Called if on a different node a new particle was added.
 *  @param part the identity of the particle added
 */
void added_particle(int part);

/** Used for example by \ref mpi_send_exclusion.
 *  Locally add a exclusion to a particle.
 *  @param part1 the identity of the first exclusion partner
 *  @param part2 the identity of the second exclusion partner
 *  @param _delete if true, delete the exclusion instead of add
 */
void local_change_exclusion(int part1, int part2, int _delete);

/** Used by \ref mpi_remove_particle, should not be used elsewhere.
 *  Remove a particle on this node.
 *  @param part the identity of the particle to remove
 */
void local_remove_particle(int part);

/** Used by \ref mpi_remove_particle, should not be used elsewhere.
 *  Locally remove all particles.
 */
void local_remove_all_particles();

/** Used by \ref mpi_rescale_particles, should not be used elsewhere.
 *  Locally rescale all particles on current node.
 *  @param dir   direction to scale (0/1/2 = x/y/z, 3 = x+y+z isotropically)
 *  @param scale factor by which to rescale (>1: stretch, <1: contract)
 */
void local_rescale_particles(int dir, double scale);

/** @brief Add bond to local particle.
 *  @param p     identity of principal atom of the bond.
 *  @param bond  field containing the bond type number and the identity
 *               of all bond partners (secondary atoms of the bond).
 */
void local_add_particle_bond(Particle &p, Utils::Span<const int> bond);

/** Synchronous send of a particle buffer to another node.
 *  The other node MUST call \ref recv_particles when this is called.
 *  The particles data is freed.
 */
void send_particles(ParticleList *particles, int node);

/** Synchronous receive of a particle buffer from another node.
 *  The other node MUST call \ref send_particles when this is called.
 *  Particles needs to initialized, it is reallocated to the correct
 *  size and the content is overwritten.
 */
void recv_particles(ParticleList *particles, int node);

#ifdef EXCLUSIONS
/** Determine if the non-bonded interactions between @p p1 and @p p2 should be
 *  calculated.
 */
inline bool do_nonbonded(Particle const *p1, Particle const *p2) {
  /* check for particle 2 in particle 1's exclusion list. The exclusion list is
     symmetric, so this is sufficient. */
  return std::none_of(p1->el.begin(), p1->el.end(),
                      [p2](int id) { return p2->p.identity == id; });
}
#endif

/** Remove bond from particle if possible */
int try_delete_bond(Particle *part, const int *bond);

/** Remove exclusion from particle if possible */
void try_delete_exclusion(Particle *part, int part2);

/** Insert an exclusion if not already set */
void try_add_exclusion(Particle *part, int part2);

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
int get_random_p_id(int type);
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
void pointer_to_virtual(Particle const *p, int const *&res);
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void pointer_to_vs_quat(Particle const *p, double const *&res);
void pointer_to_vs_relative(Particle const *p, int const *&res1,
                            double const *&res2, double const *&res3);
#endif

void pointer_to_dip(Particle const *P, double const *&res);

void pointer_to_dipm(Particle const *P, double const *&res);

#ifdef EXTERNAL_FORCES
void pointer_to_ext_force(Particle const *p, int const *&res1,
                          double const *&res2);
#ifdef ROTATION
void pointer_to_ext_torque(Particle const *p, int const *&res1,
                           double const *&res2);
#endif
void pointer_to_fix(Particle const *p, int const *&res);
#endif

#ifdef LANGEVIN_PER_PARTICLE
void pointer_to_gamma(Particle const *p, double const *&res);
void pointer_to_temperature(Particle const *p, double const *&res);
#ifdef ROTATION
void pointer_to_gamma_rot(Particle const *p, double const *&res);
#endif
#endif // LANGEVIN_PER_PARTICLE
#ifdef ROTATION
void pointer_to_rotation(Particle const *p, int const *&res);
#endif

#ifdef ENGINE
void pointer_to_swimming(Particle const *p,
                         ParticleParametersSwimming const *&swim);
#endif

#ifdef ROTATIONAL_INERTIA
void pointer_to_rotational_inertia(Particle const *p, double const *&res);
#endif
#ifdef AFFINITY
void pointer_to_bond_site(Particle const *p, double const *&res);
#endif

#ifdef MEMBRANE_COLLISION
void pointer_to_out_direction(const Particle *p, const double *&res);
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

#endif

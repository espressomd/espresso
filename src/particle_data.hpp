/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
/** \file particle_data.hpp
    For more information on particle_data,
    see \ref particle_data.cpp "particle_data.c"
*/


#include "utils.hpp"
#include "global.hpp"

/************************************************
 * defines
 ************************************************/

/// ok code for \ref place_particle
#define ES_PART_OK 0
/// error code for \ref place_particle
#define ES_PART_ERROR -1
/// ok code for \ref place_particle, particle is new
#define ES_PART_CREATED 1

/**  bonds_flag "bonds_flag" value for updating particle config without bonding information */
#define WITHOUT_BONDS 0
/**  bonds_flag "bonds_flag" value for updating particle config with bonding information */
#define WITH_BONDS 1

#ifdef EXTERNAL_FORCES
/** \ref ParticleLocal::ext_flag "ext_flag" value for particle subject to an external force. */
#define PARTICLE_EXT_FORCE 1
/** \ref ParticleLocal::ext_flag "ext_flag" value for fixed coordinate coord. */
#define COORD_FIXED(coord) (2L << coord)
/** \ref ParticleLocal::ext_flag "ext_flag" mask to check wether any of the coordinates is fixed. */
#define COORDS_FIX_MASK     (COORD_FIXED(0) | COORD_FIXED(1) | COORD_FIXED(2))

#ifdef ROTATION
/** \ref ParticleLocal::ext_flag "ext_flag" value for particle subject to an external torque. */
#define PARTICLE_EXT_TORQUE 16
#endif

#endif


/************************************************
 * data types
 ************************************************/

/** Properties of a particle which are not supposed to
    change during the integration, but have to be known
    for all ghosts. Ghosts are particles which are
    needed in the interaction calculation, but are just copies of
    particles stored on different nodes.
*/
typedef struct {
  /** unique identifier for the particle. */
  int    identity;
  /** Molecule identifier. */
  int    mol_id;
  /** particle type, used for non bonded interactions. */
  int    type;

#ifdef MASS
  /** particle mass */
  double mass;
#endif

#ifdef SHANCHEN
  double solvation[2*LB_COMPONENTS];
#endif

#ifdef ROTATIONAL_INERTIA
  /** rotational inertia */
  double rinertia[3];
#endif

#ifdef ROTATION_PER_PARTICLE
  // Determines, wether a particle's rotational degrees of freedom are
  // integrated
  int rotation;
#endif

#ifdef ELECTROSTATICS
  /** charge. */
  double q;
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  /** electrophoretic mobility times E-field: mu_0 * E */
  double mu_E[3];
#endif

#ifdef DIPOLES
  /** dipole moment (absolute value)*/
  double dipm;
#endif

#ifdef VIRTUAL_SITES
  /** is particle virual
      0 = real particle
      else = virual particle */
  int isVirtual;
  #ifdef VIRTUAL_SITES_RELATIVE
  /** In case, the "relative" implementation of virtual sites is enabled, the 
  following properties define, with respect to which real particle a virtual
  site is placed and in what distance. The relative orientation of the vector
  pointing from real particle to virtual site with respect to the orientation
  of the real particle is stored in the virtual site's quaternion attribute.
  */
  int vs_relative_to_particle_id;
  double vs_relative_distance;
  #endif
#endif

#ifdef LANGEVIN_PER_PARTICLE
  double T;
  double gamma;
#endif

#ifdef CATALYTIC_REACTIONS
  int catalyzer_count;
#endif
} ParticleProperties;

/** Positional information on a particle. Information that is
    communicated to calculate interactions with ghost particles. */
typedef struct {
  /** periodically folded position. */
  double p[3];

#ifdef ROTATION
  /** quaternions to define particle orientation */
  double quat[4];
  /** unit director calculated from the quaternions */
  double quatu[3];
#endif

#ifdef DIPOLES
  /** dipol moment. This is synchronized with quatu and quat. */
  double dip[3];
#endif

#ifdef BOND_CONSTRAINT
  /**stores the particle position at the previous time step*/
  double p_old[3];
#endif

#ifdef SHANCHEN
  double composition[LB_COMPONENTS];
#endif

} ParticlePosition;

/** Force information on a particle. Forces of ghost particles are
    collected and added up to the force of the original particle. */
typedef struct {
  /** force. */
  double f[3];

#ifdef ROTATION
  /** torque */
  double torque[3];
#endif

} ParticleForce;

/** Momentum information on a particle. Information not contained in
    communication of ghost particles so far, but a communication would
    be necessary for velocity dependend potentials. */
typedef struct {
  /** velocity. */
  double v[3];

#ifdef ROTATION
  /** angular velocity  
      ALWAYS IN PARTICLE FIXEXD, I.E., CO-ROTATING COORDINATE SYSTEM */
  double omega[3];
#endif
} ParticleMomentum;

/** Information on a particle that is needed only on the
    node the particle belongs to */
typedef struct {
  /** position in the last time step befor last Verlet list update. */
  double p_old[3];
  /** index of the simulation box image where the particle really sits. */
  int    i[3];

#ifdef EXTERNAL_FORCES
  /** flag whether to fix a particle in space.
      Values:
      <ul> <li> 0 no external influence
           <li> 1 apply external force \ref ParticleLocal::ext_force
           <li> 2,3,4 fix particle coordinate 0,1,2
           <li> 5 apply external torque \ref ParticleLocal::ext_torque
      </ul>
  */
  int ext_flag;
  /** External force, apply if \ref ParticleLocal::ext_flag == 1. */
  double ext_force[3];

  #ifdef ROTATION
  /** External torque, apply if \ref ParticleLocal::ext_flag == 16. */
  double ext_torque[3];
  #endif

#endif

#ifdef GHOST_FLAG
  /** check whether a particle is a ghost or not */
  int ghost;
#endif

#ifdef GHMC
  /** Data for the ghmc thermostat, last saved 
      position and monentum of particle */
  ParticlePosition r_ls;
  ParticleMomentum m_ls;
#endif
} ParticleLocal;

#ifdef LB
/** Data related to the Lattice Boltzmann hydrodynamic coupling */
typedef struct {
  /** fluctuating part of the coupling force */
  double f_random[3];
} ParticleLatticeCoupling;
#endif

/** Struct holding all information for one particle. */
typedef struct {
  ///
  ParticleProperties p;
  ///
  ParticlePosition r;
  ///
  ParticleMomentum m;
  ///
  ParticleForce f;
  ///
  ParticleLocal l;
  ///
#ifdef LB
  ParticleLatticeCoupling lc;
#endif
  /** bonded interactions list. The format is pretty simple:
      Just the bond type, and then the particle ids. The number of particle ids can be determined
      easily from the bonded_ia_params entry for the type. */
  IntList bl;

#ifdef EXCLUSIONS
  /** list of particles, with which this particle has no nonbonded interactions */
  IntList el;
#endif

} Particle;

/** List of particles. The particle array is resized using a sophisticated
    (we hope) algorithm to avoid unnecessary resizes.
    Access using \ref realloc_particlelist, \ref got_particle,...
*/
typedef struct {
  /** The particles payload */
  Particle *part;
  /** Number of particles contained */
  int n;
  /** Number of particles that fit in until a resize is needed */
  int max;
} ParticleList;

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
extern int  n_part;

/** Capacity of the \ref particle_node / \ref local_particles. */
extern int  max_particle_node;
/** Used only on master node: particle->node mapping. */
extern int  *particle_node;
/** id->particle mapping on all nodes. This is used to find partners
    of bonded interactions. */
extern Particle   **local_particles;

/** Particles' current configuration. Before using that
    call \ref updatePartCfg or \ref sortPartCfg to allocate
    the data if necessary (which is decided by \ref updatePartCfg). */
extern Particle *partCfg;

/** if non zero, \ref partCfg is sorted by particle order, and
    the particles are stored consecutively starting with 0. */
extern int partCfgSorted;

/** Particles' current bond partners. \ref partBondPartners is
    sorted by particle order, and the particles are stored
    consecutively starting with 0. This array is global to all nodes*/
extern int *partBondPartners;


/************************************************
 * Functions
 ************************************************/


/*       Functions acting on Particles          */
/************************************************/

/** Initialize a particle.
    This function just sets all values to the defaults (mostly zeros)!
    Do NOT use this without setting the values of the
    \ref ParticleProperties::identity "identity" and \ref ParticlePosition::p "position" to
    reasonable values. Also make sure that you update \ref local_particles.

    Add here all initializations you need to be done !!!
    If dynamic memory allocation is involved, also look at \ref free_particle.
*/
void init_particle(Particle *part);

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
    \return true iff particle adresses have changed */
int realloc_particlelist(ParticleList *plist, int size);

/** Search for a specific particle.
    \param plist the list on which to operate
    \param id the identity of the particle to search
    \return a pointer to the particle structure or NULL if particle is
    not in this list */
Particle *got_particle(ParticleList *plist, int id);

/** Append a particle at the end of a particle List.
    reallocates particles if necessary!
    This procedure does not care for \ref local_particles.
    \param plist List to append the particle to.
    \param part  Particle to append.
    \return Pointer to new location of the particle. */
Particle *append_unindexed_particle(ParticleList *plist, Particle *part);

/** Append a particle at the end of a particle List.
    reallocates particles if necessary!
    This procedure cares for \ref local_particles.
    \param plist List to append the particle to.
    \param part  Particle to append.
    \return Pointer to new location of the particle. */
Particle *append_indexed_particle(ParticleList *plist, Particle *part);

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
Particle *move_unindexed_particle(ParticleList *destList, ParticleList *sourceList, int ind);

/** Remove a particle from one particle List and append it to another.
    Refill the sourceList with last particle and update its entry in
    local_particles. Reallocates particles if necessary.  This
    procedure cares for \ref local_particles.
    \param destList   List where the particle is appended.
    \param sourceList List where the particle will be removed.
    \param ind        Index of the particle in the sourceList.
    \return Pointer to new location of the particle.
 */
Particle *move_indexed_particle(ParticleList *destList, ParticleList *sourceList, int ind);

/*    Other Functions                           */
/************************************************/

/** Update the entries in \ref local_particles for all particles in the list pl.
    @param pl the list to put in.
*/
void update_local_particles(ParticleList *pl);

/** Rebuild \ref particle_node from scratch.
    After a simulation step \ref particle_node has to be rebuild
    since the particles might have gone to a different node.
*/
void build_particle_node();

/** Invalidate \ref particle_node. This has to be done
    at the beginning of the integration.
*/
void particle_invalidate_part_node();

/** Realloc \ref local_particles. */
void realloc_local_particles();

/** Get particle data. Note that the bond intlist is
    allocated so that you are responsible to free it later.
    @param part the identity of the particle to fetch
    @param data where to store its contents.
    @return ES_OK if particle existed
*/
int get_particle_data(int part, Particle *data);

/** Call only on the master node.
    Move a particle to a new position.
    If it does not exist, it is created.
    @param part the identity of the particle to move
    @param p    its new position
    @return ES_PART_OK if particle existed, ES_PART_CREATED
    if created and ES_PART_ERROR if id is illegal
*/
int place_particle(int part, double p[3]);

/** Call only on the master node: set particle velocity.
    @param part the particle.
    @param v its new velocity.
    @return ES_OK if particle existed
*/
int set_particle_v(int part, double v[3]);

/** Call only on the master node: set particle force.
    @param part the particle.
    @param F its new force.
    @return ES_OK if particle existed
*/
int set_particle_f(int part, double F[3]);

/** Call only on the master node: set particle mass.
    @param part the particle.
    @param mass its new mass.
    @return ES_OK if particle existed
*/
int set_particle_mass(int part, double mass);

/** Call only on the master node: set particle solvation free energy.
    @param part the particle.
    @param solvation its new solvation free energy.
    @return ES_OK if particle existed
*/
int set_particle_solvation(int part, double* solvation);


#ifdef ROTATIONAL_INERTIA
/** Call only on the master node: set particle rotational inertia.
    @param part the particle.
    @param rinertia its new inertia.
    @return ES_OK if particle existed
*/
int set_particle_rotational_inertia(int part, double rinertia[3]);
#endif

#ifdef ROTATION_PER_PARTICLE
/** Call only on the master node: Specifies whether a particle's rotational
    degrees of freedom are integrated or not. If set to zero, the content of
    the torque and omega variables are meaningless
    @param part the particle.
    @param rot the degrees of freedom flag.
    @return ES_OK if particle existed
*/
int set_particle_rotation(int part, int rot);
#endif


/** Call only on the master node: set particle charge.
    @param part the particle.
    @param q its new charge.
    @return ES_OK if particle existed
*/
int set_particle_q(int part, double q);

/** Call only on the master node: set particle electrophoretic mobility.
    @param part the particle.
    @param mu_E its new mobility.
    @return ES_OK if particle existed
*/
int set_particle_mu_E(int part, double mu_E[3]);

/** Call only on the master node: set particle type.
    @param part the particle.
    @param type its new type.
    @return ES_OK if particle existed
*/
int set_particle_type(int part, int type);

/** Call only on the master node: set particle's molecule id.
    @param part the particle.
    @param mid  its new mol id.
    @return ES_OK if particle existed
*/
int set_particle_mol_id(int part, int mid);

#ifdef ROTATION
/** Call only on the master node: set particle orientation using quaternions.
    @param part the particle.
    @param quat its new value for quaternions.
    @return ES_OK if particle existed
*/
int set_particle_quat(int part, double quat[4]);

/** Call only on the master node: set particle angular velocity from lab frame.
    @param part the particle.
    @param omega its new angular velocity.
    @return ES_OK if particle existed
*/
int set_particle_omega_lab(int part, double omega[3]);

/** Call only on the master node: set particle angular velocity in body frame.
    @param part the particle.
    @param omega its new angular velocity.
    @return ES_OK if particle existed
*/
int set_particle_omega_body(int part, double omega[3]);

/** Call only on the master node: set particle torque from lab frame.
    @param part the particle.
    @param torque its new torque.
    @return ES_OK if particle existed
*/
int set_particle_torque_lab(int part, double torque[3]);

/** Call only on the master node: set particle torque in body frame.
    @param part the particle.
    @param torque its new torque.
    @return ES_OK if particle existed
*/
int set_particle_torque_body(int part, double torque[3]);
#endif

#ifdef DIPOLES
/** Call only on the master node: set particle dipole orientation.
    @param part the particle.
    @param dip its new dipole orientation.
    @return ES_OK if particle existed
*/
int set_particle_dip(int part, double dip[3]);

/** Call only on the master node: set particle dipole moment (absolut value).
    @param part the particle.
    @param dipm its new dipole moment.
    @return ES_OK if particle existed
*/
int set_particle_dipm(int part, double dipm);
#endif

#ifdef VIRTUAL_SITES
/** Call only on the master node: set particle dipole moment (absolut value).
    @param part the particle.
    @param isVirtual its new dipole moment.
    @return ES_OK if particle existed
*/
int set_particle_virtual(int part,int isVirtual);
#endif

#ifdef LANGEVIN_PER_PARTICLE
/** Call only on the master node: set particle temperature.
    @param part the particle.
    @param T its new temperature.
    @return ES_OK if particle existed
*/
int set_particle_temperature(int part, double T);

/** Call only on the master node: set particle frictional coefficient.
    @param part the particle.
    @param gamma its new frictional coefficient.
    @return ES_OK if particle existed
*/
int set_particle_gamma(int part, double gamma);
#endif

#ifdef EXTERNAL_FORCES
  #ifdef ROTATION
    /** Call only on the master node: set particle external torque.
        @param part  the particle.
        @param flag  new value for ext_flag.
        @param torque new value for ext_torque.
        @return ES_OK if particle existed
    */
    int set_particle_ext_torque(int part, int flag, double torque[3]);
  #endif
/** Call only on the master node: set particle external force.
    @param part  the particle.
    @param flag  new value for ext_flag.
    @param force new value for ext_force.
    @return ES_OK if particle existed
*/
int set_particle_ext_force(int part, int flag, double force[3]);
/** Call only on the master node: set coordinate axes for which the particles motion is fixed.
    @param part  the particle.
    @param flag new value for flagged coordinate axes to be fixed
    @return ES_OK if particle existed
*/
int set_particle_fix(int part,  int flag);
#endif

/** Call only on the master node: change particle bond.
    @param part     identity of principal atom of the bond.
    @param bond     field containing the bond type number and the
    identity of all bond partners (secundary atoms of the bond). If NULL, delete all bonds.
    @param _delete   if true, do not add the bond, rather delete it if found
    @return ES_OK on success or ES_ERROR if no success
    (e. g. particle or bond to delete does not exist)
*/
int change_particle_bond(int part, int *bond, int _delete);

#ifdef EXCLUSIONS
/** Call only on the master node: change particle constraints.
    @param part     identity of particle for which the exclusion is set.
    @param part2    identity of particle for which the exclusion is set. If -1, delete all exclusions.
    @param _delete   if true, do not add the exclusion, rather delete it if found
    @return ES_OK on success or ES_ERROR if no success
    (e. g. particles do not exist / did not have exclusion set)
*/
int change_exclusion(int part, int part2, int _delete);

/** remove all exclusions. */
void remove_all_exclusions();
#endif

/** remove particle with a given identity. Also removes all bonds to the particle.
    @param part     identity of the particle to remove
    @return ES_OK on success or ES_ERROR if particle does not exist
*/
int remove_particle(int part);

/** remove all particles.
 */
void remove_all_particles();

/** for all local particles, remove bonds incorporating the specified
    particle.
    @param part     identity of the particle to free from bonds
*/
void remove_all_bonds_to(int part);

/** Get the complete unsorted informations on all particles into \ref
    partCfg if something's changed. This is a severe performance
    drawback and might even fail for lack of memory for large systems.
    If you need the particle info sorted, call \ref sortPartCfg
    instead.  This function is lazy. If you would like the bonding
    information in \ref partCfg to be valid you should set the value
    of  to \ref WITH_BONDS.
*/
int updatePartCfg(int bonds_flag );

/** release the partCfg array. Use this function, since it also frees the
    bonds, if they are used.
*/
void freePartCfg();

/** sorts the \ref partCfg array. This is indicated by setting
    \ref partCfgSorted to 1. Note that for this to work the particles
    have to be stored consecutively starting with 0.
    This function is lazy.
    @return 1 iff sorting was possible, i. e. the particles were stored
    consecutively.
*/
int sortPartCfg();

/** Used by \ref mpi_place_particle, should not be used elsewhere.
    Move a particle to a new position.
    If it does not exist, it is created. the position must
    be on the local node!
    @param part the identity of the particle to move
    @param p    its new position
    @param _new  if true, the particle is allocated, else has to exists already
*/
void local_place_particle(int part, double p[3], int _new);

/** Used by \ref mpi_place_particle, should not be used elsewhere.
    Called if on a different node a new particle was added.
    @param part the identity of the particle added
*/
void added_particle(int part);

/** Used by \ref mpi_send_bond, should not be used elsewhere.
    Modify a bond.
    @param part the identity of the particle to change
    @param bond the bond to do
    @param _delete if true, delete the bond instead of add
    @return ES_OK for add or successful delete, ES_ERROR else
*/
int local_change_bond(int part, int *bond, int _delete);

/** Used for example by \ref mpi_send_exclusion.
    Locally add a exclusion to a particle.
    @param part1 the identity of the first exclusion partner
    @param part2 the identity of the second exclusion partner
    @param _delete if true, delete the exclusion instead of add
*/
void local_change_exclusion(int part1, int part2, int _delete);

/** Used by \ref mpi_remove_particle, should not be used elsewhere.
    Remove a particle on this node.
    @param part the identity of the particle to remove
*/
void local_remove_particle(int part);

/** Used by \ref mpi_remove_particle, should not be used elsewhere.
    Locally remove all particles.
 */
void local_remove_all_particles();

/** Used by \ref mpi_rescale_particles, should not be used elsewhere.
    Locally rescale all particles on current node.
    @param dir   direction to scale (0/1/2 = x/y/z, 3 = x+y+z isotropically)
    @param scale factor by which to rescale (>1: stretch, <1: contract)
*/
void local_rescale_particles(int dir, double scale);

/** Synchronous send of a particle buffer to another node. The other node
    MUST call \ref recv_particles when this is called. The particles data
    is freed. */
void send_particles(ParticleList *particles, int node);

/** Synchronous receive of a particle buffer from another node. The other node
    MUST call \ref send_particles when this is called. The particles are
    APPENDED to the list, so it has to be a valid one */
void recv_particles(ParticleList *particles, int node);

#ifdef EXCLUSIONS
/** Determines if the non bonded interactions between p1 and p2 should be calculated */
inline int do_nonbonded(Particle *p1, Particle *p2)
{
  int i, i2;
  /* check for particle 2 in particle 1's exclusion list. The exclusion list is
     symmetric, so this is sufficient. */
  i2  = p2->p.identity;
  for (i = 0; i < p1->el.n; i++)
    if (i2 == p1->el.e[i]) return 0;
  return 1;
}
#endif

/** Remove bond from particle if possible */
int try_delete_bond(Particle *part, int *bond);

/** Remove exclusion from particle if possible */
void try_delete_exclusion(Particle *part, int part2);

/** Insert an exclusion if not already set */
void try_add_exclusion(Particle *part, int part2);

/** Automatically add the next \<distance\> neighbors in each molecule to the exclusion list.
 This uses the bond topology obtained directly from the particles, since only this contains
 the full topology, in contrast to \ref topology::topology. To easily setup the bonds, all data
 should be on a single node, therefore the \ref partCfg array is used. With large amounts
 of particles, you should avoid this function and setup exclusions manually. */
void auto_exclusion(int distance);

/* keep a unique list for particle i. Particle j is only added if it is not i
 and not already in the list. */
void add_partner(IntList *il, int i, int j, int distance);

//value that is returned in the case there was no error, but the type was not yet indexed
#define NOT_INDEXED -3
//struct that associates the index used for the type_list and the real particle type
typedef struct {
	int max_entry;
	int * type;
} IndexOfType;

//and the other way round
typedef struct {
	int max_entry;
	int * index;
} TypeOfIndex;

typedef struct {
	int max_entry;
	int cur_size;
	int *id_list;
} TypeList;

//undefined array size
extern TypeList *type_array;
extern int number_of_type_lists;

extern TypeOfIndex Type; 
extern IndexOfType Index; 

// flag indicating init_gc was called 
extern int GC_init;

// flag that indicates that the function init_type_array was called already
extern int Type_array_init;

int init_gc(void);

/** init particle lists		*/
int init_type_array(int type);

/** resize the array for the list of ids for a certain type */
int reallocate_type_array(int type);

/** make more type_arrays available */
int reallocate_global_type_list(int size);

/** free particle lists		*/
int free_particle_lists(void);

//update particle list
int update_particle_array(int type);

/* find a particle of given type and return its id */
int find_particle_type(int type, int *id);

/** return an array with real particle id and the corresponding index of typelist */
int find_particle_type_id(int type, int *id, int *in_id );

/** delete one randomly chosen particle of given type 
 * returns ES_OK if succesful or else ES_ERROR		*/
int delete_particle_of_type(int type);

int remove_id_type_array(int part_id, int type);
int add_particle_to_list(int part_id, int type);
// print out a list of currently indexed ids
int gc_status(int type);
int number_of_particles_with_type(int type, int *number);


// The following functions are used by the python interface to obtain 
// properties of a particle, which are only compiled in in some configurations
// This is needed, because cython does not support conditional compilation 
// within a ctypedef definition


#ifdef ROTATION
void pointer_to_omega_body(Particle* p, double*& res);

void pointer_to_torque_lab(Particle* p, double*& res);

void pointer_to_quat(Particle* p, double*& res);
void pointer_to_quatu(Particle* p, double*& res);

#endif

#ifdef ELECTROSTATICS
void pointer_to_q(Particle* p, double*& res);
#endif

#ifdef VIRTUAL_SITES
void pointer_to_virtual(Particle* p, int*& res);
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void pointer_to_vs_relative(Particle* p, int*& res1,double*& res2);
#endif

#ifdef MASS
void pointer_to_mass(Particle* p, double*&  res);
#endif

void pointer_to_dip(Particle* P, double*& res);

void pointer_to_dipm(Particle* P, double*& res);
#endif

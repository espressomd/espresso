// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef PARTICLE_DATA_H
#define PARTICLE_DATA_H
/** \file particle_data.h

    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>
    For more information on particle_data,
    see \ref particle_data.c "particle_data.c"
 */

#include <tcl.h>
#include "config.h"
#include "grid.h"

/************************************************
 * defines
 ************************************************/

#define DATA_DOMAIN_DECOMPOSITION 0
#define DATA_NSQUARE              1

#ifdef EXTERNAL_FORCES
/** \ref Particle::ext_flag "ext_flag" value for unfixed particle. */
#define PARTICLE_UNFIXED   0
/** \ref Particle::ext_flag "ext_flag" value for particle subject to an external force. */
#define PARTICLE_EXT_FORCE 1
/** \ref Particle::ext_flag "ext_flag" value for fixed particle. */
#define PARTICLE_FIXED     2
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
  /** particle type, used for non bonded interactions. */
  int    type;  

#ifdef ELECTROSTATICS
  /** charge. */
  double q;
#endif
} ParticleProperties;

/// Positional information on a particle
typedef struct {
  /** periodically folded position. */
  double p[3];
  
#ifdef ROTATION
  /** quaternions to define particle orientation */
  double quat[4]; 
#endif
} ParticlePosition;

/// Force information on a particle
typedef struct {
  /** force. */
  double f[3];

#ifdef ROTATION
  /** torque */
  double torque[3];
#endif
} ParticleForce;

/// Momentum information on a particle
typedef struct {
  /** velocity. */
  double v[3];
  
#ifdef ROTATION
  /** angular velocity */
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
  /** flag wether to fix a particle in space. Values: 
      <ul> <li> 0 no external influence
      <li> 1 apply external force \ref Particle::ext_force
      <li> 2 fix particle in space </ul> */
  int ext_flag;
  /** External force, apply if \ref Particle::ext_flag == 1. */
  double ext_force[3];
#endif

} ParticleLocal;

/// Struct holding all information for one particle.
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

  /** bonded interactions list. */
  IntList bl;
} Particle;

/** List of particles. The particle array is resized using a sophisticated
    (we hope) algorithm to avoid unnecessary resizes.
    Access using \ref realloc_particles, \ref got_particle,...
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
/** total number of particles. */
extern int  n_total_particles;

/** Capacity of the \ref particle_node / \ref local_particles. */
extern int  max_particle_node;
/** Used only on master node: particle->node mapping. */
extern int  *particle_node;
/** id->particle mapping on all nodes. */
extern Particle   **local_particles;

/** Particles' current configuration. Before using that
    call \ref updatePartCfg or \ref sortPartCfg to allocate
    the data if necessary (which is decided by \ref updatePartCfg). */
extern Particle *partCfg;
/** if non zero, \ref partCfg is sorted by particle order, and
    the particles are stored consecutively starting with 0. */
extern int partCfgSorted;


/************************************************
 * functions
 ************************************************/

/** implementation of the tcl command \ref tcl_part. This command allows to
    modify particle data. */
int part(ClientData data, Tcl_Interp *interp,
	 int argc, char **argv);

/** initialize a particle list.
 *  Use with care and ONLY for initialization! */
void init_particleList(ParticleList *pList);

/** initialize a particle.
    This function just sets all values to zero!
    Do NOT use this without setting the values of the  
    \ref ParticleProperties::identity "identity" and \ref ParticlePosition::p "position" to 
    reasonable values. Also make sure that you update \ref local_particles.

    Add here all initializations you need to be done !!!
    If dynamic memory allocation is involved, also look at \ref free_particle.
*/
void init_particle(Particle *part);

/** ... and deallocate the dynamic storage of a particle. */  
void free_particle(Particle *part);

/** allocate storage for local particles and ghosts.
    \param plist the list on which to operate
    \param size the size to provide at least. It is rounded
    up to multiples of \ref PART_INCREMENT.
    \return true iff particle adresses have changed */
int realloc_particles(ParticleList *plist, int size);

/** search for a specific particle.
    \param plist the list on which to operate 
    \param id the identity of the particle to search
    \return a pointer to the particle structure or NULL if particle is
    not in this list */
Particle *got_particle(ParticleList *plist, int id);

/** append a particle at the end of a particle List.
    reallocates particles if necessary!
    This procedure does not care for \ref local_particles.
    \param plist List to append the particle to.
    \param part  Particle to append.  
    \return Pointer to new location of the particle. */
Particle *append_unindexed_particle(ParticleList *plist, Particle *part);

/** append a particle at the end of a particle List.
    reallocates particles if necessary!
    This procedure cares for \ref local_particles.
    \param plist List to append the particle to.
    \param part  Particle to append.  
    \return Pointer to new location of the particle. */
Particle *append_indexed_particle(ParticleList *plist, Particle *part);

/** remove a particle from one particle List and append it to  another.
    Refill the destList with last particle and update its entry in local_particles. 
    reallocates particles if necessary.
    This procedure does not care for \ref local_particles.
    \param destList   List where the particle is appended.
    \param sourceList List where the particle will be removed.
    \param ind        Index of the particle in the sourceList.
    \return Pointer to new location of the particle.
 */
Particle *move_unindexed_particle(ParticleList *destList, ParticleList *sourceList, int ind);

/** remove a particle from one particle List and append it to  another.
    Refill the destList with last particle and update its entry in local_particles. 
    reallocates particles if necessary.
    This procedure cares for \ref local_particles.
    \param destList   List where the particle is appended.
    \param sourceList List where the particle will be removed.
    \param ind        Index of the particle in the sourceList.
    \return Pointer to new location of the particle.
 */
Particle *move_indexed_particle(ParticleList *destList, ParticleList *sourceList, int ind);

/** update the entries in \ref local_particles for all particles in the list pl.
    @param pl the list to put in.
*/
void update_local_particles(ParticleList *pl);

/** remove bond from particle if possible */
int try_delete_bond(Particle *part, int *bond);

/** rebuild \ref particle_node from scratch.
    After a simulation step \ref particle_node has to be rebuild
    since the particles might have gone to a different node.
*/
void build_particle_node();

/** invalidate \ref particle_node. This has to be done
    at the beginning of the integration.
*/
void particle_invalidate_part_node();

/** realloc \ref local_particles. */
void realloc_local_particles();

/** get particle data. Note that the bond intlist is
    allocated so that you are responsible to free it later.
    @param part the identity of the particle to fetch
    @param data where to store its contents.
    @return TCL_OK if particle existed
*/
int get_particle_data(int part, Particle *data);

/** Call only on the master node.
    Move a particle to a new position.
    If it does not exist, it is created.
    @param part the identity of the particle to move
    @param p    its new position
    @return TCL_OK if particle existed, TCL_CONTINUE
    if created and TCL_ERROR if id is illegal
*/
int place_particle(int part, double p[3]);

/** Call only on the master node: set particle velocity.
    @param part the particle.
    @param v its new velocity.
    @return TCL_OK if particle existed
*/
int set_particle_v(int part, double v[3]);

/** Call only on the master node: set particle force.
    @param part the particle.
    @param F its new force.
    @return TCL_OK if particle existed
*/
int set_particle_f(int part, double F[3]);

/** Call only on the master node: set particle charge.
    @param part the particle.
    @param q its new charge.
    @return TCL_OK if particle existed
*/
int set_particle_q(int part, double q);

/** Call only on the master node: set particle type.
    @param part the particle.
    @param type its new type.
    @return TCL_OK if particle existed
*/
int set_particle_type(int part, int type);

#ifdef ROTATION
/** Call only on the master node: set particle orientation using quaternions.
    @param part the particle.
    @param quat its new value for quaternions.
    @return TCL_OK if particle existed
*/
int set_particle_quat(int part, double quat[4]);

/** Call only on the master node: set particle angular velocity.
    @param part the particle.
    @param omega its new angular velocity.
    @return TCL_OK if particle existed
*/
int set_particle_omega(int part, double omega[3]);

/** Call only on the master node: set particle torque.
    @param part the particle.
    @param torque its new torque.
    @return TCL_OK if particle existed
*/
int set_particle_torque(int part, double torque[3]);
#endif

#ifdef EXTERNAL_FORCES
/** Call only on the master node: set particle external forced.
    @param part  the particle.
    @param flag  new value for ext_flag.
    @param force new value for ext_force.
    @return TCL_OK if particle existed
*/
int set_particle_ext(int part, int flag, double force[3]);
#endif

/** Call only on the master node: change particle bond.
    @param part     identity of principal atom of the bond.
    @param bond     field containing the bond type number and the
    identity of all bond partners (secundary atoms of the bond).
    @param delete   if true, do not add the bond, rather delete it if found
    @return TCL_OK on success or TCL_ERROR if no success
    (e. g. particle or bond to delete does not exist)
*/
int change_particle_bond(int part, int *bond, int delete);

/** remove particle with a given identity. Also removes all bonds to the particle.
    @param part     identity of the particle to remove
    @return TCL_OK on success or TCL_ERROR if particle does not exist
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

/** Get the complete unsorted informations on all particles into
    \ref partCfg if something's changed. This is a severe performance
    drawback and might even fail for lack of memory for large systems.
    If you need the particle info sorted, call \ref sortPartCfg instead.
    This function is lazy.
*/
void updatePartCfg();

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
*/
void local_place_particle(int part, double p[3]);

/** Used by \ref mpi_place_particle, should not be used elsewhere.
    Called if on a different node a new particle was added.
    @param part the identity of the particle added
*/
void added_particle(int part);

/** Used by \ref mpi_send_bond, should not be used elsewhere.
    Modify a bond.
    @param part the identity of the particle to change
    @param bond the bond to do
    @param delete if true, delete the bond instead of add
    @return TCL_OK for add or successful delete, TCL_ERROR else
*/
int local_change_bond(int part, int *bond, int delete);

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

/** Check the existence of a bond partner on that node and return the
    corresponding particle pointer or force exit.
    @param id particle identity.
    @return Particle pointer.
 */
MDINLINE Particle *checked_particle_ptr(int id)
{
  Particle *p = local_particles[id];
  if(!p) {
    fprintf(stderr,"ERROR: Atom %d has bond to unknown particle "
	    "(probably on different node)\n", id); 
    errexit();
  }
  return p;
}

#endif

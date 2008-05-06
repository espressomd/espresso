// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef COMM_H
#define COMM_H
/** \file communication.h
    This file contains the asynchronous MPI communication.
 
    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

    It is the header file for \ref communication.c "communication.c".

    The asynchronous MPI communication is used during the script
    evaluation. Except for the master node that interpretes the Tcl
    script, all other nodes wait in \ref mpi_loop for the master node to
    issue an action using \ref mpi_issue. \ref mpi_loop immediately
    executes an MPI_Bcast and therefore waits for the master node to broadcast
    a command, which is done by \ref mpi_issue. The request consists of
    three integers, the first one describing the action issued, the second
    and third an arbitrary parameter depending on the action issued. If applicable,
    the second parameter is the node number of the slave this request is dedicated
    to.

    Adding new actions (e. g. to implement new Tcl commands) is
    simple. First, the action has to be assigned a new action number
    (the defines like \ref REQ_BCAST_PAR) by adding a new define and
    increasing \ref REQ_MAXIMUM. Then you write a mpi_* procedure that
    does a
    \verbatim mpi_issue(request, node, param)\endverbatim
    where pnode and param are arbitrary values which will be passed to the slave
    procedure.

    After this your procedure is free to do anything. However, it has
    to be in (MPI) sync with what your new mpi_*_slave does.  This
    procedure is called immediately after the broadcast with the
    arbitrary integer as parameter.  To this aim it has also to be
    added to \ref #callbacks (the array index gives your action number.
    Last but not least for debugging purposes you can add a nice name
    to \ref #names in the same way.  */

/* from here we borrow the enumeration of
   the global variables */
#include "global.h"
#include "particle_data.h"
#include "random.h"
#include "topology.h"

/**************************************************
 * exported variables
 **************************************************/

/** \name Exported Variables */
/*@{*/
/** The number of this node. */
extern int this_node;
/** The total number of nodes. */
extern int n_nodes;
/*@}*/

/**************************************************
 * for every procedure requesting a MPI negotiation
 * a slave exists which processes this request on
 * the slave nodes. It is denoted by *_slave.
 **************************************************/

/** \name Exported Functions */
/*@{*/
/** Initialize MPI and determine \ref n_nodes and \ref this_node. */
void mpi_init(int *argc, char ***argv);

/** Process requests from master node. Slave nodes main loop. */
void mpi_loop();

/** Issue REQ_TERM: stop Espresso, all slave nodes exit. */
void mpi_stop();

/** Finalize MPI. Called by all nodes upon exit */
void mpi_finalize();

/** Issue REQ_BCAST_PAR: broadcast a parameter from datafield.
    @param i the number from \ref global.h "global.h" referencing the datafield.
    @return nonzero on error
*/
int mpi_bcast_parameter(int i);

/** Issue REQ_WHO_HAS: ask nodes for their attached particles. */
void mpi_who_has();

/** Issue REQ_EVENT: tells all clients of some system change.
    The events are:
    <ul>
    <li> PARTICLE_CHANGED
    <li> INTERACTION_CHANGED
    </ul>
    Then all nodes execute the respective on_* procedure from \ref initialize.c
    Note that not all of these codes are used. Since some actions (like placing a
    particle) include communication anyways, this is handled by the way.
*/
void mpi_bcast_event(int event);

/** Issue REQ_PLACE: move particle to a position on a node.
    Also calls \ref on_particle_change.
    \param id   the particle to move.
    \param new  if non-zero, the particle is new
    \param node the node to attach it to.
    \param pos  the particles position.
*/
void mpi_place_particle(int node, int id, int new, double pos[3]);

/** Issue REQ_SET_V: send particle velocity.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param v its new velocity.
*/
void mpi_send_v(int node, int part, double v[3]);

/** Issue REQ_SET_F: send particle force.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param F its new force.
*/
void mpi_send_f(int node, int part, double F[3]);

/** Issue REQ_SET_M: send particle mass.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param mass its new mass.
*/
void mpi_send_mass(int node, int part, double mass);

/** Issue REQ_SET_Q: send particle charge.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param q its new charge.
*/
void mpi_send_q(int node, int part, double q);

#ifdef ROTATION
/** Issue REQ_SET_QUAT: send particle orientation.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param quat its new quaternions.
*/
void mpi_send_quat(int node, int part, double quat[4]);

/** Issue REQ_SET_LAMBDA: send particle angular velocity.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param omega its new angular velocity.
*/
void mpi_send_omega(int node, int part, double omega[3]);

/** Issue REQ_SET_TORQUE: send particle torque.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param torque its new torque.
*/
void mpi_send_torque(int node, int part, double torque[3]);
#endif

#ifdef DIPOLES 
/** Issue REQ_SET_DIP: send particle dipole orientation.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param dip its new dipole orientation.
*/
void mpi_send_dip(int node, int part, double dip[3]);
/** Issue REQ_SET_DIPM: send particle dipole moment.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param dipm its new dipole moment (absolut value).
*/
void mpi_send_dipm(int node, int part, double dipm);
#endif

#ifdef VIRTUAL_SITES
/** Issue REQ_SET_DIPM: send particle dipole moment.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param isVirtual its new isVirtual.
*/
void mpi_send_isVirtual(int node, int part, int isVirtual);
#endif

/** Issue REQ_SET_TYPE: send particle type.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param type its new type.
*/
void mpi_send_type(int node, int part, int type);

/** Issue REQ_SET_MOL_ID: send molecule id.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param mid its new mol_id.
*/
void mpi_send_mol_id(int node, int part, int mid);

/** Issue REQ_SET_BOND: send bond.
    Also calls \ref on_particle_change.
    \param pnode    node it is attached to.
    \param part     identity of principal atom of the bond.
    \param bond     field containing the bond type number and the identity of all bond partners (secundary atoms of the bond).
    \param delete   if true, do not add the bond, rather delete it if found
    \return 1 on success or 0 if not (e. g. bond to delete does not exist)
*/
int mpi_send_bond(int pnode, int part, int *bond, int delete);

/** Issue REQ_SET_EXCLUSION: send exclusions.
    Also calls \ref on_particle_change.
    \param part     identity of first particle of the exclusion.
    \param part2    identity of secnd particle of the exclusion.
    \param delete   if true, do not add the exclusion, rather delete it if found
*/
void mpi_send_exclusion(int part, int part2, int delete);

/** Issue REQ_REM_PART: remove a particle.
    Also calls \ref on_particle_change.
    \param id   the particle to remove.
    \param node the node it is attached to.
*/
void mpi_remove_particle(int node, int id);

/** Issue REQ_GET_PART: recv particle data. The data has to be freed later
    using \ref free_particle, otherwise the dynamically allocated parts, bonds
    and exclusions are left over.
    \param part the particle.
    \param node the node it is attached to.
    \param part_data where to store the received data.
    \note Gets a copy of the particle data not a pointer to the actual particle
    used in integration
*/
void mpi_recv_part(int node, int part, Particle *part_data);

/** Issue REQ_INTEGRATE: start integrator.
    @param n_steps how many steps to do.
    @return nonzero on error
*/
int mpi_integrate(int n_steps);

/** Issue REQ_BCAST_IA: send new ia params.
    Also calls \ref on_short_range_ia_change.

    mpi_bcast_ia_params is used for both, bonded and non-bonded
    interaction parameters. Therefor i and j are used depending on
    their value:

    \param i   particle type for non bonded interaction parameters / 
               bonded interaction type number.
    \param j   if not negative: particle type for non bonded interaction parameters /
               if negative: flag for bonded interaction */
void mpi_bcast_ia_params(int i, int j);

/** Issue REQ_BCAST_IA_SIZE: send new size of \ref ia_params.
    \param s the new size for \ref ia_params.
*/
void mpi_bcast_n_particle_types(int s);

/** Issue REQ_GATHER: gather data for analysis in \ref #analyze.
    \param job what to do:
    <ul>
	<li> 1 calculate and reduce (sum up) energies, using \ref energy_calc.
	<li> 2 calculate and reduce (sum up) pressure, stress tensor, using \ref pressure_calc.
	<li> 3 calculate and reduce (sum up) instantaneous pressure, using \ref pressure_calc.
    </ul>
    \param result where to store the gathered value(s):
    <ul><li> job=1 unused (the results are stored in a global 
	     energy array of type \ref Observable_stat)
	<li> job=2 unused (the results are stored in a global 
	     virials array of type \ref Observable_stat)
	<li> job=3 unused (the results are stored in a global 
	     virials array of type \ref Observable_stat)
    \param result_t where to store the gathered value(s):
    <ul><li> job=1 unused (the results are stored in a global 
	     energy array of type \ref Observable_stat)
	<li> job=2 unused (the results are stored in a global 
	     p_tensor tensor of type \ref Observable_stat)
	<li> job=3 unused (the results are stored in a global 
	     p_tensor tensor of type \ref Observable_stat)
    \param result_nb where to store the gathered value(s):
    <ul><li> job=1 unused (the results are stored in a global 
	     energy array of type \ref Observable_stat_non_bonded)
	<li> job=2 unused (the results are stored in a global 
	     virials_non_bonded array of type \ref Observable_stat_non_bonded)
	<li> job=3 unused (the results are stored in a global 
	     virials_non_bonded array of type \ref Observable_stat_non_bonded)
    \param result_t_nb where to store the gathered value(s):
    <ul><li> job=1 unused (the results are stored in a global 
	     energy array of type \ref Observable_stat_non_bonded)
	<li> job=2 unused (the results are stored in a global 
	     p_tensor_non_bonded tensor of type \ref Observable_stat_non_bonded)
	<li> job=3 unused (the results are stored in a global 
	     p_tensor_non_bonded tensor of type \ref Observable_stat_non_bonded)
    </ul>
*/
void mpi_gather_stats(int job, void *result, void *result_t, void *result_nb, void *result_t_nb);

/** Issue GET_LOCAL_STRESS_TENSOR: gather the contribution to the local stress tensors from
    each node.
 */

void mpi_local_stress_tensor(DoubleList *TensorInBin, int bins[3], int periodic[3], double range_start[3], double range[3]);

/** Issue REQ_GETPARTS: gather all particle informations (except bonds).
    This is slow and may use huge amounts of memory. If il is non-NULL, also
    the bonding information is also fetched and stored in a single intlist
    pointed to by il. The particles bonding information references this array,
    which is the only data you have to free later (besides the result array
    you allocated). YOU MUST NOT CALL \ref free_particle on any of these particles!
    
  \param result where to store the gathered particles
  \param il if non-NULL, the integerlist where to store the bonding info
*/
void mpi_get_particles(Particle *result, IntList *il);

/** Issue REQ_SET_TIME_STEP: send new \ref time_step and rescale the
    velocities accordingly. 
*/
void mpi_set_time_step(double time_step);

/** Issue REQ_BCAST_COULOMB: send new coulomb parameters. */
void mpi_bcast_coulomb_params();

/** Issue REQ_SEND_EXT: send nex external flag and external force. */
void mpi_send_ext(int pnode, int part, int flag, int mask, double force[3]);

/** Issue REQ_BCAST_COULOMB: send new coulomb parameters. */
void mpi_bcast_constraint(int del_num);

/** Issue REQ_RANDOM_SEED: read/set seed of random number generators on each node. */
void mpi_random_seed(int cnt, long *seed);

/** Issue REQ_RANDOM_STAT: read/set status of random number generators on each node. */
void mpi_random_stat(int cnt, RandomStatus *stat);

/** Issue REQ_BCAST_LJFORCECAP: initialize LJ force capping. */
void mpi_lj_cap_forces(double force_cap);

/** Issue REQ_BCAST_MORSEFORCECAP: initialize Morse force capping. */
void mpi_morse_cap_forces(double force_cap);

/** Issue REQ_BCAST_BUCKFORCECAP: initialize Buckingham force capping. */
void mpi_buck_cap_forces(double force_cap);

/** Issue REQ_BCAST_TABFORCECAP: initialize tabulated force capping. */
void mpi_tab_cap_forces(double force_cap);

/** Issue REQ_GET_CONSFOR: get force acting on constraint */
void mpi_get_constraint_force(int constraint, double force[3]);

/** Issue REQ_BIT_RANDOM_SEED: read/set seed of the bit random number generators on each node. */
void mpi_bit_random_seed(int cnt, int *seed);

/** Issue REQ_BIT_RANDOM_STAT: read/set status of the bit random number generators on each node. */
void mpi_bit_random_stat(int cnt, BitRandomStatus *stat);

/** Issue REQ_RESCALE_PART: rescales all particle positions in direction 'dir' by a factor 'scale'. */
void mpi_rescale_particles(int dir, double scale);

/** Issue REQ_BCAST_CS: change the cell structure on all nodes. */
void mpi_bcast_cell_structure(int cs);

/** Issue REQ_BCAST_NPTISO_GEOM: broadcast nptiso geometry parameter to all nodes. */
void mpi_bcast_nptiso_geom(void);

/** Issue REQ_BCAST_LJANGLEFORCECAP: initialize LJANGLE force capping. */
void mpi_ljangle_cap_forces(double force_cap);


/** Issue REQ_UPDATE_MOL_IDS: Update the molecule ids so that they are
    in sync with the topology.  Note that this only makes sense if you
    have a simple topology such that each particle can only belong to
    a single molecule */
void mpi_update_mol_ids(void); 

/** Issue REQ_SYNC_TOPO: Update the molecules ids to that they correspond to the topology */
int mpi_sync_topo_part_info(void);

/** Issue REQ_BCAST_LBPAR: Broadcast a parameter for Lattice Boltzmann.
 * @param field References the parameter field to be broadcasted. The references are defined in \ref lb.h "lb.h"
 */
void mpi_bcast_lb_params(int field);

/** Issue REQ_SEND_FLUID: Send a single lattice site to a processor.
 * @param node  processor to send to
 * @param index index of the lattice site
 * @param rho   local fluid density
 * @param j     local fluid velocity
 */
void mpi_send_fluid(int node, int index, double rho, double *j, double *pi);

/** Issue REQ_GET_FLUID: Receive a single lattice site from a processor.
 * @param node  processor to send to
 * @param index index of the lattice site
 * @param rho   local fluid density
 * @param j     local fluid velocity
 */
void mpi_recv_fluid(int node, int index, double *rho, double *j, double *pi);

/** Issue REQ_GET_ERRS: gather all error messages from all nodes and set the interpreter result
    to these error messages. This called only on the master node.
    The errors are append to the result, if ret_state == TCL_ERROR, otherwise the result is overwritten.
    Therefore you should end any Tcl command handler by return gather_runtime_errors(<return_value>).
    This code uses asynchronous communication.
    @param ret_state return value of the procedure
    @param interp where to put the errors
    @return new return value after the background errors, if any, have been handled
*/
int mpi_gather_runtime_errors(Tcl_Interp *interp, int ret_state);

/*@}*/

/** \name Event codes for \ref mpi_bcast_event
    These codes are used by \ref mpi_bcast_event to notify certain changes
    of doing something now.
*/
/*@{*/
#define P3M_COUNT_CHARGES 0
#define INVALIDATE_SYSTEM 1
#define CHECK_PARTICLES   2
#define MAGGS_COUNT_CHARGES 3
#define EWALD_COUNT_CHARGES 4
/*@}*/

#endif

#ifndef COMM_H
#define COMM_H
/** \file communication.h
    This file contains the random access MPI communication.
 
    <b>Responsible:</b>
    <a href="mailto:arnolda@mpip-mainz.mpg.de">Axel</a>

    It is the header file for \ref communication.c "communication.c".

    The random access MPI communication is used during script
    evaluation. Except for the master node that interpretes the Tcl
    script, all other nodes wait for the master node to issue an
    action.  That is done by \ref mpi_loop, the main loop for all
    slave nodes. It simply does an MPI_Bcast and therefore waits for
    the master node to broadcast a command. A command consists of two
    integers, the first one describing the action issued, the second
    an arbitrary parameter depending on the action issued. The actions
    are always issued by the master node, through the functions
    exported below (e. g. \ref mpi_bcast_parameter).

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

/** Issue REQ_BCAST_PAR: broadcast a parameter from datafield.
    \param i the number from \ref global.h "global.h" referencing the datafield. */
void mpi_bcast_parameter(int i);

/** Issue REQ_WHO_HAS: ask nodes for their attached particles. */
void mpi_who_has();

/** Issue REQ_EVENT: tells all clients of some system change.
    The events are:
    <ul>
    <li> PARTICLE_CHANGED
    <li> INTERACTION_CHANGED
    <li> PARAMETER_CHANGED
    <li> TOPOLOGY_CHANGED
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

/** Issue REQ_SET_Q: send particle charge.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param q its new charge.
*/
void mpi_send_q(int node, int part, double q);

/** Issue REQ_SET_TYPE: send particle type.
    Also calls \ref on_particle_change.
    \param part the particle.
    \param node the node it is attached to.
    \param type its new type.
*/
void mpi_send_type(int node, int part, int type);

/** Issue REQ_SET_BOND: send bond.
    Also calls \ref on_particle_change.
    \param pnode    node it is attached to.
    \param part     identity of principal atom of the bond.
    \param bond     field containing the bond type number and the identity of all bond partners (secundary atoms of the bond).
    \param delete   if true, do not add the bond, rather delete it if found
    \return 1 on success or 0 if not (e. g. bond to delete does not exist)
*/
int mpi_send_bond(int pnode, int part, int *bond, int delete);

/** Issue REQ_REM_PART: remove a particle.
    Also calls \ref on_particle_change.
    \param id   the particle to remove.
    \param node the node it is attached to.
*/
void mpi_remove_particle(int node, int id);

/** Issue REQ_GET_PART: recv particle data.
    \param part the particle.
    \param node the node it is attached to.
    \param part_data where to store the received data.
*/
void mpi_recv_part(int node, int part, Particle *part_data);

/** Issue REQ_INTEGRATE: start integrator.
    \param n_steps how many steps to do.
*/
void mpi_integrate(int n_steps);

/** Issue REQ_BCAST_IA: send new ia params.
    Also calls \ref on_ia_change.

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

/** Issue REQ_GATHER: gather data for analysis in \ref analyze.
    \param job what to do:
    <ul><li> 0 gather \ref minimum_part_dist
	<li> 1 calculate and reduce (sum up) energies, using \ref calc_energy.
	<li> 2 calculate and reduce (sum up) virials, using \ref calc_virials.
    </ul>
    \param result where to store the gathered value(s):
    <ul><li> job=0 a double *
             usage: double *buf; mpi_gather_stats(0, buf);
	<li> job=1 unused (the results are stored in a global 
	     energy array of type \ref Energy_stat)
	<li> job=2 unused (the results are stored in a global 
	     virials array of type \ref Energy_stat)
    </ul>
*/
void mpi_gather_stats(int job, void *result);

/** Issue REQ_GETPARTS: gather all particle informations (except bonds).
    This is slow and may use huge amounts of memory. If il is non-NULL, also
    the bonding information is also fetched and stored in a single intlist
    pointed to by il. The particles bonding information references this array,
    which is the only data you have to free later (besides the result array
    you allocated).
  \param result where to store the gathered particles
  \param il if non-NULL, the integerlist where to store the bonding info
*/
void mpi_get_particles(Particle *result, IntList *il);

/** Issue REQ_SET_TIME_STEP: send new \ref time_step and rescale the
    velocities accordingly. 
*/
void mpi_set_time_step();

/** Issue REQ_BCAST_COULOMB: send new coulomb parameters. */
void mpi_bcast_coulomb_params();

/** Issue REQ_SEND_EXT: send nex external flag and external force. */
void mpi_send_ext(int pnode, int part, int flag, double force[3]);

/** Issue REQ_BCAST_COULOMB: send new coulomb parameters. */
void mpi_bcast_constraint(int del_num);

/** Issue REQ_RANDOM_SEED: read/set seed of random number generators on each node. */
void mpi_random_seed(int cnt, long *seed);

/** Issue REQ_RANDOM_STAT: read/set status of random number generators on each node. */
void mpi_random_stat(int cnt, RandomStatus *stat);

/** Issue REQ_BCAST_LJFORCECAP: initialize LJ force capping. */
void mpi_lj_cap_forces(double force_cap);

/*@}*/

/** \name Event codes for \ref mpi_bcast_event
    These codes are used by \ref mpi_bcast_event to notify certain changes
    of the systems state to all nodes.
*/
/*@{*/
#define PARTICLE_CHANGED 0
#define INTERACTION_CHANGED 1
#define PARAMETER_CHANGED 2
#define TOPOLOGY_CHANGED 3
#define P3M_COUNT_CHARGES 4
/*@}*/

#endif

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

    \verbatim MPI_Bcast(req, 2, MPI_INT, 0, MPI_COMM_WORLD)\endverbatim 

    where 

    \verbatim req[0]=<action number>; req[1]=<arbitrary>\endverbatim. 

    After this your procedure is free to do anything. However, it has
    to be in (MPI) sync with what your new mpi_*_slave does.  This
    procedure is called immediately after the broadcast with the
    arbitrary integer as parameter.  To this aim it has also to be
    added to \ref callbacks (the array index gives your action number
    - better not choose it too high...). Last but not least for
    debugging purposes you can add a nice name to \ref names in the
    same way.  */

/* from here we borrow the enumeration of
   the global variables */
#include "global.h"
#include "particle_data.h"

/**************************************************
 * data types
 **************************************************/

typedef struct {
  int n_particles;
  float *coords;
} float_packed_particle_data;

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

/** Issue REQ_TERM: stop tcl_md, all slave nodes exit. */
void mpi_stop();

/** Issue REQ_BCAST_PAR: broadcast a parameter from datafield.
    \param i the number from \ref global.h "global.h" referencing the datafield. */
void mpi_bcast_parameter(int i);

/** Issue REQ_WHO_HAS: ask nodes for their attached particles. */
void mpi_who_has();

/** Issue REQ_CHTOPL: call this if you change the topology (box, grid, etc.).
 */
void mpi_changed_topology();

/** Issue REQ_PLACE: move particle to a position on a node.
    \param id   the particle to move.
    \param node the node to attach it to.
    \param pos  the particles position.
*/
void mpi_place_particle(int node, int id, double pos[3]);

/** Issue REQ_SET_V: send particle velocity.
    \param part the particle.
    \param node the node it is attached to.
    \param v its new velocity.
*/
void mpi_send_v(int node, int part, double v[3]);

/** Issue REQ_SET_F: send particle force.
    \param part the particle.
    \param node the node it is attached to.
    \param F its new force.
*/
void mpi_send_f(int node, int part, double F[3]);

/** Issue REQ_SET_Q: send particle charge.
    \param part the particle.
    \param node the node it is attached to.
    \param q its new charge.
*/
void mpi_send_q(int node, int part, double q);

/** Issue REQ_SET_TYPE: send particle type.
    \param part the particle.
    \param node the node it is attached to.
    \param type its new type.
*/
void mpi_send_type(int node, int part, int type);

/** Issue REQ_SET_BOND: send bond.
    \param part     identity of principal atom of the bond.
    \param node     node it is attached to.
    \param bond     field containing the bond type number and the identity of all bond partners (secundary atoms of the bond).
    \param delete   if true, do not add the bond, rather delete it if found
    \return 1 on success or 0 if not (e. g. bond to delete does not exist)
*/
int mpi_send_bond(int pnode, int part, int *bond, int delete);


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

/** Issue REQ_GATHER: gather statistics. job determines the job to
    do, at the moment either gather \ref minimum_part_dist or
    the coordinates.
    \param job what to do:
    <ul>
    <li> 0 gather \ref minimum_part_dist
    <li> 1 gather coordinates in packed float format (imd)
    </ul>
    \param result where to store the gathered value(s):
    <ul>
    <li> for job 0 a double *
    <li> for job 1 a float_packed_particle_data *
    </ul>
*/
void mpi_gather_stats(int job, void *result);

/*@}*/

#endif

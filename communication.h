#ifndef COMM_H
#define COMM_H
/** \file communication.h
    MPI communication. Header file for \ref communication.c "communication.c".
*/

/* from here we borrow the enumeration of
   the global variables */
#include "global.h"
#include "particle_data.h"

/**************************************************
 * exported variables
 **************************************************/

/** \name Exported Variables */
/*@{*/
/** the number of this node */
extern int this_node;
/** the total number of nodes */
extern int nprocs;
/*@}*/

/**************************************************
 * for every procedure requesting a MPI negotiation
 * a slave exists which processes this request on
 * the slave nodes. It is denoted by *_slave.
 **************************************************/

/** \name Exported Functions */
/*@{*/
/** Initialize MPI and determine \ref nprocs and \ref this_node. */
void mpi_init(int *argc, char ***argv);

/** Process requests from master node. Slave nodes main loop. */
void mpi_loop();

/** Issue REQ_TERM: stop MPI and determine nprocs/node. */
void mpi_stop();

/** Issue REQ_BCAST_PARM: broadcast a parameter from datafield.
    \param i the number from \ref global.h "global.h" referencing the datafield. */
void mpi_bcast_parameter(int i);

/** Issue REQ_WHO_HAS: ask nodes for their attached particles. */
void mpi_who_has();

/** Issue REQ_ATTACH: move particle to a node.
    \param part the particle to move.
    \param node the node to attach it to.
*/
void mpi_attach_particle(int part, int node);

/** Issue REQ_SET_POS: send particle position.
    \param part the particle.
    \param node the node it is attached to.
    \param pos its new position.
*/
void mpi_send_pos(int node, int part, double pos[3]);

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
    \param i,j the particle types whose parameters are to be sent.
*/
void mpi_bcast_ia_params(int i, int j);
/*@}*/

#endif

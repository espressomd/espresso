#ifndef COMM_H
#define COMM_H
#include "global.h"

/**************************************************
 * for every procedure requesting a MPI negotiation
 * a slave exists which processes this request on
 * the slave nodes. It is denoted by *_slave.
 **************************************************/

/** initialize MPI and determine nprocs/node */
void mpi_init(int *argc, char ***argv);

/** process requests from master node */
void mpi_loop();

/** REQ_TERM: stop MPI and determine nprocs/node */
void mpi_stop();
void mpi_stop_slave(int parm);

/** REQ_WHO_HAS: ask nodes for particle */
int mpi_who_has(int particle);
void mpi_who_has_slave(int parm);

/** REQ_ATTACH: move particle to a node */
void mpi_attach_particle(int part, int node);
void mpi_attach_particle_slave(int parm);

/** REQ_SET_POS: send particle position */
void mpi_send_pos(int node, int part, double pos[3]);
void mpi_send_pos_slave(int parm);

/** REQ_BCAST_PARM: broadcast a parameter from datafield */
void mpi_bcast_parameter(int i);
void mpi_bcast_parameter_slave(int parm);

/** REQ_GET_PART: recv particle data */
void mpi_recv_part(int node, int part, Particle *part_data);
void mpi_recv_part_slave(int parm);

/** REQ_INTEGRATE: start integrator */
void mpi_integrate(int n_steps);
void mpi_integrate_slave(int parm);

/*************************************
 * requests in random access mode
 * KEEP ORDER IN SYNC WITH callbacks[]
 * ALSO ADD FOR EVERY REQUEST
 * THE CALLBACK IN callbacks[]
 *************************************/
/* slave callback procedure */
typedef void (SlaveCallback)(int param);

#define REQ_TERM      0
#define REQ_BCAST_PAR 1
#define REQ_WHO_HAS   2
#define REQ_ATTACH    3
#define REQ_SET_POS   4
#define REQ_SET_V     5
#define REQ_SET_F     6
#define REQ_ADD_B2    7
#define REQ_ADD_B3    8
#define REQ_GET_PART  9
#define REQ_INTEGRATE 10
#define REQ_MAXIMUM   11

extern SlaveCallback *callbacks[];

#endif

#ifndef COMM_H
#define COMM_H
#include "global.h"

/** initialize MPI and determine nprocs/node */
void init_mpi(int *argc, char ***argv);

/** stop MPI and determine nprocs/node */
void stop_mpi();

/** process requests from master node */
void mpi_loop();

/** ask nodes for particle */
int mpi_who_has(int particle);

/** move particle to a node */
void mpi_attach_particle(int part, int node);

/** send particle data */
void mpi_send_pos(int node, int part, double pos[3]);

/* recv particle data */
void mpi_recv_part(int node, int part, Particle *part_data);

/* requests in random access mode */
#define REQ_TERM      0
#define REQ_WHO_HAS   1
#define REQ_ATTACH    2
#define REQ_SET_POS   3
#define REQ_SET_V     4
#define REQ_SET_F     5
#define REQ_ADD_B2    6
#define REQ_ADD_B3    7
#define REQ_GET_PART  8

#endif

/**************************************************/
/*******************  GHOSTS.H  *******************/
/**************************************************/
#ifndef GHOSTS_H 
#define GHOSTS_H 

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "global.h"
#include "cells.h"
#include "communication.h"
#include "grid.h"

#define PART_INCREMENT 20

/* Attention: since the ghost particle structures 
   make use of the linked cell structure, ghost_init() 
   has to be called after cells_init */

/* A word about the diretions:
   Notation X = {0,1,2,3,4,5}
   left  0    right 1
   down  2    up    3
   for   4    back  5            */

/** Structure to hold ghost particle information. 
 *  If this changes also the functions pack_gost/unpack_ghost 
 *  have to be reweritten. */
typedef struct {
  int identity;
  int type;
  double p[3];
  double q;
} Ghost;

/* Exchange Particles */

/** Buffer for particles to send. */
extern int       n_p_send_buf;
extern int       max_p_send_buf;
extern Particle *p_send_buf;
/** Buffer for particles to recieve. */
extern int       n_p_recv_buf;
extern int       max_p_recv_buf;
extern Particle *p_recv_buf;
/** Buffer for particles bonds to send. */
extern int       n_b_send_buf;
extern int       max_b_send_buf;
extern int      *b_send_buf;
/** Buffer for particles bonds recieve. */
extern int       n_b_recv_buf;
extern int       max_b_recv_buf;
extern int      *b_recv_buf;

/* exchange Ghosts */

/** maximal number of cells to send. */
extern int max_send_cells;
/** number of cells to send in direction X. */
extern int n_send_cells[6];
/** list of cell indices to send. */
extern int *send_cells;
/** number of cells to receive from direction X. */
extern int n_recv_cells[6];
/** list of cell indices to receive. */
extern int *recv_cells;
/** start indices for cells to send/recv in/from direction X. */ 
extern int cell_start[6];

/** Number of ghosts in each send cell. */ 
extern int *n_send_ghosts;
/** Number of ghosts in each recv cell. */ 
extern int *n_recv_ghosts;

/** Buffer for Ghosts to send. */
extern int   n_g_send_buf;
extern int   max_g_send_buf;
extern Ghost *g_send_buf;
/** Buffer for Ghosts to recieve. */
extern int   n_g_recv_buf;
extern int   max_g_recv_buf;
extern Ghost *g_recv_buf;


/** number of ghosts to send in direction X */
extern int ghost_send_size[6];
/** number of ghosts to recv from direction X */
extern int ghost_recv_size[6];

/** Buffer for forces/coordinates to send. */
extern double *send_buf;
extern int max_send_buf;
/** Buffer for forces/coordinates to recieve. */
extern double *recv_buf;
extern int max_recv_buf;

/** initialize ghost particle structures. */
void ghost_init();
/** exchange particles. 
    For a new verlet list setup first all local particles which have left the local box have to be send to the right processor. Procedur: direction_loop: 1) 
*/
void exchange_part();
/** exchange ghost particles. */
void exchange_ghost();
/** exchange ghost particle positions. */
void update_ghost_pos();
/** exchange ghost particle forces. */
void collect_ghost_forces();
/** exit ghost structures. */
void ghost_exit();

#endif

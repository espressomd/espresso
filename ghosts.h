/**************************************************/
/*******************  GHOSTS.H  *******************/
/**************************************************/
#ifndef GHOSTS_H 
#define GHOSTS_H 

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "cells.h"
#include "communication.h"

#define PART_INCREMENT 20

/* Attention: since the ghost particle structures 
   make use of the linked cell structure, ghost_init() 
   has to be called after cells_init */

/* A word about the diretions:
   Notation X = {1,2,3,4,5,6}
   rigth 1    left  2
   up    3    down  4
   back  5    for   6            */

/** number of cells to send in direction X. */
int n_send_cells[6];
/** list of cell indices to send. */
int *send_cells;
/** number of cells to receive from direction X. */
int n_recv_cells[6];
/** list of cell indices to receive. */
int *recv_cells;
/** start indices for cells to send/recv in/from direction X. */ 
int cell_start[6];

/** Number of ghosts in each send cell. */ 
int *n_send_ghosts;
/** Number of ghosts in each recv cell. */ 
int *n_recv_ghosts;

/** Size of send/recv buffers */
int buf_size;
/** Buffer for particles to send. */
Particle *part_send_buf;
/** Buffer for particles to recieve. */
Particle *part_recv_buf;
/** Buffer for particles to send. */
Particle *part_send_buf2;
/** Buffer for particles to recieve. */
Particle *part_recv_buf2;
/** Buffer for forces/coordinates to send. */
double *send_buf;
/** Buffer for forces/coordinates to recieve. */
double *recv_buf;

/** initialize ghost particle structures. */
void ghost_init();
/** exchange particles. 
    For a new verlet list setup first all local particles which have left the local box have to be send to the right processor. Procedur: direction_loop: 1) 
*/
void exchange_part();
/** exchange ghost particles. */
void exchange_ghost();
/** exchange ghost particle positions. */
void exchange_ghost_pos();
/** exchange ghost particle forces. */
void exchange_ghost_forces();
/** exit ghost structures. */
void ghost_exit();

#endif

/**************************************************/
/*******************  GHOSTS.H  *******************/
/**************************************************/
#ifndef GHOSTS_H 
#define GHOSTS_H 

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>

#include "global.h"
#include "communication.h"

#define PART_INCREMENT 100

/* Attention: since the ghost particle structures 
   make use of the linked cell structure, ghost_init() 
   has to be called after cells_init */

/* A word about the diretions:
   Notation X = {1,2,3,4,5,6}
   rigth 1    left  2
   up    3    down  4
   for   5    back  6            */

/** number of cells to send in direction X. */
int n_send_cells[6];
/** list of cell indices to send in direction X. */
int *send_cells[6];
/** number of cells to receive from direction X. */
int n_recv_cells[6];
/** list of cell indices to receive from direction X. */
int *recv_cells[6];

/** Buffer for particles to send. */
Particle *part_send_buf;
/** Buffer for particles to recieve. */
Particle *part_recv_buf;

/** Buffer for forces/coordinates to send. */
double *send_buf;
/** Buffer for forces/coordinates to recieve. */
double *recv_buf;

/** initialize ghost particle structures. */
void ghost_init();
/** exchange particles. */
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

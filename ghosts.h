// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef GHOSTS_H 
#define GHOSTS_H 
/** \file ghosts.h    Ghost particles and particle exchange.
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
 *
 *  In this file you find everything concerning the exchange of
 *  particle data (particles, ghosts, positions and forces) for short
 *  range interactions between the spacial domains of neighbouring
 *  nodes.
 *
 *  All structures are initialized by \ref ghost_init.
 *
 *
 */
#include <mpi.h>
#include "particle_data.h"

/** \name Defines */
/************************************************************/
/*@{*/

#define GHOST_SEND 0
#define GHOST_RECV 1
#define GHOST_BCST 2
#define GHOST_RDCE 3
#define GHOST_LOCL 4
#define GHOST_JOBMASK  15
#define GHOST_PREFETCH 16

/*@}*/

/** \name Transfer data classes */
/************************************************************/
/*@{*/

#define GHOSTTRANS_PROPRTS  1
#define GHOSTTRANS_POSITION 2
/** must be or'd together with \ref GHOSTTRANS_POSITION */
#define GHOSTTRANS_POSSHFTD 4
#define GHOSTTRANS_MOMENTUM 8
#define GHOSTTRANS_FORCE    16
#define GHOSTTRANS_PARTNUM  32

/*@}*/

/** \name Data Types */
/************************************************************/
/*@{*/

typedef struct {

  /** Communication type. */
  int type;
  /** Node to communicate with (to use with GHOST_SEND, GHOST_RECV). */
  int node;
  /** MPI communicator handle (to use with GHOST_BCST, GHOST_GATH, GHOST_RDCE). */
  MPI_Comm mpi_comm;

  /** Number of particle lists to communicate. */
  int n_part_lists;
  /** Pointer array to particle lists to communicate. */
  ParticleList **part_lists;

  /** if \ref GhostCommunicator::data_parts has \ref GHOSTTRANS_POSSHFTD, then this is the shift vector.
      Normally this a integer multiple of the box length. The shift is done on the sender side */
  double shift[3];
} GhostCommunication;

/** Properties for a ghost communication. A ghost communication is defined */
typedef struct {

  /** Particle data parts to transfer */
  int data_parts;

  /** number of communication steps. */
  int num;

  /** List of ghost communications. */
  GhostCommunication *comm;

} GhostCommunicator;

/*@}*/

/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize a communicator. */
void prepare_comm(GhostCommunicator *comm, int data_parts, int num);

/** Free a communicator. */
void free_comm(GhostCommunicator *comm);

/** Initialize ghosts. */
void ghost_init();

/** do a ghost communication */
void ghost_communicator(GhostCommunicator *gc);

/** Go through \ref ghost_cells and remove the ghost entries from \ref
    local_particles. Part of \ref dd_exchange_and_sort_particles.*/
void invalidate_ghosts();

/*@}*/

#endif

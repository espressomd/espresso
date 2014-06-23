/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
#ifndef _CELLS_H
#define _CELLS_H
/** \file cells.hpp
    This file contains everything related to the cell structure / cell
    system.
    
    The cell system (\ref Cell Structure) describes how particles are
    distributed on the cells and how particles of different cells
    (regardless if they reside on the same or different nodes)
    interact with each other. The following cell systems are implemented:
  
    <ul>
    <li> domain decomposition: The simulation box is divided spatially
    ino cells (see \ref domain_decomposition.hpp). This is suitable for
    short range interctions.
    <li> nsquare: The particles are distributed equally on all nodes
    regardless their spatial position (see \ref nsquare.hpp). This is
    suitable for long range interactions that can not be treated by a
    special method like P3M (see \ref p3m.hpp).
    <li> layered: in x and y directions, it uses a nsquared type of interaction calculation,
                  but in z it has a domain decomposition into layers.
    </ul>
  
    One can switch between different cell systems with the tcl command
    cellsystem implemented in \ref cells.cpp .
  
    Some structures are common to all cell systems: 
  
   <ul>
   <li> All cells, real cells as well as ghost cells, are stored in the array \ref cells::cells with size \ref
   n_cells. The size of this array has to be changed with \ref realloc_cells.
   <li> Their are two lists of cell pointers to acces particles and
   ghost particles on a node: \ref local_cells contains pointers to
   all cells containing the particles physically residing on that
   node. \ref ghost_cells contains pointers to all cells containing
   the ghost particles of that node. The size of these lists has to be
   changed with \ref realloc_cellplist
   <li> An example using the cell pointer lists to access particle data
   can be found in the function \ref
   print_local_particle_positions. DO NOT INVENT YOUR OWN WAY!!!
   </ul>
*/

#include "particle_data.hpp"
#include "ghosts.hpp"
#include "verlet.hpp"

/** \name Cell Structure */
/** Flag telling which cell structure is used at the moment. */
/*@{*/
/** Flag indicating that there is no cell system yet. Only at the
    VERY beginning of the code startup. */
#define CELL_STRUCTURE_NONEYET -1
/** Flag indicating that the current cell structure will be used furthor on */
#define CELL_STRUCTURE_CURRENT 0
/** cell structure domain decomposition */
#define CELL_STRUCTURE_DOMDEC  1
/** cell structure n square */
#define CELL_STRUCTURE_NSQUARE 2
/** cell structure layered */
#define CELL_STRUCTURE_LAYERED 3
/*@}*/

/** \name Flags for exchange_and_sort_particles: wether to do a global exchange
    or assume that particles did not move much (faster, used during integration,
    where moving far is a catastrophe anyways). */
/*@{*/

/** Flag for exchange_and_sort_particles : Do global exchange. */
#define CELL_GLOBAL_EXCHANGE 1
/** Flag for exchange_and_sort_particles : Do neighbor exchange. */
#define CELL_NEIGHBOR_EXCHANGE 0

/** \name Flags for cells_on_geometry_change */
/*@{*/

/** Flag for cells_on_geometry_change: the processor grid has changed. */
#define CELL_FLAG_GRIDCHANGED 1
/** Flag for cells_on_geometry_change: skip shrinking of cells. */
#define CELL_FLAG_FAST 2

/*@}*/


/************************************************/
/** \name Data Types */
/************************************************/
/*@{*/

/** A cell is a \ref ParticleList representing a particle group with
    respect to the integration algorithm.
*/
typedef ParticleList Cell;

/** List of cell pointers. */
typedef struct {
  Cell **cell;
  int n;
  int max;
} CellPList;

/** Describes a cell structure / cell system. Contains information
    about the communication of cell contents (particles, ghosts, ...) 
    between different nodes and the relation between particle
    positions and the cell system. All other properties of the cell
    system which are not common between different cell systems have to
    be stored in seperate structures. */
typedef struct {
  /** type descriptor */
  int type;

  /** Communicator to exchange ghost cell information. */
  GhostCommunicator ghost_cells_comm;
  /** Communicator to exchange ghost particles. */
  GhostCommunicator exchange_ghosts_comm;
  /** Communicator to update ghost positions. */
  GhostCommunicator update_ghost_pos_comm;
  /** Communicator to collect ghost forces. */
  GhostCommunicator collect_ghost_force_comm;
#ifdef LB
  /** Communicator for particle data used by lattice Boltzmann */
  GhostCommunicator ghost_lbcoupling_comm;
#endif

  /** Cell system dependent function to find the right node for a
      particle at position pos. 
      \param  pos Position of a particle.
      \return number of the node where to put the particle. 
  */
  int   (*position_to_node)(double pos[3]);
  /** Cell system dependent function to find the right cell for a
      particle at position pos. 
      \param  pos Position of a particle.
      \return pointer to cell  where to put the particle.
  */
  Cell *(*position_to_cell)(double pos[3]);
} CellStructure;

/*@}*/

/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** list of all cells. */
extern Cell *cells;
/** size of \ref cells::cells */
extern int n_cells;
/** list of all cells containing particles physically on the local
    node */
extern CellPList local_cells;
/** list of all cells containing ghosts */
extern CellPList ghost_cells;

/** Type of cell structure in use ( \ref Cell Structure ). */
extern CellStructure cell_structure;

/** Maximal interaction range - also the minimum cell size. Any
    cellsystem makes sure that the particle pair loop visits all pairs
    of particles that are closer than this. */
extern double max_range;

/** If non-zero, cell systems should reset the position for checking
    the Verlet criterion. Moreover, the Verlet list has to be
    rebuilt. */
extern int rebuild_verletlist;

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Reinitialize the cell structures.
    @param new_cs gives the new topology to use afterwards. May be set to
    \ref CELL_STRUCTURE_CURRENT for not changing it.
*/
void cells_re_init(int new_cs);

/** Reallocate the list of all cells (\ref cells::cells). */
void realloc_cells(int size);

/** Initialize a list of cell pointers */
inline void init_cellplist(CellPList *cpl) {
  cpl->n    = 0;
  cpl->max  = 0;
  cpl->cell = NULL;
}

/** Reallocate a list of cell pointers */
inline void realloc_cellplist(CellPList *cpl, int size)
{
  if(size != cpl->max) {
    cpl->max = size;
    cpl->cell = (Cell **) realloc(cpl->cell, sizeof(Cell *)*cpl->max);
  }
}

/** sort the particles into the cells and initialize the ghost particle structures.
    @param global_flag if this is CELLS_GLOBAL_EXCHANGE, particle positions can have changed
    arbitrarly, otherwise the change should have been smaller then skin.  */
void cells_resort_particles(int global_flag);

/** this function is called whenever the cell system has to be
    reinitialized, e.g. if cutoffs have changed, or the skin, grid,
    .... It calculates the maximal interaction range, and
    as said reinitializes the cells structure if something significant
    has changed.

    If bit CELL_FLAG_FAST is set, the routine should try to save time.
    Currently this means that if the maximal range decreased, it does
    not reorganize the particles. This is used in the NpT algorithm to
    avoid frequent reorganizaion of particles.

    If bit CELL_FLAG_GRIDCHANGED is set, it means the nodes' topology
    has changed, i. e. the grid or periodicity. In this case a full
    reorganization is due.

    @param flags a combination of CELL_FLAG_GRIDCHANGED and
    CELL_FLAG_FAST, see above.

*/
void cells_on_geometry_change(int flags);

/** update ghost information. If \ref resort_particles == 1,
    also a resorting of the particles takes place. */
void cells_update_ghosts();

/** Calculate and return the total number of particles on this
    node. */
int cells_get_n_particles();

/** Debug function to print particle positions. */
void print_local_particle_positions();

/** Debug function to print ghost positions. */
void print_ghost_positions();

/** spread the particle resorting criterion across the nodes. */
void announce_resort_particles();

/* Checks if a particle resorting is required. */
void check_resort_particles();

/* Do a strict particle sorting, including order in the cells. */
void local_sort_particles();

/*@}*/

#endif

#ifndef DOMAIN_DECOMP_H
#define DOMAIN_DECOMP_H
#include "cells.h"

/** Structure containing the information about the cell grid used for domain decomposition. */
typedef struct {
  /** linked cell grid in nodes spatial domain. */
  int cell_grid[3];
  /** linked cell grid with ghost frame. */
  int ghost_cell_grid[3];
  /** cell size. 
      Def: \verbatim cell_grid[i] = (int)(local_box_l[i]/max_range); \endverbatim */
  double cell_size[3];
  /** inverse cell size = \see cell_size ^ -1. */
  double inv_cell_size[3];
}  DomainDecomposition;


/** Structure containing information about non bonded interactions
    with particles in a neighbor cell. */
typedef struct {
  /** Just for transparency the index of the neighbor cell. */
  int cell_ind;
  /** Pointer to particle list of neighbor cell. */
  ParticleList *pList;
  /** Verlet list for non bonded interactions of a cell with a neighbor cell. */
  PairList vList;
} IA_Neighbor;


typedef struct {
  /** number of interacting neighbor cells . 

      A word about the interacting neighbor cells:

      In a 3D lattice each cell has 27 neighbors (including
      itself!). Since we deal with pair forces, it is sufficient to
      calculate only half of the interactions (Newtons law: actio =
      reactio). For each cell 13+1=14 neighbors. This has only to be
      done for the inner cells. 

      Caution: This implementation needs double sided ghost
      communication! For single sided ghost communication one would
      need some ghost-ghost cell interaction as well, which we do not
      need! 

      It follows: inner cells: n_neighbors = 14
      ghost cells:             n_neighbors = 0
  */
  int n_neighbors;
  /** Interacting neighbor cell list  */
  IA_Neighbor *nList;
} IA_Neighbor_List;


/************************************************************/
/** \name Exported Variables */
/************************************************************/
/*@{*/

/** Information about the domain decomposition. */
extern DomainDecomposition dd;

/** maximal skin size. */
extern double max_skin;
/** Maximal number of cells per node. In order to avoid memory
 *  problems due to the cell grid one has to specify the maximal
 *  number of \ref #cells . The corresponding callback function is
 *  \ref max_num_cells_callback. If the number of cells \ref n_cells,
 *  defined by \ref ghost_cell_grid is larger than max_num_cells the
 *  cell grid is reduced. max_num_cells has to be larger than 27, e.g
 *  one inner cell.  max_num_cells is initialized with the default
 *  value specified in \ref config.h: \ref CELLS_MAX_NUM_CELLS.
 */
extern int max_num_cells;

/*@}*/

/************************************************************/
/** \name Exported Functions */
/************************************************************/
/*@{*/

/** Initialize the topology. The argument is list of cells, which particles have to be
    sorted into their cells. The particles might not belong to this node.
    This procedure is used when particle data or cell structure has changed and
    the cell structure has to be reinitialized. This also includes setting up the
    cell_structure array. */
void dd_topology_init(CellPList *cl);

/** Called when the current cell structure is invalidated because for example the
    box length has changed. This procedure may NOT destroy the old inner and ghost
    cells, but it should free all other organizational data. Note that parameters
    like the box length or the node_grid may already have changed. Therefore
    organizational data has to be stored independently from variables
    that may be changed from outside. */
void dd_topology_release();

/** Just resort the particles. Used during integration. The particles are stored in
    the cell structure. Domain decomposition can assume for example that particles
    only have to be sent to neighboring nodes. */
void  dd_exchange_and_sort_particles();

/** implements \ref cell_structure::position_to_cell. */
Cell *dd_position_to_cell(double pos[3]);

/** implements \ref cell_structure::position_to_cell. */
Cell *dd_position_to_node(double pos[3]);

/** Callback for setmd maxnumcells (maxnumcells >= 27). 
    see also \ref max_num_cells */
int max_num_cells_callback(Tcl_Interp *interp, void *_data);


/*@}*/

#endif

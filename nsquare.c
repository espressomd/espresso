#include "nsquare.h"
#include "communication.h"
#include "debug.h"

Cell *local;
CellPList me_ghosts;
CellPList other_ghosts;

Cell *nsq_position_to_cell(double pos[3])
{
  return local;
}

void nsq_topology_release()
{
  realloc_cellplist(&me_ghosts, 0);
  realloc_cellplist(&other_ghosts, 0);
}


void nsq_topology_init(CellPList *old)
{
  Particle *part;
  int c, p, np;

  CELL_TRACE(fprintf(stderr, "%d: nsq_topology_init, %d\n", this_node, old->n));
  cell_structure.type = CELL_STRUCTURE_NSQUARE;
  cell_structure.position_to_node = map_position_node_array;
  cell_structure.position_to_cell = nsq_position_to_cell;

  realloc_cells(n_nodes);

  local = &cells[0];

  /* mark cells */
  realloc_cellplist(&local_cells, local_cells.n = 1);
  local_cells.cell[0] = local;
  realloc_cellplist(&ghost_cells, ghost_cells.n = n_nodes - 1);
  for (c = 1; c < n_nodes; c++)
    ghost_cells.cell[c - 1] = &cells[c];

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np   = old->cell[c]->n;
    for (p = 0; p < np; p++)
      append_unindexed_particle(local, &part[p]);
  }
  update_local_particles(local);

  // ghosts not yet done
  init_cellplist(&me_ghosts);
  init_cellplist(&other_ghosts);
}


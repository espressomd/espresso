#include <mpi.h>
#include <string.h>
#include "nsquare.h"
#include "communication.h"
#include "debug.h"
#include "ghosts.h"

Cell *local;
CellPList me_do_ghosts;

Cell *nsq_position_to_cell(double pos[3])
{
  return local;
}

void nsq_topology_release()
{
  CELL_TRACE(fprintf(stderr,"%d: nsq_topology_release:\n",this_node));
  /* free ghost cell pointer list */
  realloc_cellplist(&me_do_ghosts, 0);
  free_comm(&cell_structure.ghost_cells_comm);
  free_comm(&cell_structure.exchange_ghosts_comm);
  free_comm(&cell_structure.update_ghost_pos_comm);
  free_comm(&cell_structure.collect_ghost_force_comm);
}

static void nsq_prepare_comm(GhostCommunicator *comm, int data_parts)
{
  int n;
  prepare_comm(comm, data_parts, n_nodes);
  /* every node has its dedicated comm step */
  for(n = 0; n < n_nodes; n++) {
    comm->comm[n].part_lists = malloc(sizeof(ParticleList *));
    comm->comm[n].part_lists[0] = &cells[n];
    comm->comm[n].n_part_lists = 1;
    comm->comm[n].mpi_comm = MPI_COMM_WORLD;
  }
}

void nsq_topology_init(CellPList *old)
{
  Particle *part;
  int n, c, p, np, ntodo, diff;

  CELL_TRACE(fprintf(stderr, "%d: nsq_topology_init, %d\n", this_node, old->n));

  if (cell_structure.type != CELL_STRUCTURE_NSQUARE) {
    cell_structure.type = CELL_STRUCTURE_NSQUARE;
    cell_structure.position_to_node = map_position_node_array;
    cell_structure.position_to_cell = nsq_position_to_cell;

    realloc_cells(n_nodes);

    /* mark cells */
    local = &cells[this_node];
    realloc_cellplist(&local_cells, local_cells.n = 1);
    local_cells.cell[0] = local;

    realloc_cellplist(&ghost_cells, ghost_cells.n = n_nodes - 1);
    c = 0;
    for (n = 0; n < n_nodes; n++)
      if (n != this_node)
	ghost_cells.cell[c++] = &cells[n];

    /* distribute force calculation work  */
    ntodo = (n_nodes + 3)/2;
    init_cellplist(&me_do_ghosts);
    realloc_cellplist(&me_do_ghosts, ntodo);
    for (n = 0; n < n_nodes; n++) {
      diff = n - this_node;
      /* simple load balancing formula. Basically diff % n, where n >= n_nodes, n odd. */
      if (((diff >= 0 && (diff % 2) == 0) ||
	   (diff < 0 && ((-diff) % 2) == 1))) {
	CELL_TRACE(fprintf(stderr, "%d: doing interactions with %d\n", this_node, n));
	me_do_ghosts.cell[me_do_ghosts.n++] = &cells[n];
      }
    }

    /* create communicators */
    nsq_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);
    nsq_prepare_comm(&cell_structure.exchange_ghosts_comm,     GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
    nsq_prepare_comm(&cell_structure.update_ghost_pos_comm,    GHOSTTRANS_POSITION);
    nsq_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);

    /* here we just decide what to transfer where */
    /* first round: all nodes prefetch their send data */
    cell_structure.ghost_cells_comm.comm[n].type         = GHOST_BCST | GHOST_PREFETCH;
    cell_structure.exchange_ghosts_comm.comm[n].type     = GHOST_BCST | GHOST_PREFETCH;
    cell_structure.update_ghost_pos_comm.comm[n].type    = GHOST_BCST | GHOST_PREFETCH;
    /* here all nodes send every comm step, no prefetch */
    cell_structure.collect_ghost_force_comm.comm[n].type = GHOST_RDCE;

    for (n = 1; n < n_nodes; n++) {
      /* use the prefetched send buffers */
      if (this_node != n) {
	cell_structure.ghost_cells_comm.comm[n].type         = GHOST_BCST;
	cell_structure.exchange_ghosts_comm.comm[n].type     = GHOST_BCST;
	cell_structure.update_ghost_pos_comm.comm[n].type    = GHOST_BCST;
      }
      else {
	cell_structure.ghost_cells_comm.comm[n].type         = GHOST_BCST | GHOST_PREFETCH;
	cell_structure.exchange_ghosts_comm.comm[n].type     = GHOST_BCST | GHOST_PREFETCH;
	cell_structure.update_ghost_pos_comm.comm[n].type    = GHOST_BCST | GHOST_PREFETCH;
      }
      cell_structure.collect_ghost_force_comm.comm[n].type = GHOST_RDCE;
    }
  }

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np   = old->cell[c]->n;
    for (p = 0; p < np; p++)
      append_unindexed_particle(local, &part[p]);
  }
  update_local_particles(local);
}

void nsq_balance_particles()
{
  int i, n, surplus, s_node, tmp, lack, l_node, transfer;

  int pp = cells_get_n_particles();
  int *ppnode = malloc(n_nodes*sizeof(int));
  /* minimal difference between node shares */
  int minshare = n_total_particles/n_nodes;
  int maxshare = minshare + 1;

  
  MPI_Allgather(&pp, 1, MPI_INT, ppnode, 1, MPI_INT, MPI_COMM_WORLD);
  for (;;) {
    /* find node with most excessive particles */
    surplus = 0;
    s_node = -1;
    for (n = 0; n < n_nodes; n++) {
      tmp = ppnode[n] - maxshare;
      if (tmp > surplus) {
	surplus = tmp;
	s_node = n;
      }
    }

    /* find node with least excessive particles */
    lack = 0;
    l_node = -1;
    for (n = 0; n < n_nodes; n++) {
      tmp = minshare - ppnode[n];
      if (tmp < lack) {
	lack = tmp;
	l_node = n;
      }
    }

    /* all nodes ok */
    if (s_node == -1 && l_node == -1)
      break;

    /* should not happen: minshare or maxshare wrong or more likely,
       the algorithm */
    if (s_node == -1 || l_node == -1) {
      fprintf(stderr, "%d: Particle load balancing failed\n", this_node);
      break;
    }

    transfer = lack < surplus ? lack : surplus;

    if (s_node == this_node) {
      ParticleList send_buf;
      init_particleList(&send_buf);
      realloc_particles(&send_buf, send_buf.n = transfer);
      for (i = 0; i < transfer; i++) {
	memcpy(&send_buf.part[i], &local->part[--local->n], sizeof(Particle));
      }
      realloc_particles(local, local->n);
      send_particles(&send_buf, l_node);
    }
    else if (l_node == this_node) {
      recv_particles(local, n);
    }

    ppnode[s_node] -= transfer;
    ppnode[l_node] += transfer;
  }
  free(ppnode);
}

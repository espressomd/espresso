#include <mpi.h>
#include <string.h>
#include "layered.h"
#include "communication.h"
#include "debug.h"
#include "ghosts.h"
#include "forces.h"
#include "pressure.h"
#include "energy.h"

/* Organization: Layers only in one direction.
   ghost_bottom
   c1
   c2
   c3
   .
   .
   .
   cn
   ghost_top

   First, all nodes send downwards, then upwards. Within these actions,
   first the odd nodes send. For even n_nodes, this algorithm is straight
   forward: first the odd nodes send, the even receive, then vice versa.
   For odd n_nodes, we have
   1) 1->2 3->4 5->1
   2) 2->3 4->5
   So in the first round node 5 has to wait for node 1 to
   complete the send and get ready to receive. In other words,
   what physically happens is:
   1) 1->2 3->4 5->*
   2) *->1 2->3 4->*
   3) *->5
   This means that one pending receive has to be done in addition
   provided that all other transmissions can happen in parallel.

*/

/** wether we are the lowest node */
#define LAYERED_BOTTOM 1
/** wether we are the highest node */
#define LAYERED_TOP    2
/** same as PERIODIC(2) */
#define LAYERED_PERIODIC 4
#define LAYERED_BTM_MASK (LAYERED_BOTTOM|LAYERED_PERIODIC)
#define LAYERED_TOP_MASK (LAYERED_TOP|LAYERED_PERIODIC)
/** node has a neighbor above (modulo n_nodes) */
#define LAYERED_TOP_NEIGHBOR ((layered_flags & LAYERED_TOP_MASK) != LAYERED_TOP)
/** node has a neighbor below (modulo n_nodes) */
#define LAYERED_BTM_NEIGHBOR ((layered_flags & LAYERED_BTM_MASK) != LAYERED_BOTTOM)

int layered_flags = 0;
int n_layers = 0, btm, top;
double layer_h = 0, layer_h_i = 0;

Cell *layered_position_to_cell(double pos[3])
{
  int cpos = (int)((pos[2] - my_left[2])*layer_h_i) + 1;
  if (cpos < 1) {
    if (!LAYERED_BTM_NEIGHBOR)
      cpos = 1;
    else
      return NULL;
  }
  else if (cpos > n_layers) {
    /* not periodic, but at top */
    if (!LAYERED_TOP_NEIGHBOR)
      cpos = n_layers;
    else
      return NULL;
  }
  return &(cells[cpos]);  
}

void layered_topology_release()
{
  CELL_TRACE(fprintf(stderr,"%d: layered_topology_release:\n", this_node));
  free_comm(&cell_structure.ghost_cells_comm);
  free_comm(&cell_structure.exchange_ghosts_comm);
  free_comm(&cell_structure.update_ghost_pos_comm);
  free_comm(&cell_structure.collect_ghost_force_comm);
}

static void layered_prepare_comm(GhostCommunicator *comm, int data_parts)
{
  GhostCommunication *gc;
  int even_odd;
  int c, n = 4; /* how many communications to do by default: up even/odd, down even/odd */
  
  /* if periodic and bottom or top, send shifted */
  if (((layered_flags & LAYERED_TOP_MASK) == LAYERED_TOP_MASK) ||
      ((layered_flags & LAYERED_BTM_MASK) == LAYERED_BTM_MASK))
    data_parts |= GHOSTTRANS_POSSHFTD;

  /* one communication missing if not periodic but on border */
  if (!LAYERED_TOP_NEIGHBOR)
    n -= 2;
  if (!LAYERED_BTM_NEIGHBOR)
    n -= 2;

  prepare_comm(comm, data_parts, n);
  /* always sending/receiving 1 cell per time step */
  for(c = 0; c < n; c++) {
    comm->comm[c].part_lists = malloc(sizeof(ParticleList *));
    comm->comm[c].n_part_lists = 1;
    comm->comm[c].mpi_comm = MPI_COMM_WORLD;
  }

  gc = comm->comm;

  /* downwards */
  for (even_odd = 0; even_odd < 2; even_odd++) {
    /* send */
    if (this_node % 2 == even_odd && LAYERED_TOP_NEIGHBOR) {
      gc->type = GHOST_SEND;
      gc->node = btm;
      if (data_parts == GHOSTTRANS_FORCE)
	gc->part_lists[0] = &cells[0];
      else {
	gc->part_lists[0] = &cells[1];
	if (data_parts & GHOSTTRANS_POSSHFTD) {
	  gc->shift[0] = gc->shift[1] = 0;
	  gc->shift[2] = box_l[2];
	}
      }
      gc++;
    }
    /* recv. Note we test r_node as we always have to test the sender
       as for odd n_nodes maybe we send AND receive. */
    if (top % 2 == even_odd && LAYERED_BTM_NEIGHBOR) {
      gc->type = GHOST_RECV;
      gc->node = top;
      if (data_parts == GHOSTTRANS_FORCE)
	gc->part_lists[0] = &cells[n_layers];
      else
	gc->part_lists[0] = &cells[n_layers + 1];
      gc++;
    }
  }

  /* upwards */
  for (even_odd = 0; even_odd < 2; even_odd++) {
    /* send */
    if (this_node % 2 == even_odd && LAYERED_TOP_NEIGHBOR) {
      gc->type = GHOST_SEND;
      gc->node = top;
      if (data_parts == GHOSTTRANS_FORCE)
	gc->part_lists[0] = &cells[n_layers + 1];
      else {
	gc->part_lists[0] = &cells[n_layers];
	if (data_parts & GHOSTTRANS_POSSHFTD) {
	  gc->shift[0] = gc->shift[1] = 0;
	  gc->shift[2] = -box_l[2];
	}
      }
      gc++;
    }
    /* recv. Note we test r_node as we always have to test the sender
       as for odd n_nodes maybe we send AND receive. */
    if (btm % 2 == even_odd && LAYERED_BTM_NEIGHBOR) {
      gc->type = GHOST_RECV;
      gc->node = btm;
      if (data_parts == GHOSTTRANS_FORCE)
	gc->part_lists[0] = &cells[1];
      else
	gc->part_lists[0] = &cells[0];
      gc++;
    }
  }
}

void layered_topology_init(CellPList *old)
{
  Particle *part;
  int c, p, np;

  CELL_TRACE(fprintf(stderr, "%d: layered_topology_init, %d old particles\n", this_node, old->n));

  if (cell_structure.type != CELL_STRUCTURE_LAYERED) {
    cell_structure.type = CELL_STRUCTURE_LAYERED;
    cell_structure.position_to_node = map_position_node_array;
    cell_structure.position_to_cell = layered_position_to_cell;
  }

  top = (this_node + 1) % n_nodes;
  btm = (this_node - 1) % n_nodes;

  layer_h = box_l[2]/(double)(n_layers);
  layer_h_i = 1/layer_h;

  if (this_node == 0)
    layered_flags = LAYERED_BOTTOM;
  else if (this_node == n_nodes - 1)
    layered_flags = LAYERED_TOP;

  if (PERIODIC(2))
    layered_flags |= LAYERED_PERIODIC;

  /* allocate cells and mark them */
  realloc_cells(n_layers + 2);
  realloc_cellplist(&local_cells, local_cells.n = n_layers);
  for (c = 0; c < n_layers; c++)
    local_cells.cell[c] = &cells[c + 1];
  realloc_cellplist(&ghost_cells, ghost_cells.n = 2);
  ghost_cells.cell[0] = &cells[0];
  ghost_cells.cell[1] = &cells[n_layers + 1];

  /* create communicators */
  layered_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);
  layered_prepare_comm(&cell_structure.exchange_ghosts_comm,     GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
  layered_prepare_comm(&cell_structure.update_ghost_pos_comm,    GHOSTTRANS_POSITION);
  layered_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np   = old->cell[c]->n;
    for (p = 0; p < np; p++)
      append_unindexed_particle(layered_position_to_cell(part[p].r.p), &part[p]);
  }
  for (c = 1; c <= n_layers; c++)
    update_local_particles(&cells[c]);
  CELL_TRACE(fprintf(stderr, "%d: layered_topology_init done\n", this_node));
}
 
static void layered_append_particles(ParticleList *pl)
{
  int p;

  for(p = 0; p < pl->n; p++) {
    if (LAYERED_BTM_NEIGHBOR && pl->part[p].r.p[2] < my_left[2])
      continue;
    if (LAYERED_TOP_NEIGHBOR && pl->part[p].r.p[2] > my_right[2])
      continue;

    move_indexed_particle(layered_position_to_cell(pl->part[p].r.p), pl, p);
  }
}

void layered_exchange_and_sort_particles(int global_flag)
{
  Particle *part;
  Cell *nc, *oc;
  int c, p, flag, redo;
  ParticleList send_buf_dn, send_buf_up;
  ParticleList recv_buf;

  init_particleList(&send_buf_dn);
  init_particleList(&send_buf_up);

  init_particleList(&recv_buf);

  for (;;) {
    /* sort local particles */
    for (c = 1; c <= n_layers; c++) {
      oc = &cells[c];

      for (p = 0; p < oc->n; p++) {
	part = &oc->part[p];
	nc = layered_position_to_cell(part->r.p);

	if (LAYERED_BTM_NEIGHBOR && part->r.p[2] < my_left[2])
	  move_indexed_particle(&send_buf_dn, oc, p);
	else if (LAYERED_TOP_NEIGHBOR && part->r.p[2] > my_right[2])
	  move_indexed_particle(&send_buf_up, oc, p);
	else {
	  if (nc != oc) {
	    move_indexed_particle(nc, oc, p);
	    if (p < oc->n) p--;
	  }
	}
      }
    }

    /* transfer */
    if (this_node % 2 == 0) {
      /* send down */
      if (LAYERED_BTM_NEIGHBOR)
	send_particles(&send_buf_dn, btm);
      if (LAYERED_TOP_NEIGHBOR)
	recv_particles(&recv_buf, top);
      /* send up */
      if (LAYERED_BTM_NEIGHBOR)
	send_particles(&send_buf_up, top);
      if (LAYERED_TOP_NEIGHBOR)
	recv_particles(&recv_buf, btm);
    }
    else {
      if (LAYERED_TOP_NEIGHBOR)
	recv_particles(&recv_buf, top);
      if (LAYERED_BTM_NEIGHBOR)
	send_particles(&send_buf_dn, btm);
      if (LAYERED_BTM_NEIGHBOR)
	recv_particles(&recv_buf, btm);
      if (LAYERED_TOP_NEIGHBOR)
	send_particles(&send_buf_up, top);
    }
    
    layered_append_particles(&recv_buf);

    /* handshake redo */
    flag = (recv_buf.n != 0);

    if (global_flag == LAYERED_FULL_EXCHANGE) {
      MPI_Allreduce(&flag, &redo, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
      if (!redo)
	break;
    }
    else {
      if (flag) {
	fprintf(stderr,"%d: layered_exchange_and_sort_particles: particle moved to far\n", this_node);
	errexit();
      }
    }
  }
}

/** nonbonded and bonded force calculation using the verlet list */
void layered_calculate_ia()
{
  int c, i, j;
  Cell *celll, *cellt, *cellb;
  int      npl, npt, npb;
  Particle *pl, *pt, *pb, *p1;
  double dist2, d[3];
 
  for (c = 1; c <= n_layers; c++) {
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;

    if (LAYERED_BTM_NEIGHBOR) {
      cellb = &cells[c-1];
      pb    = cellb->part;
      npb   = cellb->n;
    }
    else {
      npb   = 0;
      pb    = NULL;
    }

    if (LAYERED_TOP_NEIGHBOR) {
      cellt = &cells[c+1];
      pt    = cellt->part;
      npt   = cellt->n;
    }
    else {
      npt   = 0;
      pt    = NULL;
    }

    for(i = 0; i < npl; i++) {
      p1 = &pl[i];

      add_bonded_force(p1);
#ifdef CONSTRAINTS
      add_constraints_forces(p1);
#endif

      /* cell itself and bonded / constraints */
      for(j = i+1; j < npl; j++) {
	dist2 = distance2vec(p1->r.p, pl[j].r.p, d);
	add_non_bonded_pair_force(p1, &pl[j], d, sqrt(dist2), dist2);
      }

      /* top neighbor */
      for(j = 0; j < npt; j++) {
	dist2 = distance2vec(p1->r.p, pt[j].r.p, d);
	add_non_bonded_pair_force(p1, &pt[j], d, sqrt(dist2), dist2);
      }

      /* bottom neighbor */
      for(j = 0; j < npb; j++) {
	dist2 = distance2vec(p1->r.p, pb[j].r.p, d);
	add_non_bonded_pair_force(p1, &pb[j], d, sqrt(dist2), dist2);
      }
    }
  }
}

void layered_calculate_energies()
{
}

void layered_calculate_virials()
{
}

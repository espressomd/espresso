/** \file ghosts.c
 *
 *  Ghost particles and particle exchange.
 *
 *  In this file you find everything concerning the exchange of
 *  particle data (particles, ghosts, positions and forces) for short
 *  range interactions between the spacial domains of neighbouring
 *  nodes.
 *
 *  All structures are initialized by \ref ghost_init.  Exchanging
 *  particles that move from one to another node domain is done by
 *  \ref exchange_part. Setup of the ghost particle frame is done by
 *  \ref exchange_ghost. Call \ref sort_particles_into_cells before
 *  each use of \ref exchange_ghost. During integration using a verlet
 *  list ghost positions and forces are exchanged by \ref
 *  update_ghost_pos and \ref collect_ghost_forces. Before calling
 *  \ref ghost_init again tou should clean up with \ref ghost_exit.
 *
 *  Communication \anchor communication is done only between
 *  neighboring nodes. The communication scheme can be senn in the
 *  figure below.
 *
 *  \image html ghost_communication.gif "Scheme of ghost/particle communication"
 *  \image latex ghost_communication.eps "Scheme of ghost/particle communication" width=8cm 
 *
 *  To reduce the number of communications per node from 26 (number of
 *  neighbor nodes in 3 dimensions) to 6 the order of the
 *  communication is important:
 *  <ol> 
 *      <li> x-direction: left and right - MPI_Barrier
 *      <li> y-direction: forth and back - MPI_Barrier
 *      <li> z-direction: down and up - MPI_Barrier
 *  </ol>
 *  In this way also edges and corners are communicated
 *  (See also \ref directions for our conventions).
 *
 *  \warning \b Since the ghost particle structures make use of the
 *  linked cell structure, \ref ghost_init has to be called after \ref
 *  cells_init */
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "ghosts.h"
#include "debug.h"
#include "global.h"
#include "cells.h"
#include "communication.h"
#include "grid.h"
#include "utils.h"

/************************************************
 * defines
 ************************************************/

/**  granularity of the communication buffers. */
#define PART_INCREMENT 20

/************************************************
 * variables
 ************************************************/

/** \name Variables for particle exchange */
/*@{*/
/** Buffer for particles to send. */
Particle *p_send_buf;
int       n_p_send_buf;
int       max_p_send_buf;
/** Buffer for particles to recieve. */
Particle *p_recv_buf;
int       n_p_recv_buf;
int       max_p_recv_buf;
/** Buffer for particles bonds to send. */
int      *b_send_buf;
int       n_b_send_buf;
int       max_b_send_buf;
/** Buffer for particles bonds recieve. */
int      *b_recv_buf;
int       n_b_recv_buf;
int       max_b_recv_buf;
/*@}*/

/** \name Variables for ghost particle exchange */
/*@{*/
/** Buffer for Ghosts to send. */
Ghost *g_send_buf;
int   n_g_send_buf;
int   max_g_send_buf;
/** Buffer for Ghosts to recieve. */
Ghost *g_recv_buf;
int   n_g_recv_buf;
int   max_g_recv_buf;

/** list of cell indices to send. */
int *send_cells;
/** list of cell indices to receive. */
int *recv_cells;
/** total number of send/recv cells */
int ntot_send_cells;

/** maximal number of cells to send in one direction. */
int max_send_cells;
/** number of cells to send in direction X. */
int n_send_cells[6];
/** number of cells to receive from direction X. */
int n_recv_cells[6];
/** start indices for cells to send/recv in/from direction X. */ 
int cell_start[6];

/** Number of ghosts in each send cell. */ 
int *n_send_ghosts;
/** Number of ghosts in each recv cell. */ 
int *n_recv_ghosts;

/** number of ghosts to send in direction X */
int ghost_send_size[6];
/** number of ghosts to recv from direction X */
int ghost_recv_size[6];
/*@}*/

/** \name Variables for ghost force/position exchange */
/*@{*/
/** Buffer for forces/coordinates to send. */
double *send_buf;
int max_send_buf;
/** Buffer for forces/coordinates to recieve. */
double *recv_buf;
int max_recv_buf;
/*@}*/

/************************************************
 * privat functions
 ************************************************/

/** Creates a linear index list of a sub grid.
 *
 *  The sub grid is defined by its lower and upper corner:
 *  from (lc[0],lc[1],lc[2]) to (hc[0],hc[1],hc[2])
 *  The grid dimension is given with (gs[0],gs[1],gs[2])
 *  The linear index list of length <returnvalue> is stored 
 *  in list starting from position start. max should be the 
 *  total length of list to ensure that the indices fit into list.
 *
 * @return         size of the subgrid.
 * @param *list    array to store the indices of the sub grid.
 * @param start    start index for sub grid indices.
 * @param max      size of the array list.
 * @param lc[3]    lower corner of sub grid.
 * @param hc[3]    upper corner of sub grid.
 * @param gs[3]    grid dimension .
 */  
int sub_grid_indices(int* list, int start, int max,
		     int lc[3], int hc[3], int gs[3]);

/** moves particle (ind) to the send buffers. 
 * 
 *  subroutine of \ref exchange_part .  
 *
 *  Moves one particle (struct: Particle, local index: ind) to the
 *  send buffer (p_send_buf, b_send_buf) and removes it from the local
 *  particle field. The empty space in the local particle field is
 *  then filled with the last particle (See illustration).
 *
 *  \image html move_to_p_buf.gif "Move particle to send buffer" 
 *  \image latex move_to_p_buf.eps "Move particle to send buffer" \width=10cm
 *  
 *  Variable changes:
 *  <ul> 
 *    <li> n_p_send_buf += 1
 *    <li> n_b_send_buf += number of bond partners
 *    <li> remove particle from local_index field (set to -1)
 *    <li> actualize local_index entry for last particle
 *    <li> n_particles -= 1
 *  </ul>
 *  Reallocation of p_send_buf and b_send_buf if necessary.  
 *
 * \warning \b Supports only two particle bonds at the moment
 *
 * @return    local index to continue the particle loop.     
 * @param ind local index of the particle. */
int move_to_p_buf(int ind);

/** send particles in direction s_dir.
 *
 *  subroutine of \refexchange_part. 
 *
 *  Check if communication goes to a different node.
 *  <ol>
 *    <li> Communication to a different node:
           Two step communication (first all even PEs send than all odd)
	   Send buffer sizes.
	   Reallocate recieve buffers if necessary.
	   Send particle and bond buffers.
      <li> Communication inside one node
           Exchange adress and sizes of send and recieve buffers.
 *  </ol>
 *
 * \warning \b Supports only two particle bonds at the moment
 *
 * @param s_dir    send direction. 
*/
void send_particles(int s_dir);

/** appends recieved particles to local particle array.
 *
 *  subroutine of \ref exchange_part. 
 *
 *  Reallocate particle buffer (particles) if necessary. Copy
 *  particles and their bonds from recieve buffers to local particle
 *  buffer.
 *
 * \warning \b Supports only two particle bonds at the moment.
 *
*/
void append_particles(void);

/** Send ghost particles in direction s_dir.
 *
 *  subroutine of \ref exchange_ghosts. 
 *
 *  Does an unbuffered communication from all nodes to their neighbor
 *  node in direction s_dir.
 *  <ol>
 *      <li> send number of ghost in each send cell for this direction:
 *           \ref n_send_ghosts to \ref n_recv_ghosts.
 *           The total number of ghosts to send/receive is in 
 *           \ref n_send_ghosts[\ref max_send_cells].
 *      <li> send ghost particle information:
 *           \ref g_send_buf to \ref g_recv_buf
 *  </ol>
 *  \ref g_recv_buf is reallocated if necessary.
 *
 *  If communication goes to the same node, just the pointers of the
 *  array pairs (\ref g_send_buf, \ref g_recv_buf) and the variable
 *  pairs (\ref max_g_send_buf, \ref max_g_recv_buf) and (\ref
 *  n_g_send_buf, \ref n_g_recv_buf) are exchanged.
 *
 * @param s_dir send direction.  
*/
void send_ghosts(int s_dir);

/** Copy identity, type, position and charge from particle to ghost array.
 * 
 *  Attention: no array checks are performed.
 *
 * @param g_array   ghost particle array.
 * @param g_ind     index for ghost particle.
 * @param p_array   particle array.
 * @param p_ind     index of particle.
 */
void pack_ghost(Ghost *g_array, int g_ind, Particle *p_array, int p_ind);

/** Copy identity, type, position and charge from ghost to particle array.
 * 
 *  Attention: no array checks are performed.
 *
 * @param p_array   particle array.
 * @param p_ind     index for particle.
 * @param g_array   ghost particle array.
 * @param g_ind     index of ghost particle.
 */
void unpack_ghost(Particle *p_array, int p_ind, Ghost *g_array, int g_ind);

/** send positions/forces in direction s_dir.
 *
 *  Does an unbuffered communication from all nodes to their neighbor
 *  node in direction s_dir.
 *
 *  send positions / forces: \ref send_buf to \ref recv_buf.
 *  If communication goes to the same node just the pointers 
 *  of \ref send_buf and \ref recv_buf are exchanged.
 *
 * @param s_dir       direction to send to/recv from.
 * @param send_size   number of positions/forces to send. 
 * @param recv_size   number of positions/forces to receive.
 */
void send_posforce(int s_dir,int send_size, int recv_size);

/************************************************
 * public functions
 ************************************************/

void ghost_init()
{
  int i;
  /* ghost cell grid, cell grid */
  int gcg[3],cg[3];
  /* cell list sizes for the 3 direction, total number of sed cellsm. */
  int anz[3];
  /* send/recv frames start indizes, end indizes */
  int lc[3],hc[3],done[3]={0,0,0};

  GHOST_TRACE(fprintf(stderr,"%d: ghost_init:\n",this_node));

  /* Init PE node_neighbors */
  calc_node_neighbors(this_node);

  /* Init exchange particles */
  max_p_send_buf = max_p_recv_buf = PART_INCREMENT;
  p_send_buf = (Particle *)malloc(max_p_send_buf*sizeof(Particle));
  p_recv_buf = (Particle *)malloc(max_p_recv_buf*sizeof(Particle));
  max_b_send_buf = max_b_recv_buf = PART_INCREMENT;
  b_send_buf = (int *)malloc(max_b_send_buf*sizeof(int));
  b_recv_buf = (int *)malloc(max_b_recv_buf*sizeof(int));


  /* Init ghost exchange */
  /* preparation of help variables */
  for(i=0;i<3;i++) {
    gcg[i] = ghost_cell_grid[i];
    cg[i]  = cell_grid[i];
  }
  anz[0] = cg[1] *cg[2];
  anz[1] = cg[2] *gcg[0];
  anz[2] = gcg[0]*gcg[1];
  ntot_send_cells = 2*(anz[0] + anz[1] + anz[2]);
  cell_start[0] = 0;
  for(i=1;i<6;i++) cell_start[i] = cell_start[i-1] + anz[(i-1)/2];

  /* create send/recv cell index lists */
  recv_cells = (int *)malloc(ntot_send_cells*sizeof(int));
  send_cells = (int *)malloc(ntot_send_cells*sizeof(int));
  /* direction loop (sorry, it looks nasty, and it is!!!). */
  for(i=0;i<3;i++) {
    lc[(i+1)%3] = 1-done[(i+1)%3]; hc[(i+1)%3] = cg[(i+1)%3]+done[(i+1)%3];
    lc[(i+2)%3] = 1-done[(i+2)%3]; hc[(i+2)%3] = cg[(i+2)%3]+done[(i+2)%3];
    /* send to :   left, down, for */
    lc[(i+0)%3] = 1;               hc[(i+0)%3] = 1;
    n_send_cells[(2*i)] = sub_grid_indices(send_cells, cell_start[(2*i)], 
					     ntot_send_cells, lc, hc, gcg);
     /* recv from : right, up, back */
    lc[(i+0)%3] = 0;               hc[(i+0)%3] = 0;
    n_recv_cells[(2*i)+1] = sub_grid_indices(recv_cells, cell_start[(2*i)+1], 
					     ntot_send_cells, lc, hc, gcg);
    /* send to :   right, up, back */
    lc[(i+0)%3] = cg[(i+0)%3];     hc[(i+0)%3] = cg[(i+0)%3];
    n_send_cells[(2*i)+1] = sub_grid_indices(send_cells, cell_start[(2*i)+1], 
					 ntot_send_cells, lc, hc, gcg);
    /* recv from : left, down, for */
    lc[(i+0)%3] = cg[(i+0)%3]+1;   hc[(i+0)%3] = cg[(i+0)%3]+1;
    n_recv_cells[(2*i)] = sub_grid_indices(recv_cells, cell_start[(2*i)], 
					 ntot_send_cells, lc, hc, gcg);
    done[i] = 1;
  }
  
  /* allocation of ghost cell information arrays */
  max_send_cells = 0;
  for(i=0;i<3;i++) if(anz[i]>max_send_cells) max_send_cells = anz[i];  
  n_send_ghosts = (int *)malloc((max_send_cells+1)*sizeof(int));
  n_recv_ghosts = (int *)malloc((max_send_cells+1)*sizeof(int));
  /* ghost exchange buffers */
  max_g_send_buf = max_g_recv_buf = PART_INCREMENT;
  g_send_buf = (Ghost *)malloc(max_g_send_buf*sizeof(Ghost));
  g_recv_buf = (Ghost *)malloc(max_g_recv_buf*sizeof(Ghost));

  /* init exchange forces/positions  */
  max_send_buf = max_recv_buf = PART_INCREMENT;
  send_buf       = (double *)malloc(3*max_send_buf*sizeof(double));
  recv_buf       = (double *)malloc(3*max_recv_buf*sizeof(double));
}

void exchange_part()
{
  int d, lr, dir, n;
  GHOST_TRACE(fprintf(stderr,"%d: exchange_part:\n",this_node));

  /* fold coordinates to primary simulation box */
  for(n=0;n<n_particles;n++) fold_particle(particles[n].p,particles[n].i);


  /* check part array */
#ifdef ADDITIONAL_CHECKS
  for(n=0;n<n_particles;n++) {
    if(particles[n].identity <0 || particles[n].identity > max_seen_particle) {
      fprintf(stderr,"%d: illigal identity %d of part %d\n",
	      this_node,particles[n].identity,n);
      errexit();
    }
    for(dir=0;dir<3;dir++) {
      if(periodic[dir] && (particles[n].p[dir] < 0 || particles[n].p[dir] > box_l[dir])) {
	fprintf(stderr,"%d: illigal position[%d] = %f of part %d\n",
		this_node,dir,particles[n].p[dir],n);
	errexit();
      }
    }
  }
#endif

  for(d=0; d<3; d++) {                           /* direction loop */  
    if(node_grid[d] > 1) { /* catch single node case for direction d! */
      for(lr=0; lr<2; lr++) {
	dir = 2*d + lr;
	n_p_send_buf = n_p_recv_buf = 0;
	n_b_send_buf = n_b_recv_buf = 0;
	if(lr==0) {                                /* left */
	  for(n=0;n<n_particles;n++)
	    if(particles[n].p[d] <   my_left[d] ) {
	      n = move_to_p_buf(n); 
	    }
	}
	else {                                     /* right */
	  for(n=0;n<n_particles;n++)                  
	    if(particles[n].p[d] >=  my_right[d]) {
	      n = move_to_p_buf(n);
	    }
	}
	send_particles(dir);
	append_particles();
      }
      MPI_Barrier(MPI_COMM_WORLD);
    }
  }

}

void exchange_ghost()
{
  int i,dir,c,n,m,c_ind, d;
  int old_max_send, old_max_recv;

  GHOST_TRACE(fprintf(stderr,"%d: exchange_ghost:\n",this_node));

  /* remove ghosts */
  for(n=n_particles; n<n_particles+n_ghosts; n++) local_index[particles[n].identity] = -1;

  old_max_send = max_send_buf;
  old_max_recv = max_recv_buf;
  n_ghosts = 0;
  for(dir=0; dir<6; dir++) {                      /* direction loop */
    n_send_ghosts[max_send_cells] = 0;
    n_g_send_buf = 0;
    /* send cell loop - number of ghosts to send */
    for(c=0; c<n_send_cells[dir]; c++) {
      c_ind = send_cells[cell_start[dir]+c];
      n_send_ghosts[c]               = cells[c_ind].n_particles;
      n_send_ghosts[max_send_cells] += cells[c_ind].n_particles;
    }
    /* check buffer size */
    if(n_send_ghosts[max_send_cells] > max_g_send_buf) {
      max_g_send_buf = n_send_ghosts[max_send_cells];
      g_send_buf = (Ghost *)realloc(g_send_buf,max_g_send_buf*sizeof(Ghost));
    }
    /* send cell loop - pack ghosts */
    for(c=0; c<n_send_cells[dir]; c++) {
      c_ind = send_cells[cell_start[dir]+c];
      for(n=0;n<cells[c_ind].n_particles;n++) {
	if(n_g_send_buf > max_g_send_buf) {
	  fprintf(stderr,"%d: g_send_buf overflow (%d > %d)\n",
		  this_node,max_g_send_buf,n_g_send_buf);
	}
	pack_ghost(g_send_buf,n_g_send_buf,particles,cells[c_ind].particles[n]);
	n_g_send_buf++;
      }
    }
    /* fold ghost coordinates if they cross the boundary */
    switch(boundary[dir]) {
    case 0: break;
    case 1:
      d=dir/2;
      for(n=0;n<n_g_send_buf;n++) g_send_buf[n].p[d] += box_l[dir];
      break;
    case -1:
      d=dir/2;
      for(n=0;n<n_g_send_buf;n++) g_send_buf[n].p[d] -= box_l[dir];
      break;
    }

    /* send ghosts */
    send_ghosts(dir);
    /* sort recieved ghosts into cells */

    /* index of first ghost of that direction in local particle array */
    m = n_particles + n_ghosts;
    /* check particle array */
    n_ghosts += n_recv_ghosts[max_send_cells];
    if(n_particles+n_ghosts >= max_particles) 
      realloc_particles(n_particles+n_ghosts);
    /* index of actual ghost in g_recv_buf */
    n=0;
    for(c=0; c<n_recv_cells[dir]; c++) {   /* cell loop - unpack ghosts */
      c_ind = recv_cells[cell_start[dir] + c];
      if(n_recv_ghosts[c] >= cells[c_ind].max_particles) {
	realloc_cell_particles(c_ind, n_recv_ghosts[c]+1);
      }
      for(i=0; i<n_recv_ghosts[c];i++) {
	unpack_ghost(particles,m,g_recv_buf,n);
	local_index[particles[m].identity] = m;
	cells[c_ind].particles[i] = m;
	GHOST_TRACE(fprintf(stderr,"%d: Ghost %d at (%.2f, %.2f, %.2f) sorted in cell %d as No. %d\n",this_node,m,particles[m].p[0],particles[m].p[1],particles[m].p[2],c_ind,i));
	m++; n++; 
      }
      cells[c_ind].n_particles = n_recv_ghosts[c];
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }                                               /* END direction loop */


  /* resize pos/force buffer if necessary */
  if(max_send_buf > max_recv_buf) max_recv_buf = max_send_buf;
  else                            max_send_buf = max_recv_buf;
  if(max_send_buf > old_max_send)
    send_buf = (double *)realloc(send_buf,3*max_send_buf*sizeof(double));
  if(max_recv_buf > old_max_recv)    
    recv_buf = (double *)realloc(recv_buf,3*max_recv_buf*sizeof(double));
}

void update_ghost_pos()
{
  int dir,c,c_ind,n,g,i;
  GHOST_TRACE(fprintf(stderr,"%d: update_ghost_pos:\n",this_node));
  MPI_Barrier(MPI_COMM_WORLD);

  for(dir=0; dir<6; dir++) {          /* direction loop forward */
    /* loop send cells -  copy positions to buffer*/
    g=0;
    for(c=0; c<n_send_cells[dir]; c++) {
      c_ind = send_cells[cell_start[dir]+c];
      for(n=0; n<cells[c_ind].n_particles; n++) {
	for(i=0;i<3;i++) send_buf[g+i] = particles[cells[c_ind].particles[n]].p[i];
	g += 3;
      }
    }
    /* fold positions if they cross the boundary */
    switch(boundary[dir]) {
    case 0: break;
    case 1:
      for(n=dir/2; n<g; n+=3) send_buf[n] += box_l[dir];
      break;
    case -1:
      for(n=dir/2; n<g; n+=3) send_buf[n] -= box_l[dir];
      break;
    }
    /* send buffer */
    send_posforce(dir,3*ghost_send_size[dir],3*ghost_recv_size[dir]);

    /* loop recv cells - copy positions from buffer */
    g=0;
    for(c=0; c<n_recv_cells[dir]; c++) {
      c_ind = recv_cells[cell_start[dir]+c];
      for(n=0; n<cells[c_ind].n_particles; n++) {
	for(i=0;i<3;i++) particles[cells[c_ind].particles[n]].p[i] = recv_buf[g+i];
	g += 3;
      }
    }
    MPI_Barrier(MPI_COMM_WORLD); 
  }
}

void collect_ghost_forces()
{
  int dir,c,c_ind,n,g,i;
  GHOST_TRACE(fprintf(stderr,"%d: collect_ghost_forces:\n",this_node));
  MPI_Barrier(MPI_COMM_WORLD);

  for(dir=5; dir>=0; dir--) {         /* direction loop backward */
    /* loop recv cells - copy forces to buffer*/
    g=0;
    for(c=0; c<n_recv_cells[dir]; c++) {
      c_ind = recv_cells[cell_start[dir]+c];
      for(n=0; n<cells[c_ind].n_particles; n++) {
	for(i=0;i<3;i++) {
	  send_buf[g+i] = particles[cells[c_ind].particles[n]].f[i];
	}
	g += 3;
      }
    }
    /* send forces */
    send_posforce(dir,3*ghost_recv_size[dir],3*ghost_send_size[dir]);
    /* loop send cells - add buffer forces to local forces */
    g=0;
    for(c=0; c<n_send_cells[dir]; c++) {
      c_ind = send_cells[cell_start[dir]+c];
      for(n=0;n<cells[c_ind].n_particles;n++) {  
	for(i=0; i<3; i++) {
	  particles[cells[c_ind].particles[n]].f[i] += recv_buf[g+i];
	}
	g += 3;
      }
    } 
    MPI_Barrier(MPI_COMM_WORLD); 
  }
  MPI_Barrier(MPI_COMM_WORLD); 
}

void ghost_exit()
{
  GHOST_TRACE(fprintf(stderr,"%d: ghost_exit:\n",this_node));
  if(max_p_send_buf>0)  free(p_send_buf);
  if(max_p_recv_buf>0)  free(p_recv_buf);
  if(max_b_send_buf>0)  free(b_send_buf);
  if(max_b_recv_buf>0)  free(b_recv_buf);
  if(ntot_send_cells>0) free(send_cells);
  if(ntot_send_cells>0) free(recv_cells);
  if(max_send_cells>0)  free(n_send_ghosts);
  if(max_send_cells>0)  free(n_recv_ghosts);
  if(max_g_send_buf>0)  free(g_send_buf);
  if(max_g_recv_buf>0)  free(g_recv_buf);
  if(max_send_buf>0)    free(send_buf);
  if(max_recv_buf>0)    free(recv_buf);
}

/*******************  privat functions  *******************/

/** Creates an linear index list of a sub grid.
    The sub grid is defined by its lower and upper corner:\\
    from (lc[0],lc[1],lc[2]) to (hc[0],hc[1],hc[2])\\
    The grid dimension is given with (gs[0],gs[1],gs[2])\\
    The linear index list of length <returnvalue> is stored 
    in list starting from position start. max should be the 
    total length of list to ensure that the indices fit into list.
 */  
int sub_grid_indices(int* list, int start, int max,
		     int lc[3], int hc[3], int gs[3])
{
  int i,size,p0,p1,p2;
  /* sanity check */
  for(i=0;i<3;i++) {
    if(lc[i]<0 || lc[i] >= gs[i]) return 0;
    if(hc[i]<0 || hc[i] >= gs[i]) return 0;
    if(lc[i] > hc[i]) return 0;
  }

  size = (hc[0]+1-lc[0]);
  for(i=1;i<3;i++) size *= (hc[i]+1-lc[i]);
  /* check array size */
  if(size+start>max) return -1;

  i=start;
  for(p0=lc[0];p0<=hc[0];p0++)
    for(p1=lc[1];p1<=hc[1];p1++)
      for(p2=lc[2];p2<=hc[2];p2++) {
	if(i>max) fprintf(stderr,"%d: sub_grid_indices: Array overflow: %d>%d\n",this_node,i,max); 
	list[i] = get_linear_index(p0,p1,p2,gs);
	i++;
      }

  return size;
}

int move_to_p_buf(int ind)
{
  int bonds=0,i;
  /* count bond partners */
  for(i=0;i<particles[ind].n_bonds;i++) {
    /* ATTENTION: HERE YOU HAVE TO CHECK THE NUMBER OF BOND PARTNERS FOR EACH BOND ! */
    bonds++;
  }
    
  /* check buffer sizes */
  if(n_p_send_buf >= max_p_send_buf) {
    max_p_send_buf += PART_INCREMENT;
    p_send_buf = (Particle *)realloc(p_send_buf, max_p_send_buf*sizeof(Particle));
  }
  if(bonds>0 && n_b_send_buf+bonds >= max_b_send_buf) {
    max_b_send_buf = n_b_send_buf + bonds +  PART_INCREMENT;
    b_send_buf = (int *)realloc(b_send_buf, max_b_send_buf*sizeof(int));
  }

  /* copy particle from local array to send buffer */
  if(bonds>0) {
    memcpy(&b_send_buf[n_b_send_buf], particles[ind].bonds, bonds*sizeof(int));
    n_b_send_buf += bonds;
    free(particles[ind].bonds);
  }
  memcpy(&p_send_buf[n_p_send_buf], &particles[ind], sizeof(Particle));
  n_p_send_buf++;
  /* delete particle = 
     1: update the local_index list 
     2: move last particle (n_particles-1) to the free position (ind)
        ... if there is anything to move */
  local_index[particles[ind].identity] = -1;
  if(ind < n_particles-1) {
    local_index[particles[n_particles-1].identity] = ind;
    memcpy(&particles[ind],&particles[n_particles-1],sizeof(Particle));
    /* decrease return value in order to 
       continue the particle loop at the right place */
    ind--;
  }
  n_particles--;
  return ind;
}

void send_particles(int s_dir)
{
  int evenodd;
  int r_dir;
  int send_sizes[2],recv_sizes[2];
  MPI_Status status;

  /* check if communication goes to the very same node */
  if(node_neighbors[s_dir] != this_node) {
    send_sizes[0]=n_p_send_buf;
    send_sizes[1]=n_b_send_buf;
    /* calc recv direction (r_dir) from send direction (s_dir) */
    if(s_dir%2 == 0) r_dir = s_dir+1;
    else             r_dir = s_dir-1;
    /* two step communication: first all even positions than all odd */
    for(evenodd=0; evenodd<2;evenodd++) {
      if((node_pos[s_dir/2]+evenodd)%2==0) {
	MPI_Send(send_sizes,2,MPI_INT,node_neighbors[s_dir],0,MPI_COMM_WORLD);
	MPI_Send(p_send_buf,n_p_send_buf*sizeof(Particle), MPI_BYTE, 
		 node_neighbors[s_dir],0,MPI_COMM_WORLD);
	if(n_b_send_buf>0)
	  MPI_Send(b_send_buf, n_b_send_buf, MPI_INT, node_neighbors[s_dir], 0, MPI_COMM_WORLD);
      }
      else {
	MPI_Recv(recv_sizes,2,MPI_INT,node_neighbors[r_dir],0,MPI_COMM_WORLD,&status);
	n_p_recv_buf = recv_sizes[0];
	n_b_recv_buf = recv_sizes[1];
	if(n_p_recv_buf >= max_p_recv_buf) {
	  max_p_recv_buf = n_p_recv_buf;
	  p_recv_buf = (Particle *)realloc(p_recv_buf,max_p_recv_buf*sizeof(Particle));
	}
	if(n_b_recv_buf >= max_b_recv_buf) {
	  max_b_recv_buf = n_b_recv_buf;
	  b_recv_buf = (int *)realloc(b_recv_buf,max_b_recv_buf*sizeof(int));
	}
	MPI_Recv(p_recv_buf,n_p_recv_buf*sizeof(Particle), MPI_BYTE,
		 node_neighbors[r_dir],0,MPI_COMM_WORLD,&status);
	if(n_b_recv_buf>0) 
	  MPI_Recv(b_recv_buf,n_b_recv_buf, MPI_INT, node_neighbors[r_dir],0,MPI_COMM_WORLD,&status);
      }
    }
  }
  else {                 /* communication goes to the same node! */ 
    int tmp_size;
    Particle *tmp_pp;
    int *tmp_ip;
    /* particle buffer */
    tmp_size = n_p_send_buf; 
    n_p_send_buf = n_p_recv_buf; n_p_recv_buf = tmp_size;
    tmp_size = max_p_send_buf; 
    max_p_send_buf = max_p_recv_buf; max_p_recv_buf = tmp_size;
    tmp_pp = p_send_buf;
    p_send_buf = p_recv_buf; p_recv_buf = tmp_pp;
    /* bond buffer */
    tmp_size = n_b_send_buf; 
    n_b_send_buf = n_b_recv_buf; n_b_recv_buf = tmp_size;
    tmp_size = max_b_send_buf; 
    max_b_send_buf = max_b_recv_buf; max_b_recv_buf = tmp_size;
    tmp_ip = b_send_buf;
    b_send_buf = b_recv_buf; b_recv_buf = tmp_ip;
  }
}

void append_particles(void)
{
  int n,i,bonds,b_ind=0;
  /* check particle array size */
  if(n_particles + n_p_recv_buf >= max_particles) 
    realloc_particles(n_particles + n_p_recv_buf);
  /* copy particle data */
  memcpy(&particles[n_particles],p_recv_buf,n_p_recv_buf*sizeof(Particle));
  /* update local_index list */
  for(n=n_particles; n < n_particles+n_p_recv_buf; n++) {
    local_index[particles[n].identity] = n;
    bonds=0;
    /* ATTENTION: HERE YOU HAVE TO CHECK THE NUMBER OF BOND PARTNERS FOR EACH BOND ! */
    for(i=0; i<particles[n].n_bonds; i++) { bonds++; }
    if(bonds>0) {
      particles[n].bonds = (int *)malloc(bonds*sizeof(int));
      memcpy(particles[n].bonds, &(b_recv_buf[b_ind]), bonds*sizeof(int));
      b_ind += bonds;
    }
  }
  /* update particle number */
  n_particles += n_p_recv_buf;
}

void send_ghosts(int s_dir)
{
  int evenodd;
  int r_dir;
  MPI_Status status;


  /* check if communication goes to the very same node */
  if(node_neighbors[s_dir] != this_node) {
    /* calc recv direction (r_dir) from send direction (s_dir) */
    if(s_dir%2 == 0) r_dir = s_dir+1;
    else             r_dir = s_dir-1;
    /* two step communication: first all even positions than all odd */
    for(evenodd=0; evenodd<2;evenodd++) {
      if((node_pos[s_dir/2]+evenodd)%2==0) {
	MPI_Send(n_send_ghosts, max_send_cells+1, MPI_INT,
		 node_neighbors[s_dir], 0, MPI_COMM_WORLD);
	MPI_Send(g_send_buf,n_send_ghosts[max_send_cells]*sizeof(Ghost), MPI_BYTE, 
		 node_neighbors[s_dir],0,MPI_COMM_WORLD);
      }
      else {
	MPI_Recv(n_recv_ghosts, max_send_cells+1, MPI_INT,
		 node_neighbors[r_dir],0,MPI_COMM_WORLD,&status);
	if(n_recv_ghosts[max_send_cells] >= max_g_recv_buf) {
	  max_g_recv_buf = n_recv_ghosts[max_send_cells];
	  g_recv_buf = (Ghost *)realloc(g_recv_buf,max_g_recv_buf*sizeof(Ghost));
	}
	MPI_Recv(g_recv_buf,n_recv_ghosts[max_send_cells]*sizeof(Ghost), MPI_BYTE,
		 node_neighbors[r_dir],0,MPI_COMM_WORLD,&status);
      }
    }
  }
  else {                 /* communication goes to the same node! */ 
    int tmp;
    Ghost *tmp_gp;

    memcpy(n_recv_ghosts,n_send_ghosts,(max_send_cells+1)*sizeof(int));
    tmp_gp = g_send_buf;
    g_send_buf = g_recv_buf; g_recv_buf = tmp_gp;
    tmp = max_g_send_buf;
    max_g_send_buf = max_g_recv_buf;
    max_g_recv_buf = tmp;
    tmp = n_g_send_buf;
    n_g_send_buf = n_g_recv_buf;
    n_g_recv_buf = tmp;
  }

  /* store number of ghosts to send in direction s_dir */
  ghost_send_size[s_dir] = n_send_ghosts[max_send_cells];
  if(ghost_send_size[s_dir]>max_send_buf) 
    max_send_buf = ghost_send_size[s_dir];
  /* store number of ghosts to recieve from direction s_dir */
  ghost_recv_size[s_dir] = n_recv_ghosts[max_send_cells];
  if(ghost_recv_size[s_dir]>max_recv_buf) 
    max_recv_buf = ghost_recv_size[s_dir];
}

void pack_ghost(Ghost *g_array, int g_ind, Particle *p_array, int p_ind)
{
  g_array[g_ind].identity = p_array[p_ind].identity ;
  g_array[g_ind].type     = p_array[p_ind].type ;
  memcpy(g_array[g_ind].p, p_array[p_ind].p, 3*sizeof(double));
  g_array[g_ind].q        = p_array[p_ind].q ;
}

void unpack_ghost(Particle *p_array, int p_ind, Ghost *g_array, int g_ind)
{
  p_array[p_ind].identity = g_array[g_ind].identity ;
  p_array[p_ind].type     = g_array[g_ind].type ;
  memcpy(p_array[p_ind].p, g_array[g_ind].p, 3*sizeof(double));
  p_array[p_ind].q        = g_array[g_ind].q ;
}

void send_posforce(int s_dir,int send_size, int recv_size)
{
  int evenodd;
  int r_dir;
  MPI_Status status;

  /* check if communication goes to the very same node */
  if(node_neighbors[s_dir] != this_node) {
    /* calc recv direction (r_dir) from send direction (s_dir) */
    if(s_dir%2 == 0) r_dir = s_dir+1;
    else             r_dir = s_dir-1;
    /* two step communication: first all even positions than all odd */
    for(evenodd=0; evenodd<2;evenodd++) {
      if((node_pos[s_dir/2]+evenodd)%2==0) 
	MPI_Send(send_buf, send_size, MPI_DOUBLE, 
		 node_neighbors[s_dir],0,MPI_COMM_WORLD);
      else 
	MPI_Recv(recv_buf, recv_size, MPI_DOUBLE,
		 node_neighbors[r_dir],0,MPI_COMM_WORLD,&status);    
    }
  }
  else {                  /* communication goes to the same node! */ 
    double *tmp_dp;
    tmp_dp = send_buf;
    send_buf = recv_buf; recv_buf = tmp_dp;        
  }
}

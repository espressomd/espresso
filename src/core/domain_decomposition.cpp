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
/** \file domain_decomposition.cpp
 *
 *  This file contains everything related to the cell system: domain decomposition.
 *  See also \ref domain_decomposition.hpp
 */

#include "domain_decomposition.hpp"
#include "errorhandling.hpp"
#include "forces.hpp"
#include "pressure.hpp"
#include "energy.hpp"
#include "constraint.hpp"
#include "initialize.hpp"
#include "external_potential.hpp"

/************************************************/
/** \name Defines */
/************************************************/
/*@{*/

/** half the number of cell neighbors in 3 Dimensions. */
#define CELLS_MAX_NEIGHBORS 14

/*@}*/

/************************************************/
/** \name Variables */
/************************************************/
/*@{*/

DomainDecomposition dd = { 1, {0,0,0}, {0,0,0}, {0,0,0}, {0,0,0}, NULL };

int max_num_cells = CELLS_MAX_NUM_CELLS;
int min_num_cells = 1;
double max_skin   = 0.0;

/*@}*/

/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Convenient replace for loops over all cells. */
#define DD_CELLS_LOOP(m,n,o) \
  for(o=0; o<dd.ghost_cell_grid[2]; o++) \
    for(n=0; n<dd.ghost_cell_grid[1]; n++) \
      for(m=0; m<dd.ghost_cell_grid[0]; m++)

/** Convenient replace for loops over Local cells. */
#define DD_LOCAL_CELLS_LOOP(m,n,o) \
  for(o=1; o<dd.cell_grid[2]+1; o++) \
    for(n=1; n<dd.cell_grid[1]+1; n++) \
      for(m=1; m<dd.cell_grid[0]+1; m++)

/** Convenient replace for inner cell check. usage: if(DD_IS_LOCAL_CELL(m,n,o)) {...} */
#define DD_IS_LOCAL_CELL(m,n,o) \
  ( m > 0 && m < dd.ghost_cell_grid[0] - 1 && \
    n > 0 && n < dd.ghost_cell_grid[1] - 1 && \
    o > 0 && o < dd.ghost_cell_grid[2] - 1 ) 

/** Convenient replace for ghost cell check. usage: if(DD_IS_GHOST_CELL(m,n,o)) {...} */
#define DD_IS_GHOST_CELL(m,n,o) \
  ( m == 0 || m == dd.ghost_cell_grid[0] - 1 || \
    n == 0 || n == dd.ghost_cell_grid[1] - 1 || \
    o == 0 || o == dd.ghost_cell_grid[2] - 1 ) 

/** Calculate cell grid dimensions, cell sizes and number of cells.
 *  Calculates the cell grid, based on \ref local_box_l and \ref
 *  max_range. If the number of cells is larger than \ref
 *  max_num_cells, it increases max_range until the number of cells is
 *  smaller or equal \ref max_num_cells. It sets: \ref
 *  DomainDecomposition::cell_grid, \ref
 *  DomainDecomposition::ghost_cell_grid, \ref
 *  DomainDecomposition::cell_size, \ref
 *  DomainDecomposition::inv_cell_size, and \ref n_cells.
 */
void dd_create_cell_grid()
{
  int i,n_local_cells,new_cells,min_ind;
  double cell_range[3], min_size, scale, volume;
  CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid: max_range %f\n",this_node,max_range));
  CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid: local_box %f-%f, %f-%f, %f-%f,\n",this_node,my_left[0],my_right[0],my_left[1],my_right[1],my_left[2],my_right[2]));
  
  /* initialize */
  cell_range[0]=cell_range[1]=cell_range[2] = max_range;

  if (max_range < ROUND_ERROR_PREC*box_l[0]) {
    /* this is the initialization case */
    n_local_cells = dd.cell_grid[0] = dd.cell_grid[1] = dd.cell_grid[2]=1;
  }
  else {
    /* Calculate initial cell grid */
    volume = local_box_l[0];
    for(i=1;i<3;i++) volume *= local_box_l[i];
    scale = pow(max_num_cells/volume, 1./3.);
    for(i=0;i<3;i++) {
      /* this is at least 1 */
      dd.cell_grid[i] = (int)ceil(local_box_l[i]*scale);
      cell_range[i] = local_box_l[i]/dd.cell_grid[i];

      if ( cell_range[i] < max_range ) {
	/* ok, too many cells for this direction, set to minimum */
	dd.cell_grid[i] = (int)floor(local_box_l[i]/max_range);
	if ( dd.cell_grid[i] < 1 ) {
	  char *error_msg = runtime_error(ES_INTEGER_SPACE + 2*ES_DOUBLE_SPACE + 128);
	  ERROR_SPRINTF(error_msg, "{002 interaction range %g in direction %d is larger than the local box size %g} ",
			max_range, i, local_box_l[i]);
	  dd.cell_grid[i] = 1;
	}
	cell_range[i] = local_box_l[i]/dd.cell_grid[i];
      }
    }

    /* It may be necessary to asymmetrically assign the scaling to the coordinates, which the above approach will not do.
       For a symmetric box, it gives a symmetric result. Here we correct that. */
    for (;;) {
      n_local_cells = dd.cell_grid[0];
      for (i = 1; i < 3; i++)
	n_local_cells *= dd.cell_grid[i];

      /* done */
      if (n_local_cells <= max_num_cells)
	break;

      /* find coordinate with the smallest cell range */
      min_ind = 0;
      min_size = cell_range[0];
      for (i = 1; i < 3; i++)
	if (dd.cell_grid[i] > 1 && cell_range[i] < min_size) {
	  min_ind = i;
	  min_size = cell_range[i];
	}
      CELL_TRACE(fprintf(stderr, "%d: minimal coordinate %d, size %f, grid %d\n", this_node,min_ind, min_size, dd.cell_grid[min_ind]));

      dd.cell_grid[min_ind]--;
      cell_range[min_ind] = local_box_l[min_ind]/dd.cell_grid[min_ind];
    }
    CELL_TRACE(fprintf(stderr, "%d: final %d %d %d\n", this_node, dd.cell_grid[0], dd.cell_grid[1], dd.cell_grid[2]));

    /* sanity check */
    if (n_local_cells < min_num_cells) {
      char *error_msg = runtime_error(ES_INTEGER_SPACE + 2*ES_DOUBLE_SPACE + 128);
      ERROR_SPRINTF(error_msg, "{001 number of cells %d is smaller than minimum %d (interaction range too large or min_num_cells too large)} ",
		    n_local_cells, min_num_cells);
    }
  }

  /* quit program if unsuccesful */
  if(n_local_cells > max_num_cells) {
    char *error_msg = runtime_error(128);
    ERROR_SPRINTF(error_msg, "{003 no suitable cell grid found} ");
  }

  /* now set all dependent variables */
  new_cells=1;
  for(i=0;i<3;i++) {
    dd.ghost_cell_grid[i] = dd.cell_grid[i]+2;	
    new_cells              *= dd.ghost_cell_grid[i];
    dd.cell_size[i]       = local_box_l[i]/(double)dd.cell_grid[i];
    dd.inv_cell_size[i]   = 1.0 / dd.cell_size[i];
  }
  max_skin = dmin(dmin(dd.cell_size[0],dd.cell_size[1]),dd.cell_size[2]) - max_cut;

  /* allocate cell array and cell pointer arrays */
  realloc_cells(new_cells);
  realloc_cellplist(&local_cells, local_cells.n = n_local_cells);
  realloc_cellplist(&ghost_cells, ghost_cells.n = new_cells-n_local_cells);

  CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid, n_cells=%d, local_cells.n=%d, ghost_cells.n=%d, dd.ghost_cell_grid=(%d,%d,%d)\n", this_node, n_cells,local_cells.n,ghost_cells.n,dd.ghost_cell_grid[0],dd.ghost_cell_grid[1],dd.ghost_cell_grid[2]));
}

/** Fill local_cells list and ghost_cells list for use with domain
    decomposition.  \ref cells::cells is assumed to be a 3d grid with size
    \ref DomainDecomposition::ghost_cell_grid . */
void dd_mark_cells()
{
  int m,n,o,cnt_c=0,cnt_l=0,cnt_g=0;
  DD_CELLS_LOOP(m,n,o) {
    if(DD_IS_LOCAL_CELL(m,n,o)) local_cells.cell[cnt_l++] = &cells[cnt_c++]; 
    else                        ghost_cells.cell[cnt_g++] = &cells[cnt_c++];
  } 
}

/** Fill a communication cell pointer list. Fill the cell pointers of
    all cells which are inside a rectangular subgrid of the 3D cell
    grid (\ref DomainDecomposition::ghost_cell_grid) starting from the
    lower left corner lc up to the high top corner hc. The cell
    pointer list part_lists must already be large enough.
    \param part_lists  List of cell pointers to store the result.
    \param lc          lower left corner of the subgrid.
    \param hc          high up corner of the subgrid.
 */
int dd_fill_comm_cell_lists(Cell **part_lists, int lc[3], int hc[3])
{
  int i,m,n,o,c=0;
  /* sanity check */
  for(i=0; i<3; i++) {
    if(lc[i]<0 || lc[i] >= dd.ghost_cell_grid[i]) return 0;
    if(hc[i]<0 || hc[i] >= dd.ghost_cell_grid[i]) return 0;
    if(lc[i] > hc[i]) return 0;
  }

  for(o=lc[0]; o<=hc[0]; o++) 
    for(n=lc[1]; n<=hc[1]; n++) 
      for(m=lc[2]; m<=hc[2]; m++) {
	i = get_linear_index(o,n,m,dd.ghost_cell_grid);
	CELL_TRACE(fprintf(stderr,"%d: dd_fill_comm_cell_list: add cell %d\n",this_node,i));
	part_lists[c] = &cells[i];
	c++;
      }
  return c;
}

/** Create communicators for cell structure domain decomposition. (see \ref GhostCommunicator) */
void  dd_prepare_comm(GhostCommunicator *comm, int data_parts)
{
  int dir,lr,i,cnt, num, n_comm_cells[3];
  int lc[3],hc[3],done[3]={0,0,0};

  /* calculate number of communications */
  num = 0;
  for(dir=0; dir<3; dir++) { 
    for(lr=0; lr<2; lr++) {
#ifdef PARTIAL_PERIODIC
      /* No communication for border of non periodic direction */
      if( PERIODIC(dir) || (boundary[2*dir+lr] == 0) ) 
#endif
	{
	  if(node_grid[dir] == 1 ) num++;
	  else num += 2;
	}
    }
  }

  /* prepare communicator */
  CELL_TRACE(fprintf(stderr,"%d Create Communicator: prep_comm data_parts %d num %d\n",this_node,data_parts,num));
  prepare_comm(comm, data_parts, num);

  /* number of cells to communicate in a direction */
  n_comm_cells[0] = dd.cell_grid[1]       * dd.cell_grid[2];
  n_comm_cells[1] = dd.cell_grid[2]       * dd.ghost_cell_grid[0];
  n_comm_cells[2] = dd.ghost_cell_grid[0] * dd.ghost_cell_grid[1];

  cnt=0;
  /* direction loop: x, y, z */
  for(dir=0; dir<3; dir++) {
    lc[(dir+1)%3] = 1-done[(dir+1)%3]; 
    lc[(dir+2)%3] = 1-done[(dir+2)%3];
    hc[(dir+1)%3] = dd.cell_grid[(dir+1)%3]+done[(dir+1)%3];
    hc[(dir+2)%3] = dd.cell_grid[(dir+2)%3]+done[(dir+2)%3];
    /* lr loop: left right */
    /* here we could in principle build in a one sided ghost
       communication, simply by taking the lr loop only over one
       value */
    for(lr=0; lr<2; lr++) {
      if(node_grid[dir] == 1) {
	/* just copy cells on a single node */
#ifdef PARTIAL_PERIODIC
	if( PERIODIC(dir ) || (boundary[2*dir+lr] == 0) ) 
#endif
	  {
	    comm->comm[cnt].type          = GHOST_LOCL;
	    comm->comm[cnt].node          = this_node;
	    /* Buffer has to contain Send and Recv cells -> factor 2 */
	    comm->comm[cnt].part_lists    = (ParticleList**)malloc(2*n_comm_cells[dir]*sizeof(ParticleList *));
	    comm->comm[cnt].n_part_lists  = 2*n_comm_cells[dir];
	    /* prepare folding of ghost positions */
	    if((data_parts & GHOSTTRANS_POSSHFTD) && boundary[2*dir+lr] != 0) 
	      comm->comm[cnt].shift[dir] = boundary[2*dir+lr]*box_l[dir];
	    /* fill send comm cells */
	    lc[(dir+0)%3] = hc[(dir+0)%3] = 1+lr*(dd.cell_grid[(dir+0)%3]-1);  
	    dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
	    CELL_TRACE(fprintf(stderr,"%d: prep_comm %d copy to          grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,
			       lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	    /* fill recv comm cells */
	    lc[(dir+0)%3] = hc[(dir+0)%3] = 0+(1-lr)*(dd.cell_grid[(dir+0)%3]+1);
	    /* place recieve cells after send cells */
	    dd_fill_comm_cell_lists(&comm->comm[cnt].part_lists[n_comm_cells[dir]],lc,hc);
	    CELL_TRACE(fprintf(stderr,"%d: prep_comm %d copy from        grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	    cnt++;
	  }
      }
      else {
	/* i: send/recv loop */
	for(i=0; i<2; i++) {  
#ifdef PARTIAL_PERIODIC
	  if( PERIODIC(dir) || (boundary[2*dir+lr] == 0) ) 
#endif
	    if((node_pos[dir]+i)%2==0) {
	      comm->comm[cnt].type          = GHOST_SEND;
	      comm->comm[cnt].node          = node_neighbors[2*dir+lr];
	      comm->comm[cnt].part_lists    = (ParticleList**)malloc(n_comm_cells[dir]*sizeof(ParticleList *));
	      comm->comm[cnt].n_part_lists  = n_comm_cells[dir];
	      /* prepare folding of ghost positions */
	      if((data_parts & GHOSTTRANS_POSSHFTD) && boundary[2*dir+lr] != 0) 
		comm->comm[cnt].shift[dir] = boundary[2*dir+lr]*box_l[dir];
	      
	      lc[(dir+0)%3] = hc[(dir+0)%3] = 1+lr*(dd.cell_grid[(dir+0)%3]-1);  
	      dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
	      
	      CELL_TRACE(fprintf(stderr,"%d: prep_comm %d send to   node %d grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,
				 comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	      cnt++;
	    }
#ifdef PARTIAL_PERIODIC
	  if( PERIODIC(dir) || (boundary[2*dir+(1-lr)] == 0) ) 
#endif
	    if((node_pos[dir]+(1-i))%2==0) {
	      comm->comm[cnt].type          = GHOST_RECV;
	      comm->comm[cnt].node          = node_neighbors[2*dir+(1-lr)];
	      comm->comm[cnt].part_lists    = (ParticleList**)malloc(n_comm_cells[dir]*sizeof(ParticleList *));
	      comm->comm[cnt].n_part_lists  = n_comm_cells[dir];
	      
	      lc[(dir+0)%3] = hc[(dir+0)%3] = 0+(1-lr)*(dd.cell_grid[(dir+0)%3]+1);
	      dd_fill_comm_cell_lists(comm->comm[cnt].part_lists,lc,hc);
	      
	      CELL_TRACE(fprintf(stderr,"%d: prep_comm %d recv from node %d grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,
				 comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]));
	      cnt++;
	    }
	}
      }
      done[dir]=1;
    }
  }
}

/** Revert the order of a communicator: After calling this the
    communicator is working in reverted order with exchanged
    communication types GHOST_SEND <-> GHOST_RECV. */
void dd_revert_comm_order(GhostCommunicator *comm)
{
  int i,j,nlist2;
  GhostCommunication tmp;
  ParticleList *tmplist;

  CELL_TRACE(fprintf(stderr,"%d: dd_revert_comm_order: anz comm: %d\n",this_node,comm->num));

  /* revert order */
  for(i=0; i<(comm->num/2); i++) {
    tmp = comm->comm[i];
    comm->comm[i] = comm->comm[comm->num-i-1];
    comm->comm[comm->num-i-1] = tmp;
  }
  /* exchange SEND/RECV */
  for(i=0; i<comm->num; i++) {
    if(comm->comm[i].type == GHOST_SEND) comm->comm[i].type = GHOST_RECV;
    else if(comm->comm[i].type == GHOST_RECV) comm->comm[i].type = GHOST_SEND;
    else if(comm->comm[i].type == GHOST_LOCL) {
      nlist2=comm->comm[i].n_part_lists/2;
      for(j=0;j<nlist2;j++) {
	tmplist = comm->comm[i].part_lists[j];
	comm->comm[i].part_lists[j] = comm->comm[i].part_lists[j+nlist2];
	comm->comm[i].part_lists[j+nlist2] = tmplist;
      }
    }
  }
}

/** Of every two communication rounds, set the first receivers to prefetch and poststore */
void dd_assign_prefetches(GhostCommunicator *comm)
{
  int cnt;

  for(cnt=0; cnt<comm->num; cnt += 2) {
    if (comm->comm[cnt].type == GHOST_RECV && comm->comm[cnt + 1].type == GHOST_SEND) {
      comm->comm[cnt].type |= GHOST_PREFETCH | GHOST_PSTSTORE;
      comm->comm[cnt + 1].type |= GHOST_PREFETCH | GHOST_PSTSTORE;
    }
  }
}

/** update the 'shift' member of those GhostCommunicators, which use
    that value to speed up the folding process of its ghost members
    (see \ref dd_prepare_comm for the original), i.e. all which have
    GHOSTTRANS_POSSHFTD or'd into 'data_parts' upon execution of \ref
    dd_prepare_comm. */
void dd_update_communicators_w_boxl()
{
  int cnt=0;
  /* direction loop: x, y, z */
  for(int dir=0; dir<3; dir++) {
    /* lr loop: left right */
    for(int lr=0; lr<2; lr++) {
      if(node_grid[dir] == 1) {
#ifdef PARTIAL_PERIODIC
	if( PERIODIC(dir ) || (boundary[2*dir+lr] == 0) ) 
#endif
	  {
	    /* prepare folding of ghost positions */
	    if(boundary[2*dir+lr] != 0) {
	      cell_structure.exchange_ghosts_comm.comm[cnt].shift[dir] = boundary[2*dir+lr]*box_l[dir]; 
	      cell_structure.update_ghost_pos_comm.comm[cnt].shift[dir] = boundary[2*dir+lr]*box_l[dir]; 
	    }
	    cnt++;
	  }
      }
      else {
	/* i: send/recv loop */
	for(int i=0; i<2; i++) {  
#ifdef PARTIAL_PERIODIC
	  if( PERIODIC(dir) || (boundary[2*dir+lr] == 0) ) 
#endif
	    if((node_pos[dir]+i)%2==0) {
	      /* prepare folding of ghost positions */
	      if(boundary[2*dir+lr] != 0) {
		cell_structure.exchange_ghosts_comm.comm[cnt].shift[dir] = boundary[2*dir+lr]*box_l[dir];
		cell_structure.update_ghost_pos_comm.comm[cnt].shift[dir] = boundary[2*dir+lr]*box_l[dir];
	      }
	      cnt++;
	    }
#ifdef PARTIAL_PERIODIC
	  if( PERIODIC(dir) || (boundary[2*dir+(1-lr)] == 0) ) 
#endif
	    if((node_pos[dir]+(1-i))%2==0) {
	      cnt++;
	    }
	}
      }
    }
  }
}

/** Init cell interactions for cell system domain decomposition.
 * initializes the interacting neighbor cell list of a cell The
 * created list of interacting neighbor cells is used by the verlet
 * algorithm (see verlet.cpp) to build the verlet lists.
 */
void dd_init_cell_interactions()
{
  int m,n,o,p,q,r,ind1,ind2,c_cnt=0,n_cnt;
 
  /* initialize cell neighbor structures */
  dd.cell_inter = (IA_Neighbor_List *) realloc(dd.cell_inter,local_cells.n*sizeof(IA_Neighbor_List));
  for(m=0; m<local_cells.n; m++) { 
    dd.cell_inter[m].nList = NULL; 
    dd.cell_inter[m].n_neighbors=0; 
  }

  /* loop all local cells */
  DD_LOCAL_CELLS_LOOP(m,n,o) {
    dd.cell_inter[c_cnt].nList = (IA_Neighbor *) realloc(dd.cell_inter[c_cnt].nList, CELLS_MAX_NEIGHBORS*sizeof(IA_Neighbor));
    dd.cell_inter[c_cnt].n_neighbors = CELLS_MAX_NEIGHBORS;
 
    n_cnt=0;
    ind1 = get_linear_index(m,n,o,dd.ghost_cell_grid);
    /* loop all neighbor cells */
    for(p=o-1; p<=o+1; p++)	
      for(q=n-1; q<=n+1; q++)
	for(r=m-1; r<=m+1; r++) {   
	  ind2 = get_linear_index(r,q,p,dd.ghost_cell_grid);
	  if(ind2 >= ind1) {
	    dd.cell_inter[c_cnt].nList[n_cnt].cell_ind = ind2;
	    dd.cell_inter[c_cnt].nList[n_cnt].pList    = &cells[ind2];
	    init_pairList(&dd.cell_inter[c_cnt].nList[n_cnt].vList);
	    n_cnt++;
	  }
	}
    c_cnt++;
  }
}

/*************************************************/

/** Returns pointer to the cell which corresponds to the position if
    the position is in the nodes spatial domain otherwise a NULL
    pointer. */
Cell *dd_save_position_to_cell(double pos[3]) 
{
  int i,cpos[3];
  double lpos;

  for(i=0;i<3;i++) {
    lpos = pos[i] - my_left[i];

    cpos[i] = (int)(lpos*dd.inv_cell_size[i])+1;
    
    /* particles outside our box. Still take them if
       VERY close or nonperiodic boundary */
    if (cpos[i] < 1) {
      if (lpos > -ROUND_ERROR_PREC*box_l[i]
#ifdef PARTIAL_PERIODIC
	  || (!PERIODIC(i) && boundary[2*i])
#endif
	  )
	cpos[i] = 1;
      else
	return NULL;
    }
    else if (cpos[i] > dd.cell_grid[i]) {
      if (lpos < local_box_l[i] + ROUND_ERROR_PREC*box_l[i]
#ifdef PARTIAL_PERIODIC
	  || (!PERIODIC(i) && boundary[2*i+1])
#endif
	  )
	cpos[i] = dd.cell_grid[i];
      else
	return NULL;
    }
  }
  i = get_linear_index(cpos[0],cpos[1],cpos[2], dd.ghost_cell_grid); 
  return &(cells[i]);  
}

Cell *dd_position_to_cell(double pos[3])
{
  int i,cpos[3];
  double lpos;

  for(i=0;i<3;i++) {
    lpos = pos[i] - my_left[i];

    cpos[i] = (int)(lpos*dd.inv_cell_size[i])+1;

    if (cpos[i] < 1) {
      cpos[i] = 1;
#ifdef ADDITIONAL_CHECKS
      if (PERIODIC(i) && lpos < -ROUND_ERROR_PREC*box_l[i]) {
	char *errtext = runtime_error(128 + 3*ES_DOUBLE_SPACE);
	ERROR_SPRINTF(errtext, "{005 particle @ (%g, %g, %g) is outside of the allowed cell grid} ", pos[0], pos[1], pos[2]);
      }
#endif
    }
    else if (cpos[i] > dd.cell_grid[i]) {
      cpos[i] = dd.cell_grid[i];
#ifdef ADDITIONAL_CHECKS
      if (PERIODIC(i) && lpos > local_box_l[i] + ROUND_ERROR_PREC*box_l[i]) {
	char *errtext = runtime_error(128 + 3*ES_DOUBLE_SPACE);
	ERROR_SPRINTF(errtext, "{005 particle @ (%g, %g, %g) is outside of the allowed cell grid} ", pos[0], pos[1], pos[2]);
      }
#endif
    }
  }
  i = get_linear_index(cpos[0],cpos[1],cpos[2], dd.ghost_cell_grid);  
  return &cells[i];
}

/*************************************************/

/** Append the particles in pl to \ref local_cells and update \ref local_particles.  
    @return 0 if all particles in pl reside in the nodes domain otherwise 1.*/
int dd_append_particles(ParticleList *pl, int fold_dir)
{
  int p, dir, c, cpos[3], flag=0, fold_coord=fold_dir/2;

  CELL_TRACE(fprintf(stderr, "%d: dd_append_particles %d\n", this_node, pl->n));

  for(p=0; p<pl->n; p++) {
    if(boundary[fold_dir] != 0)
      fold_coordinate(pl->part[p].r.p, pl->part[p].l.i, fold_coord);
    
    for(dir=0;dir<3;dir++) {
      cpos[dir] = (int)((pl->part[p].r.p[dir]-my_left[dir])*dd.inv_cell_size[dir])+1;

      if (cpos[dir] < 1) { 
	cpos[dir] = 1;
#ifdef PARTIAL_PERIODIC 
	if( PERIODIC(dir) ) 
#endif
	  {
	    flag=1;
	    CELL_TRACE(if(fold_coord==2){fprintf(stderr, "%d: dd_append_particles: particle %d (%f,%f,%f) not inside node domain.\n", this_node,pl->part[p].p.identity,pl->part[p].r.p[0],pl->part[p].r.p[1],pl->part[p].r.p[2]);});
	  }
      }
      else if (cpos[dir] > dd.cell_grid[dir]) {
	cpos[dir] = dd.cell_grid[dir];
#ifdef PARTIAL_PERIODIC 
	if( PERIODIC(dir) ) 
#endif
	  {
	    flag=1;
	    CELL_TRACE(if(fold_coord==2){fprintf(stderr, "%d: dd_append_particles: particle %d (%f,%f,%f) not inside node domain.\n", this_node,pl->part[p].p.identity,pl->part[p].r.p[0],pl->part[p].r.p[1],pl->part[p].r.p[2]);});
	  }
      }
    }
    c = get_linear_index(cpos[0],cpos[1],cpos[2], dd.ghost_cell_grid);
    CELL_TRACE(fprintf(stderr,"%d: dd_append_particles: Appen Part id=%d to cell %d\n",this_node,pl->part[p].p.identity,c));
    append_indexed_particle(&cells[c],&pl->part[p]);
  }
  return flag;
}
 
/*@}*/

/************************************************************/
/* Public Functions */
/************************************************************/

void dd_on_geometry_change(int flags) {
  /* check that the CPU domains are still sufficiently large. */
  for (int i = 0; i < 3; i++)
    if (local_box_l[i] < max_range) {
      char *errtext = runtime_error(128 + ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtext,"{013 box_l in direction %d is too small} ", i);
    }

  /* A full resorting is necessary if the grid has changed. We simply
     don't have anything fast for this case. Probably also not
     necessary. */
  if (flags & CELL_FLAG_GRIDCHANGED) {
    CELL_TRACE(fprintf(stderr,"%d: dd_on_geometry_change full redo\n",
		       this_node));
    cells_re_init(CELL_STRUCTURE_CURRENT);
    return;
  }

  /* otherwise, re-set our geometrical dimensions which have changed
     (in addition to the general ones that \ref grid_changed_box_l
     takes care of) */
  for(int i=0; i<3; i++) {
    dd.cell_size[i]       = local_box_l[i]/(double)dd.cell_grid[i];
    dd.inv_cell_size[i]   = 1.0 / dd.cell_size[i];
  }
  double min_cell_size = dmin(dmin(dd.cell_size[0],dd.cell_size[1]),dd.cell_size[2]);
  max_skin = min_cell_size - max_cut;

  CELL_TRACE(fprintf(stderr, "%d: dd_on_geometry_change: max_range = %f, min_cell_size = %f, max_skin = %f\n", this_node, max_range, min_cell_size, max_skin));
  
  if (max_range > min_cell_size) {
    /* if new box length leads to too small cells, redo cell structure
       using smaller number of cells. */
    cells_re_init(CELL_STRUCTURE_DOMDEC);
    return;
  }

  /* If we are not in a hurry, check if we can maybe optimize the cell
     system by using smaller cells. */
  if (!(flags & CELL_FLAG_FAST)) {
    int i;
    for(i=0; i<3; i++) {
      int poss_size = (int)floor(local_box_l[i]/max_range);
      if (poss_size > dd.cell_grid[i])
	break;
    }
    if (i < 3) {
      /* new range/box length allow smaller cells, redo cell structure,
	 possibly using smaller number of cells. */
      cells_re_init(CELL_STRUCTURE_DOMDEC);
      return;
    }
  }

  dd_update_communicators_w_boxl();

  /* tell other algorithms that the box length might have changed. */
  on_boxl_change();
}

/************************************************************/
void dd_topology_init(CellPList *old)
{
  int c,p,np;
  int exchange_data, update_data;
  Particle *part;

  CELL_TRACE(fprintf(stderr, "%d: dd_topology_init: Number of recieved cells=%d\n", this_node, old->n));

  /** broadcast the flag for using verlet list */
  MPI_Bcast(&dd.use_vList, 1, MPI_INT, 0, comm_cart);
 
  cell_structure.type             = CELL_STRUCTURE_DOMDEC;
  cell_structure.position_to_node = map_position_node_array;
  cell_structure.position_to_cell = dd_position_to_cell;

  /* set up new domain decomposition cell structure */
  dd_create_cell_grid();
  /* mark cells */
  dd_mark_cells();
  /* create communicators */
  dd_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);

  exchange_data = (GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION | GHOSTTRANS_POSSHFTD);
  update_data   = (GHOSTTRANS_POSITION | GHOSTTRANS_POSSHFTD);
  dd_prepare_comm(&cell_structure.exchange_ghosts_comm,  exchange_data);
  dd_prepare_comm(&cell_structure.update_ghost_pos_comm, update_data);

  dd_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);
  /* collect forces has to be done in reverted order! */
  dd_revert_comm_order(&cell_structure.collect_ghost_force_comm);

  dd_assign_prefetches(&cell_structure.ghost_cells_comm);
  dd_assign_prefetches(&cell_structure.exchange_ghosts_comm);
  dd_assign_prefetches(&cell_structure.update_ghost_pos_comm);
  dd_assign_prefetches(&cell_structure.collect_ghost_force_comm);

#ifdef LB
  dd_prepare_comm(&cell_structure.ghost_lbcoupling_comm, GHOSTTRANS_COUPLING) ;
  dd_assign_prefetches(&cell_structure.ghost_lbcoupling_comm) ;
#endif

  /* initialize cell neighbor structures */
  dd_init_cell_interactions();

  /* copy particles */
  for (c = 0; c < old->n; c++) {
    part = old->cell[c]->part;
    np   = old->cell[c]->n;
    for (p = 0; p < np; p++) {
      Cell *nc = dd_save_position_to_cell(part[p].r.p);
      /* particle does not belong to this node. Just stow away
	 somewhere for the moment */
      if (nc == NULL)
	nc = local_cells.cell[0];
      append_unindexed_particle(nc, &part[p]);
    }
  }
  for(c=0; c<local_cells.n; c++) {
    update_local_particles(local_cells.cell[c]);
  }
  CELL_TRACE(fprintf(stderr,"%d: dd_topology_init: done\n",this_node));
}

/************************************************************/
void dd_topology_release()
{
  int i,j;
  CELL_TRACE(fprintf(stderr,"%d: dd_topology_release:\n",this_node));
  /* release cell interactions */
  for(i=0; i<local_cells.n; i++) {
    for(j=0; j<dd.cell_inter[i].n_neighbors; j++) 
      free_pairList(&dd.cell_inter[i].nList[j].vList);
    dd.cell_inter[i].nList = (IA_Neighbor *) realloc(dd.cell_inter[i].nList,0);
  }
  dd.cell_inter = (IA_Neighbor_List *) realloc(dd.cell_inter,0);
  /* free ghost cell pointer list */
  realloc_cellplist(&ghost_cells, ghost_cells.n = 0);
  /* free ghost communicators */
  free_comm(&cell_structure.ghost_cells_comm);
  free_comm(&cell_structure.exchange_ghosts_comm);
  free_comm(&cell_structure.update_ghost_pos_comm);
  free_comm(&cell_structure.collect_ghost_force_comm);
#ifdef LB
  free_comm(&cell_structure.ghost_lbcoupling_comm);
#endif
}

/************************************************************/
void  dd_exchange_and_sort_particles(int global_flag)
{
  int dir, c, p, i, finished=0;
  ParticleList *cell,*sort_cell, send_buf_l, send_buf_r, recv_buf_l, recv_buf_r;
  Particle *part;
  CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles(%d):\n",this_node,global_flag));

  init_particlelist(&send_buf_l);
  init_particlelist(&send_buf_r);
  init_particlelist(&recv_buf_l);
  init_particlelist(&recv_buf_r);
  while(finished == 0 ) {
    finished=1;
    /* direction loop: x, y, z */  
    for(dir=0; dir<3; dir++) { 
      if(node_grid[dir] > 1) {
	/* Communicate particles that have left the node domain */
	/* particle loop */
	for(c=0; c<local_cells.n; c++) {
	  cell = local_cells.cell[c];
	  for (p = 0; p < cell->n; p++) {
	    part = &cell->part[p];
	    /* Move particles to the left side */
	    // Without the factor 0.5 in front of ROUND_ERROR_PREC, particles sitting exactly on the boundary 
	    // may be accepted (i.e. not sent) here and rejected later on by dd_save_position_to_cell
	    if(part->r.p[dir] - my_left[dir] < -0.5*ROUND_ERROR_PREC*box_l[dir]) {
#ifdef PARTIAL_PERIODIC 
	      if( PERIODIC(dir) || (boundary[2*dir]==0) ) 
#endif
		{
		  CELL_TRACE(fprintf(stderr,"%d: dd_ex_and_sort_p: send part left %d\n",this_node,part->p.identity));
		  local_particles[part->p.identity] = NULL;
		  move_indexed_particle(&send_buf_l, cell, p);
		  if(p < cell->n) p--;
		}
	    }
	    /* Move particles to the right side */
	    // Factor 0.5 see above
	    else if(part->r.p[dir] - my_right[dir] >= 0.5*ROUND_ERROR_PREC*box_l[dir]) {
#ifdef PARTIAL_PERIODIC 
	      if( PERIODIC(dir) || (boundary[2*dir+1]==0) ) 
#endif
		{
		  CELL_TRACE(fprintf(stderr,"%d: dd_ex_and_sort_p: send part right %d\n",this_node,part->p.identity));
		  local_particles[part->p.identity] = NULL;
		  move_indexed_particle(&send_buf_r, cell, p);
		  if(p < cell->n) p--;
		}
	    }
	    /* Sort particles in cells of this node during last direction */
	    else if(dir==2) {
	      sort_cell = dd_save_position_to_cell(part->r.p);
	      if(sort_cell != cell) {
		if(sort_cell==NULL) {
		  CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles: Take another loop",this_node));
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: CP1 Particle %d (%f,%f,%f) not inside node domain.\n",
				     this_node,part->p.identity,part->r.p[0],part->r.p[1],part->r.p[2]));		 
		  finished=0;
		  sort_cell = local_cells.cell[0];
		  if(sort_cell != cell) {
		    move_indexed_particle(sort_cell, cell, p);
		    if(p < cell->n) p--;
		  }
		}
		else {
		  move_indexed_particle(sort_cell, cell, p);
		  if(p < cell->n) p--;
		}
	      }
	    }
	  }
	}

	/* Exchange particles */
	if(node_pos[dir]%2==0) {
	  send_particles(&send_buf_l, node_neighbors[2*dir]);
	  recv_particles(&recv_buf_r, node_neighbors[2*dir+1]);
	  send_particles(&send_buf_r, node_neighbors[2*dir+1]);
	  recv_particles(&recv_buf_l, node_neighbors[2*dir]);
	}
	else {
	  recv_particles(&recv_buf_r, node_neighbors[2*dir+1]);
	  send_particles(&send_buf_l, node_neighbors[2*dir]);
	  recv_particles(&recv_buf_l, node_neighbors[2*dir]);
	  send_particles(&send_buf_r, node_neighbors[2*dir+1]);
	}
	/* sort received particles to cells */
	if(dd_append_particles(&recv_buf_l, 2*dir  ) && dir == 2) finished = 0;
	if(dd_append_particles(&recv_buf_r, 2*dir+1) && dir == 2) finished = 0; 
	/* reset send/recv buffers */
	send_buf_l.n = 0;
	send_buf_r.n = 0;
	recv_buf_l.n = 0;
	recv_buf_r.n = 0;
      }
      else {
	/* Single node direction case (no communication) */
	/* Fold particles that have left the box */
	/* particle loop */
	for(c=0; c<local_cells.n; c++) {
	  cell = local_cells.cell[c];
	  for (p = 0; p < cell->n; p++) {
	    part = &cell->part[p];
#ifdef PARTIAL_PERIODIC 
	    if( PERIODIC(dir) ) 
#endif
	      {
		fold_coordinate(part->r.p, part->l.i, dir);
	      }
	    if (dir==2) {
	      sort_cell = dd_save_position_to_cell(part->r.p);
	      if(sort_cell != cell) {
		if(sort_cell==NULL) {
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: CP2 Particle %d (%f,%f,%f) not inside node domain.\n",
				     this_node,part->p.identity,part->r.p[0],part->r.p[1],part->r.p[2]));
		  finished=0;
		  sort_cell = local_cells.cell[0];
		  if(sort_cell != cell) {
		    move_indexed_particle(sort_cell, cell, p);
		    if(p < cell->n) p--;
		  }      
		}
		else {
		  CELL_TRACE(fprintf(stderr, "%d: dd_exchange_and_sort_particles: move particle id %d\n", this_node,part->p.identity));
		  move_indexed_particle(sort_cell, cell, p);
		  if(p < cell->n) p--;
		}
	      }
	    }
	  }
	}
      }
    }

    /* Communicate wether particle exchange is finished */
    if(global_flag == CELL_GLOBAL_EXCHANGE) {
      if(this_node==0) {
	int sum;
	MPI_Reduce(&finished, &sum, 1, MPI_INT, MPI_SUM, 0, comm_cart);
	if( sum < n_nodes ) finished=0; else finished=sum; 
      } else {
	MPI_Reduce(&finished, NULL, 1, MPI_INT, MPI_SUM, 0, comm_cart);
      }
      MPI_Bcast(&finished, 1, MPI_INT, 0, comm_cart);
    } else {
      if(finished == 0) {
	char *errtext = runtime_error(128);
	ERROR_SPRINTF(errtext,"{004 some particles moved more than min_local_box_l, reduce the time step} ");
	/* the bad guys are all in cell 0, but probably their interactions are of no importance anyways.
	   However, their positions have to be made valid again. */
	finished = 1;
	/* all out of range coordinates in the left overs cell are moved to (0,0,0) */
	cell = local_cells.cell[0];
	for (p = 0; p < cell->n; p++) {
	  part = &cell->part[p];
	  if(dir < 3 && (part->r.p[dir] < my_left[dir] || part->r.p[dir] > my_right[dir]))
	    for (i = 0; i < 3; i++)
	      part->r.p[i] = 0;
	}
      }
    }
    CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles: finished value: %d\n",this_node,finished));
  }

  realloc_particlelist(&send_buf_l, 0);
  realloc_particlelist(&send_buf_r, 0);
  realloc_particlelist(&recv_buf_l, 0);
  realloc_particlelist(&recv_buf_r, 0);

#ifdef ADDITIONAL_CHECKS
  check_particle_consistency();
#endif

  CELL_TRACE(fprintf(stderr,"%d: dd_exchange_and_sort_particles finished\n",this_node));
}

/*************************************************/

int calc_processor_min_num_cells()
{
  int i, min = 1;
  /* the minimal number of cells can be lower if there are at least two nodes serving a direction,
     since this also ensures that the cell size is at most half the box length. However, if there is
     only one processor for a direction, there have to be at least two cells for this direction. */
  for (i = 0; i < 3; i++) 
    if (node_grid[i] == 1) 
      min *= 2;
  return min;
}

void calc_link_cell()
{
  int c, np1, n, np2, i ,j, j_start;
  Cell *cell;
  IA_Neighbor *neighbor;
  Particle *p1, *p2;
  double dist2, vec21[3];

  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {

    cell = local_cells.cell[c];
    p1   = cell->part;
    np1  = cell->n;

    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = &dd.cell_inter[c].nList[n];
      p2  = neighbor->pList->part;
      np2 = neighbor->pList->n;
      /* Loop cell particles */
      for(i=0; i < np1; i++) {
	j_start = 0;
	/* Tasks within cell: bonded forces */
	if(n == 0) {
	  add_bonded_force(&p1[i]);
#ifdef CONSTRAINTS
	  add_constraints_forces(&p1[i]);
#endif
    add_external_potential_forces(&p1[i]);
	  if (rebuild_verletlist)
	    memcpy(p1[i].l.p_old, p1[i].r.p, 3*sizeof(double));

	  j_start = i+1;
	}
	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {
#ifdef EXCLUSIONS
          if(do_nonbonded(&p1[i], &p2[j]))
#endif
	    {
	      dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);
	      add_non_bonded_pair_force(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
	    }
	}
      }
    }
  }
  rebuild_verletlist = 0;
}

/************************************************************/

void calculate_link_cell_energies()
{
  int c, np1, np2, n, i, j, j_start;
  Cell *cell;
  IA_Neighbor *neighbor;
  Particle *p1, *p2;
  double dist2, vec21[3];

  CELL_TRACE(fprintf(stderr,"%d: calculate link-cell energies\n",this_node));

  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p1   = cell->part;
    np1  = cell->n;
    /* calculate bonded interactions (loop local particles) */
    for(i = 0; i < np1; i++)  {
      add_kinetic_energy(&p1[i]);
      add_bonded_energy(&p1[i]);
#ifdef CONSTRAINTS
      add_constraints_energy(&p1[i]);
#endif
      add_external_potential_energy(&p1[i]);

      if (rebuild_verletlist)
        memcpy(p1[i].l.p_old, p1[i].r.p, 3*sizeof(double));
    }

    CELL_TRACE(fprintf(stderr,"%d: cell %d with %d neighbors\n",this_node,c, dd.cell_inter[c].n_neighbors));
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = &dd.cell_inter[c].nList[n];
      p2  = neighbor->pList->part;
      np2 = neighbor->pList->n;
      /* Loop cell particles */
      for(i=0; i < np1; i++) {
	j_start = 0;
	if(n == 0) j_start = i+1;
	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {	
#ifdef EXCLUSIONS
          if(do_nonbonded(&p1[i], &p2[j]))
#endif
	    {
	      dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);
	      add_non_bonded_pair_energy(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
	    }
	}
      }
    }
  }
  rebuild_verletlist = 0;
}

/************************************************************/

void calculate_link_cell_virials(int v_comp)
{
  int c, np1, np2, n, i, j, j_start;
  Cell *cell;
  IA_Neighbor *neighbor;
  Particle *p1, *p2;
  double dist2, vec21[3];

  CELL_TRACE(fprintf(stderr,"%d: calculate link-cell energies\n",this_node));

  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p1   = cell->part;
    np1  = cell->n;
    /* calculate bonded interactions (loop local particles) */
    for(i = 0; i < np1; i++)  {
      add_kinetic_virials(&p1[i], v_comp);
      add_bonded_virials(&p1[i]);
#ifdef BOND_ANGLE_OLD
      add_three_body_bonded_stress(&p1[i]);
#endif
#ifdef BOND_ANGLE
      add_three_body_bonded_stress(&p1[i]);
#endif
      if (rebuild_verletlist)
        memcpy(p1[i].l.p_old, p1[i].r.p, 3*sizeof(double));
    }

    CELL_TRACE(fprintf(stderr,"%d: cell %d with %d neighbors\n",this_node,c, dd.cell_inter[c].n_neighbors));
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = &dd.cell_inter[c].nList[n];
      p2  = neighbor->pList->part;
      np2 = neighbor->pList->n;
      /* Loop cell particles */
      for(i=0; i < np1; i++) {
	j_start = 0;
	if(n == 0) j_start = i+1;
	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {	
#ifdef EXCLUSIONS
          if(do_nonbonded(&p1[i], &p2[j]))
#endif
	    {
	      dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);
	      add_non_bonded_pair_virials(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
	    }
	}
      }
    }
  }
  rebuild_verletlist = 0;
}

/************************************************************/

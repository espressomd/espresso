#include "integrate.h"
#include "communication.h"
#include "domain_decomposition.h"
#include "debug.h"

/************************************************/
/** \name Defines */
/************************************************/
/*@{*/

/** half the number of cell neighbors in 3 Dimensions*/
#define CELLS_MAX_NEIGHBORS 14

/*@}*/

/************************************************/
/** \name Variables */
/************************************************/
/*@{*/

DomainDecomposition dd;

int max_num_cells = CELLS_MAX_NUM_CELLS;
double max_skin   = 0.0;

/** cell size. 
    Def: \verbatim cell_grid[i] = (int)(local_box_l[i]/max_range); \endverbatim */
double cell_size[3];
/** inverse cell size = \see cell_size ^ -1. */
double inv_cell_size[3];


/*@}*/


/************************************************************/
/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Calculate cell grid dimensions, cell sizes and number of cells.
 *  Calculates the cell grid, based on \ref local_box_l and \ref
 *  max_range. If the number of cells is larger than \ref
 *  max_num_cells, it increases max_range until the number of cells is
 *  smaller or equal \ref max_num_cells. It sets: \ref cell_grid, \ref
 *  ghost_cell_grid, \ref cell_size, \ref inv_cell_size, and \ref
 *  n_cells.
 */
void dd_create_cell_grid();

/** create communicators */
void  dd_prepare_comm(GhostCommunicator *comm, int data_parts);

/** initializes the interacting neighbor cell list of a cell
 *  (\ref #cells::neighbors).  the created list of interacting neighbor
 *  cells is used by the verlet algorithm (see verlet.c) to build the
 *  verlet lists.
 *
 *  @param i linear index of the cell.  */
void init_cell_neighbors(int i);
 
/** Calculate grid dimensions for the largest 3d grid in box where the
    cell size is given by size. 
    Help routine for \ref dd_ccreate_cell_grid.
    @param grid  resulting 3d grid.
    @param box   box where the cell grid should fit in.
    @param size  desired cell size.
    @param max   maximal number of cells.
*/
int calc_3d_cell_grid(int grid[3], double box[3], double size, int max);
/*@}*/

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

/** Convenient replace for inner cell check. usage: if(IS_INNER_CELL(m,n,o)) {...} */
#define DD_IS_LOCAL_CELL(m,n,o) \
  ( m > 0 && m < dd.ghost_cell_grid[0] - 1 && \
    n > 0 && n < dd.ghost_cell_grid[1] - 1 && \
    o > 0 && o < dd.ghost_cell_grid[2] - 1 ) 

/** Convenient replace for ghost cell check. usage: if(IS_GHOST_CELL(m,n,o)) {...} */
#define IS_GHOST_CELL(m,n,o) \
  ( m == 0 || m == dd.ghost_cell_grid[0] - 1 || \
    n == 0 || n == dd.ghost_cell_grid[1] - 1 || \
    o == 0 || o == dd.ghost_cell_grid[2] - 1 ) 

/************************************************************/
void dd_mark_cells()
{
  int m,n,o,cnt_c=0,cnt_l=0,cnt_g=0;
  DD_CELLS_LOOP(m,n,o)
    if(DD_IS_LOCAL_CELL(m,n,o)) local_cells.cell[cnt_l++] = &cells[cnt_c++]; 
    else                        ghost_cells.cell[cnt_g++] = &cells[cnt_c++];   
}

/************************************************************/
void dd_revert_comm_order(GhostCommunicator *comm)
{
  int i;
  GhostCommunication tmp;
  /* revret order */
  for(i=0; i<(comm->num/2); i++) {
    tmp = comm->comm[i];
    comm->comm[i] = comm->comm[comm->num-i-1];
    comm->comm[comm->num-i-1] = tmp;
  }
  /* exchange SEND/RECV */
  for(i=0; i<comm->num; i++) {
    if(comm->comm[i].type == GHOST_SEND) comm->comm[i].type = GHOST_RECV;
    else comm->comm[i].type = GHOST_SEND;
  }
}

/************************************************************/
void dd_topology_init(CellPList *cl)
{
  int c,p,np;
  Particle *part;

  CELL_TRACE(fprintf(stderr, "%d: dd_topology_init, Number of recieved cells=%d\n", this_node, cl->n));
  cell_structure.type             = CELL_STRUCTURE_DOMDEC;
  cell_structure.position_to_node = map_position_node_array;
  cell_structure.position_to_cell = dd_position_to_cell;

  /* not yet fully initialized */
  if(max_range <= 0) {
    max_range  = min_local_box_l/2.0; 
    max_range2 = SQR(max_range); 
  }

  /* set up new domain decomposition cell structure */
  dd_create_cell_grid();
  /* mark cells */
  dd_mark_cells();
  /* create communicators */
  dd_prepare_comm(&cell_structure.ghost_cells_comm,         GHOSTTRANS_PARTNUM);
  dd_prepare_comm(&cell_structure.exchange_ghosts_comm,     GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION);
  dd_prepare_comm(&cell_structure.update_ghost_pos_comm,    GHOSTTRANS_POSITION);
  dd_prepare_comm(&cell_structure.collect_ghost_force_comm, GHOSTTRANS_FORCE);
  /* collect forces has to be done in reverted order! */
  dd_revert_comm_order(&cell_structure.collect_ghost_force_comm);

  /* copy particles */
  fprintf(stderr,"%d: copy particles from %d cells!\n",this_node,cl->n);
  for (c = 0; c < cl->n; c++) {
    part = cl->cell[c]->part;
    np   = cl->cell[c]->n;
    for (p = 0; p < np; p++) {
      fprintf(stderr,"%d: copy part %d of old cell %d to new cs\n",this_node,p,c);
      // fprintf(stderr,"%d: copy part %d  to new cs\n",this_node,part[p].p.identity);
      append_unindexed_particle(cell_structure.position_to_cell(part[p].r.p), &part[p]);
    }
  }
  fprintf(stderr,"%d: update_local_particles for %d local cells\n",this_node,local_cells.n);
  for(c=0; c<n_cells; c++) fprintf(stderr,"%p ",&cells[c]);
  for(c=0; c<local_cells.n; c++) {
    fprintf(stderr,"%d: update_local_particles: cell %d with %p particles\n",this_node,c,local_cells.cell[c]);
    update_local_particles(local_cells.cell[c]);
  }
}

/************************************************************/
void dd_topology_release()
{
}

int dd_fill_pl_pointers(ParticleList **part_lists, int lc[3], int hc[3])
{
  int i,m,n,o,c=0;
  /* sanity check */
  for(i=0; i<3; i++) {
    if(lc[i]<0 || lc[i] >= dd.ghost_cell_grid[i]) return 0;
    if(hc[i]<0 || hc[i] >= dd.ghost_cell_grid[i]) return 0;
    if(lc[i] > hc[i]) return 0;
  }

  for(o=lc[2]; o<=hc[2]; o++) 
    for(n=lc[1]; n<=hc[1]; n++) 
      for(m=lc[0]; m<=hc[0]; m++) {
	i = get_linear_index(o,n,m,dd.ghost_cell_grid);
	part_lists[c] = &cells[i];
	c++;
      }
  return c;
}

/************************************************************/
void  dd_prepare_comm(GhostCommunicator *comm, int data_parts)
{
  int dir,lr,i,cnt, num=12, n_comm_cells[3];
  int lc[3],hc[3],done[3]={0,0,0};

#ifdef PARTIAL_PERIODIC
  /* calculate number of communications */
  num = 0;
  for(dir=0; dir<3; dir++) 
    for(lr=0; lr<2; lr++) 
    if( (periodic[dir] == 1) || (boundary[2*dir+lr] == 0) ) num += 2;
#endif

  /* prepare communicator */
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
      /* i: send/recv loop */
      for(i=0; i<2; i++) {  
#ifdef PARTIAL_PERIODIC
	if( (periodic[dir] == 1) || (boundary[2*dir+lr] == 0) ) 
#endif
	  if((node_pos[dir]+i)%2==0) {
	    comm->comm[cnt].type          = GHOST_SEND;
	    comm->comm[cnt].node          = node_neighbors[2*dir+lr];
	    comm->comm[cnt].part_lists    = malloc(n_comm_cells[dir]*sizeof(ParticleList *));
	    comm->comm[cnt].n_part_lists  = n_comm_cells[dir];

	    lc[(dir+0)%3] = hc[(dir+0)%3] = 1+lr*(dd.cell_grid[(dir+0)%3]-1);  
	    dd_fill_pl_pointers(comm->comm[cnt].part_lists,lc,hc);

	    //fprintf(stderr,"%d: comm %d send to   node %d grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]);
	    cnt++;
	  }
#ifdef PARTIAL_PERIODIC
	if( (periodic[dir] == 1) || (boundary[2*dir+(1-lr)] == 0) ) 
#endif
	  if((node_pos[dir]+(1-i))%2==0) {
	    comm->comm[cnt].type          = GHOST_RECV;
	    comm->comm[cnt].node          = node_neighbors[2*dir+(1-lr)];
	    comm->comm[cnt].part_lists    = malloc(n_comm_cells[dir]*sizeof(ParticleList *));
	    comm->comm[cnt].n_part_lists  = n_comm_cells[dir];

	    lc[(dir+0)%3] = hc[(dir+0)%3] = 0+(1-lr)*(dd.cell_grid[(dir+0)%3]+1);
	    dd_fill_pl_pointers(comm->comm[cnt].part_lists,lc,hc);

	    //fprintf(stderr,"%d: comm %d recv from node %d grid (%d,%d,%d)-(%d,%d,%d)\n",this_node,cnt,comm->comm[cnt].node,lc[0],lc[1],lc[2],hc[0],hc[1],hc[2]);
	    cnt++;
	  }
      }
    }
    done[dir]=1;
  }

}


/************************************************************/
Cell *dd_position_to_cell(double pos[3]) 
{
  int i,cpos[3];
  for(i=0;i<3;i++) {
    cpos[i] = (int)((pos[i]-my_left[i])*dd.inv_cell_size[i])+1;

#ifdef PARTIAL_PERIODIC
    if(periodic[i] == 0) {
      if (cpos[i] < 1)                 cpos[i] = 1;
      else if (cpos[i] > dd.cell_grid[i]) cpos[i] = dd.cell_grid[i];
    }
#endif

#ifdef ADDITIONAL_CHECKS
    if(cpos[i] < 1 || cpos[i] >  dd.cell_grid[i]) {
      fprintf(stderr,"%d: illegal cell position cpos[%d]=%d, ghost_grid[%d]=%d for pos[%d]=%f\n",this_node,i,cpos[i],i,dd.ghost_cell_grid[i],i,pos[i]);
      errexit();
    }
#endif

  }
  i = get_linear_index(cpos[0],cpos[1],cpos[2], dd.ghost_cell_grid); 
  return &(cells[i]);  
}

/************************************************************/
void dd_create_cell_grid()
{
  int i,n_local_cells,new_cells;
  double cell_range, min_box_l;
  CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid, max_range=%f local_box_l=(%f,%f,%f)\n", this_node,max_range,local_box_l[0],local_box_l[1],local_box_l[2]));
  
  /* initialize */
  min_box_l = dmin(dmin(local_box_l[0],local_box_l[1]),local_box_l[2]);
  cell_range = max_range;
  n_local_cells = max_num_cells+1;
  /* try finding a suitable cell grid. */
  while(n_local_cells > max_num_cells && 2.0*cell_range <= min_box_l) {
    n_local_cells = calc_3d_cell_grid(dd.cell_grid, local_box_l, cell_range, max_num_cells);
    cell_range *= 1.05;
  }
  /* quit program if unsuccesful */
  if(n_local_cells > max_num_cells) {
      fprintf(stderr, "%d: dd_create_cell_grid: grid (%d,%d,%d), n_local_cells=%d\n",
	      this_node,dd.cell_grid[0],dd.cell_grid[1],dd.cell_grid[2],n_local_cells);
      fprintf(stderr, "    Error: no suitable cell grid found (max_num_cells = %d)\n",
	      max_num_cells);
      fflush(stderr);
      errexit();
  } 
  /* now set all dependent variables */
  new_cells=1;
  for(i=0;i<3;i++) {
    dd.ghost_cell_grid[i] = dd.cell_grid[i]+2;	
    new_cells              *= dd.ghost_cell_grid[i];
    dd.cell_size[i]       = local_box_l[i]/(double)dd.cell_grid[i];
    dd.inv_cell_size[i]   = 1.0 / dd.cell_size[i];
  }
  cell_range = dmin(dmin(dd.cell_size[0],dd.cell_size[1]),dd.cell_size[2]);
  max_skin = cell_range - max_cut;

  /* allocate cell array and cell pointer arrays */
  realloc_cells(new_cells);
  realloc_cellplist(&local_cells, local_cells.n = n_local_cells);
  realloc_cellplist(&ghost_cells, ghost_cells.n = new_cells-n_local_cells);

  CELL_TRACE(fprintf(stderr, "%d: dd_create_cell_grid, n_cells=%d, local_cells.n=%d, ghost_cells.n=%d, dd.ghost_cell_grid=(%d,%d,%d)\n", this_node, n_cells,local_cells.n,ghost_cells.n,dd.ghost_cell_grid[0],dd.ghost_cell_grid[1],dd.ghost_cell_grid[2]));
}

/************************************************************/
int calc_3d_cell_grid(int grid[3], double box[3], double size, int max)
{
  int i,n_cells=1;
  for(i=0;i<3;i++) {
    grid[i] = (int)(box[i]/size);
    /* catch large integer case */
    if( grid[i] > max) {
      grid[0] = grid[1] = grid[2] = -1;
      return (max+1);
    }
    n_cells *= grid[i];
  }
  return n_cells;
}

#if 0

/************************************************************/
void init_cell_neighbors(int i)
{
  int j,m,n,o,cnt=0;
  int p1[3],p2[3];

  if(is_inner_cell(i,ghost_cell_grid)) { 
    cells[i].nList = (IA_Neighbor *) realloc(cells[i].nList,CELLS_MAX_NEIGHBORS*sizeof(IA_Neighbor));    
    get_grid_pos(i,&p1[0],&p1[1],&p1[2], ghost_cell_grid);
    /* loop through all neighbors */
    for(m=-1;m<2;m++)
      for(n=-1;n<2;n++)
	for(o=-1;o<2;o++) {
	  p2[0] = p1[0]+o;   p2[1] = p1[1]+n;   p2[2] = p1[2]+m;
	  j = get_linear_index(p2[0],p2[1],p2[2], ghost_cell_grid);
	  /* take the upper half of all neighbors 
	     and add them to the neighbor list */
	  if(j >= i) {
	    //CELL_TRACE(fprintf(stderr,"%d: cell %d neighbor %d\n",this_node,i,j));
	    cells[i].nList[cnt].cell_ind = j;
	    cells[i].nList[cnt].pList = &(cells[j].pList);
	    init_pairList(&(cells[i].nList[cnt].vList));
	    cnt++;
	  }
	}
    cells[i].n_neighbors = cnt;
  }
  else { 
    cells[i].n_neighbors = 0;
  }   
}


/************************************************************/
void cells_pre_init()
{
  int i;
  CELL_TRACE(fprintf(stderr,"%d: cells_pre_init():\n",this_node));

  /* set cell grid variables to a (1,1,1) grid */
  for(i=0;i<3;i++) {
    cell_grid[i] = 1;
    ghost_cell_grid[i] = 3;
  }
  n_cells = 27;

  /* allocate space */
  cells = (Cell *)malloc(n_cells*sizeof(Cell));
  for(i=0; i<n_cells; i++) init_cell(&cells[i]);
}

/*************************************************/
int pos_to_cell_grid_ind(double pos[3])
{
  int i,cpos[3];
  
  for(i=0;i<3;i++) {
    cpos[i] = (int)((pos[i]-my_left[i])*inv_cell_size[i])+1;

#ifdef PARTIAL_PERIODIC
    if(periodic[i] == 0) {
      if (cpos[i] < 1)                 cpos[i] = 1;
      else if (cpos[i] > cell_grid[i]) cpos[i] = cell_grid[i];
    }
#endif

#ifdef ADDITIONAL_CHECKS
    if(cpos[i] < 1 || cpos[i] >  cell_grid[i]) {
      fprintf(stderr,"%d: illegal cell position cpos[%d]=%d, ghost_grid[%d]=%d for pos[%d]=%f\n",this_node,i,cpos[i],i,ghost_cell_grid[i],i,pos[i]);
      errexit();
    }
#endif

  }
  return get_linear_index(cpos[0],cpos[1],cpos[2], ghost_cell_grid);  
}

/*************************************************/
int pos_to_capped_cell_grid_ind(double pos[3])
{
  int i,cpos[3];
  
  for(i=0;i<3;i++) {
    cpos[i] = (int)((pos[i]-my_left[i])*inv_cell_size[i])+1;

    if (cpos[i] < 1)
      cpos[i] = 1;
    else if (cpos[i] > cell_grid[i])
      cpos[i] = cell_grid[i];
  }
  return get_linear_index(cpos[0],cpos[1],cpos[2], ghost_cell_grid);  
}

void cells_changed_topology()
{
  int i;
  if(max_range <= 0) {
    /* not yet fully initialized */
    max_range = min_local_box_l/2.0;
  }

  for(i=0;i<3;i++) {
    cell_size[i] =  local_box_l[i];
    inv_cell_size[i] = 1.0 / cell_size[i];
  }
  dd_create_cell_grid() etc.;
}
#endif

/*************************************************/

int max_num_cells_callback(Tcl_Interp *interp, void *_data)
{
  int data = *(int *)_data;
  if (data < 27) {
    Tcl_AppendResult(interp, "WARNING: max_num_cells has to be at least 27. Set max_num_cells = 27!", (char *) NULL);
    data = 27;
  }
  max_num_cells = data;
  mpi_bcast_parameter(FIELD_MAXNUMCELLS);
  return (TCL_OK);
}

int find_node(double pos[3])
{
  int i, im[3];
  double f_pos[3];

  for (i = 0; i < 3; i++)
    f_pos[i] = pos[i];
  
  fold_position(f_pos, im);

  for (i = 0; i < 3; i++) {
    im[i] = (int)floor(node_grid[i]*f_pos[i]*box_l_i[i]);
#ifdef PARTIAL_PERIODIC
    if (!periodic[i]) {
      if (im[i] < 0)
	im[i] = 0;
      else if (im[i] >= node_grid[i])
	im[i] = node_grid[i] - 1;
    }
#endif
  }
  return map_array_node(im);
}


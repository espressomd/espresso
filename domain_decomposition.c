#include "domain_decomposition.h"
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

int cell_grid[3];
int ghost_cell_grid[3];
int max_num_cells = CELLS_MAX_NUM_CELLS;

/** cell size. 
    Def: \verbatim cell_grid[i] = (int)(local_box_l[i]/max_range); \endverbatim */
double cell_size[3];
/** inverse cell size = \see cell_size ^ -1. */
double inv_cell_size[3];

double max_skin=0.0;

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
void calc_cell_grid();

/** initializes the interacting neighbor cell list of a cell
 *  (\ref #cells::neighbors).  the created list of interacting neighbor
 *  cells is used by the verlet algorithm (see verlet.c) to build the
 *  verlet lists.
 *
 *  @param i linear index of the cell.  */
void init_cell_neighbors(int i);
 
/*@}*/


/** Convenient replace for inner cell check. usage: if(IS_INNER_CELL(m,n,o)) {...} */
#define IS_INNER_CELL(m,n,o) \
  ( m > 0 && m < ghost_cell_grid[0] - 1 && \
    n > 0 && n < ghost_cell_grid[1] - 1 && \
    o > 0 && o < ghost_cell_grid[2] - 1 ) 

/** Convenient replace for ghost cell check. usage: if(IS_GHOST_CELL(m,n,o)) {...} */
#define IS_GHOST_CELL(m,n,o) \
  ( m == 0 || m == ghost_cell_grid[0] - 1 || \
    n == 0 || n == ghost_cell_grid[1] - 1 || \
    o == 0 || o == ghost_cell_grid[2] - 1 ) 


void calc_cell_grid()
{
  int i;
  double cell_range;

  /* normal case */
  n_cells=1;
  for(i=0;i<3;i++) {
    ghost_cell_grid[i] = (int)(local_box_l[i]/max_range) + 2;
    /* make sure that at least one inner cell exists in all directions.
       Helps with highly anisotropic systems. */
    if (ghost_cell_grid[i] <= 2)
      ghost_cell_grid[i] = 3;
    /* for very large systems we might get integer multiplication overflow.
       but normally max_num_cells^3 is still ok */
    if (ghost_cell_grid[i] >= max_num_cells) {
      n_cells = max_num_cells + 1;
      break;
    }
    n_cells *= ghost_cell_grid[i];
  }

  /* catch case, n_cells > max_num_cells */
  if(n_cells > max_num_cells) {
    int count;
    double step;
    double min_box_l;
    double max_box_l;

    min_box_l = dmin(dmin(local_box_l[0],local_box_l[1]),local_box_l[2]);
    max_box_l = dmax(dmax(local_box_l[0],local_box_l[1]),local_box_l[2]);
    step = ((max_box_l/2.0)-max_range)/100; /* Maximal 100 trials! */
    if(step<0.0) {
      fprintf(stderr,"%d: calc_cell_grid: Error: negative step! Ask your local Guru\n",this_node);
      errexit();
    }
    cell_range = max_range;

    count = 100;
    while(n_cells > max_num_cells && --count > 0) {
      cell_range += step;
      n_cells=1;
      for(i=0;i<3;i++) {
	ghost_cell_grid[i] = (int)(local_box_l[i]/cell_range) + 2;
	/* make sure that at least one inner cell exists in all directions.
	   Helps with highly anisotropic systems. */
	if (ghost_cell_grid[i] <= 2)
	  ghost_cell_grid[i] = 3;
	/* for very large systems we might get integer multiplication overflow.
	   but normally max_num_cells^3 is still ok */
	if (ghost_cell_grid[i] >= max_num_cells) {
	  n_cells = max_num_cells + 1;
	  break;
	}
	n_cells *= ghost_cell_grid[i];
      }
    }
    if (count == 0) {
      fprintf(stderr, "%d: calc_cell_grid: Error: no suitable cell grid found (max_num_cells was %d)\n",
	      this_node,max_num_cells);
      errexit();
    }
    /* Store information about possible larger skin. */
  }
  cell_range = dmin(dmin(cell_size[0],cell_size[1]),cell_size[2]);
  max_skin = cell_range - max_cut;

  /* now set all dependent variables */
  for(i=0;i<3;i++) {
    cell_grid[i] = ghost_cell_grid[i]-2;	
    cell_size[i]     = local_box_l[i]/(double)cell_grid[i];
    inv_cell_size[i] = 1.0 / cell_size[i];
  }
}

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
  calc_cell_grid() etc.;
}

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
  mpi_bcast_event(PARAMETER_CHANGED);
  mpi_bcast_event(TOPOLOGY_CHANGED);
  return (TCL_OK);
}

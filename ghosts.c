/**************************************************/
/*******************  GHOSTS.C  *******************/
/**************************************************/

#include "ghosts.h"
#include "debug.h"

/** Creates an linear index list of a sub grid.
    The sub grid is defined by its lower and upper corner:\\
    from (lc[0],lc[1],lc[2]) to (hc[0],hc[1],hc[2])\\
    The grid dimension is given with (gs[0],gs[1],gs[2])\\
    The linear index list of length <returnvalue> is stored 
    in list starting from position start. max should be the 
    total length of list to ensure that the indices fit into list.
 */  
int sub_grid_indices(int* list, int start, int max,
		     int lc[3], int hc[3], int gs[3]);


void ghost_init()
{
  int i,j;
  /* ghost cell grid, cell grid */
  int gcg[3],cg[3];
  /* cell list sizes for the 3 direction, total number of sed cellsm, max(anz[i]). */
  int anz[3],n_total_cells,max_anz=0;
  /* send/recv frames start indizes, end indizes */
  int lc[3],hc[3],done[3]={0,0,0};

  GHOST_TRACE(fprintf(stderr,"%d: ghost_init:\n",this_node));

  /* preparation of help variables */
  for(i=0;i<3;i++) {
    gcg[i] = ghost_cell_grid[i];
    cg[i]  = cell_grid[i];
  }
  anz[0] = cg[1] *cg[2];
  anz[1] = cg[2] *gcg[0];
  anz[2] = gcg[0]*gcg[1];
  n_total_cells = 2*(anz[0] + anz[1] + anz[2]);
  cell_start[0] = 0;
  for(i=1;i<6;i++) cell_start[i] = cell_start[i-1] + anz[(i-1)/2];

  GHOST_TRACE(fprintf(stderr,"%d: cells: (%d, %d, %d) ghostcells (%d, %d, %d)\n",
		      this_node,cg[0],cg[1],cg[2],gcg[0],gcg[1],gcg[2]));
  GHOST_TRACE(fprintf(stderr,"%d: anz(%d,%d,%d) tot %d c_start(%d,%d,%d,%d,%d,%d)\n",
		      this_node,anz[0],anz[1],anz[2],n_total_cells,cell_start[0],cell_start[1],
		      cell_start[2],cell_start[3],cell_start[4],cell_start[5]));

  /* create send/recv cell index lists */
  send_cells = (int *)malloc(n_total_cells*sizeof(int));
  recv_cells = (int *)malloc(n_total_cells*sizeof(int));
  /* direction loop (sorry, it looks nasty, and it is!!!). */
  for(i=0;i<3;i++) {
    lc[(i+1)%3] = 1-done[(i+1)%3]; hc[(i+1)%3] = cg[(i+1)%3]+done[(i+1)%3];
    lc[(i+2)%3] = 1-done[(i+2)%3]; hc[(i+2)%3] = cg[(i+2)%3]+done[(i+2)%3];
    /* send to :   right, up, back */
    lc[(i+0)%3] = cg[(i+0)%3];     hc[(i+0)%3] = cg[(i+0)%3];
    n_send_cells[2*i] = sub_grid_indices(send_cells, cell_start[2*i], 
					 n_total_cells, lc, hc, gcg);
    /* recv from : right, up, back */
    lc[(i+0)%3] = cg[(i+0)%3]+1;   hc[(i+0)%3] = cg[(i+0)%3]+1;
    n_recv_cells[2*i] = sub_grid_indices(recv_cells, cell_start[2*i], 
					 n_total_cells, lc, hc, gcg);
    /* send to :   left, done, for*/
    lc[(i+0)%3] = 1;               hc[(i+0)%3] = 1;
    n_send_cells[(2*i)+1] = sub_grid_indices(send_cells, cell_start[(2*i)+1], 
					     n_total_cells, lc, hc, gcg);
     /* recv from : left, done, for*/
    lc[(i+0)%3] = 0;               hc[(i+0)%3] = 0;
    n_recv_cells[(2*i)+1] = sub_grid_indices(recv_cells, cell_start[(2*i)+1], 
					     n_total_cells, lc, hc, gcg);
    done[i] = 1;
  }
  
  for(i=0;i<6;i++) {
    GHOST_TRACE(fprintf(stderr,"%d: dir %d n_send_c %d: {",
			this_node,i,n_send_cells[i]));
    for(j=0;j<n_send_cells[i];j++)
      GHOST_TRACE(fprintf(stderr,"%d ",send_cells[j+cell_start[i]]));
    GHOST_TRACE(fprintf(stderr,"}\n"));
    GHOST_TRACE(fprintf(stderr,"%d: dir %d n_recv_c %d: {",
			this_node,i,n_recv_cells[i]));
    for(j=0;j<n_recv_cells[i];j++)
      GHOST_TRACE(fprintf(stderr,"%d ",recv_cells[j+cell_start[i]]));
    GHOST_TRACE(fprintf(stderr,"}\n"));
  }
    
  /* allocation of ghost cell information arrays */
  for(i=0;i<3;i++) if(anz[i]>max_anz) max_anz = anz[i];   
  n_send_ghosts = (int *)malloc(max_anz*sizeof(int));
  n_recv_ghosts = (int *)malloc(max_anz*sizeof(int));
 

  /* particle, force/pos buffers */
  buf_size = PART_INCREMENT;
  /*
  part_send_buf  = (Particle *)malloc(buf_size*sizeof(Particle));
  part_recv_buf  = (Particle *)malloc(buf_size*sizeof(Particle));
  part_send_buf2 = (Particle *)malloc(buf_size*sizeof(Particle));
  part_recv_buf2 = (Particle *)malloc(buf_size*sizeof(Particle));
  send_buf       = (double *)malloc(3*buf_size*sizeof(double));
  recv_buf       = (double *)malloc(3*buf_size*sizeof(double));
  */
  GHOST_TRACE(fprintf(stderr,"allocation done (exit ghost_init)\n"));
}

void exchange_part()
{
  GHOST_TRACE(fprintf(stderr,"%d: exchange_part:\n",this_node));
}

void exchange_ghost()
{
  GHOST_TRACE(fprintf(stderr,"%d: exchange_ghost:\n",this_node));
}

void exchange_ghost_pos()
{
  GHOST_TRACE(fprintf(stderr,"%d: exchange_ghost_pos:\n",this_node));
}

void exchange_ghost_forces()
{
  GHOST_TRACE(fprintf(stderr,"%d: exchange_ghost_forces:\n",this_node));
}

void ghost_exit()
{
  GHOST_TRACE(fprintf(stderr,"%d: ghost_exit:\n",this_node));
  /*
  free(send_cells);
  free(recv_cells);
  free(n_send_ghosts);
  free(n_recv_ghosts);
  if(buf_size>0) {
    free(part_send_buf);
    free(part_recv_buf);
    free(send_buf);
    free(recv_buf);
  }
  */
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
  int i;
  int size=0;
  int p0,p1,p2;
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

  // GHOST_TRACE(fprintf(stderr,"lilosg: lc (%d,%d,%d) hc (%d,%d,%d) gs (%d,%d,%d) size %d\n",
  //		      lc[0],lc[1],lc[2],hc[0],hc[1],hc[2],gs[0],gs[1],gs[2],size));

  i=start;
  for(p0=lc[0];p0<=hc[0];p0++)
    for(p1=lc[1];p1<=hc[1];p1++)
      for(p2=lc[2];p2<=hc[2];p2++) {
	list[i] = get_linear_index(p0,p1,p2,gs[0],gs[1],gs[2]);
	//	GHOST_TRACE(fprintf(stderr,"%d ",list[i]));
	i++;
      }
  //   GHOST_TRACE(fprintf(stderr,"\n"));

  return size;
}

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

  /* Init exchange particles */
  max_send_le = max_recv_le = max_send_ri = max_recv_ri = PART_INCREMENT;
  part_send_le_buf = (Particle *)malloc(max_send_le*sizeof(Particle));
  part_recv_le_buf = (Particle *)malloc(max_recv_le*sizeof(Particle));
  part_send_ri_buf = (Particle *)malloc(max_send_ri*sizeof(Particle));
  part_recv_ri_buf = (Particle *)malloc(max_recv_ri*sizeof(Particle));

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
 

  /* init exchange forces/positions  */
  max_send_buf = max_recv_buf = PART_INCREMENT;
  send_buf       = (double *)malloc(3*max_send_buf*sizeof(double));
  recv_buf       = (double *)malloc(3*max_recv_buf*sizeof(double));
  
  GHOST_TRACE(fprintf(stderr,"allocation done (exit ghost_init)\n"));
}

void exchange_part()
{
  int d,i,n;
  int le_ind,ri_ind;
  int recv_le,recv_ri;
  int pe_pos[3];

  GHOST_TRACE(fprintf(stderr,"%d: exchange_part:\n",this_node));

  map_node_array(this_node,&pe_pos[0],&pe_pos[1],&pe_pos[2]);

  /* fold coordinates to primary simulation box */
  for(n=0;n<n_particles;n++) fold_particle(particles[n].p,particles[n].i);

  for(d=0;d<3;d++) {                           /* direction loop */  
    le_ind=0; ri_ind=0;
    for(n=0;n<n_particles;n++) {               /* particle loop */
      if(particles[n].p[d] <  my_left[d] ) {   /* to the left */
	if(le_ind > max_send_le) {       
	  max_send_le += PART_INCREMENT;
	  part_send_le_buf = (Particle *)realloc(part_send_le_buf,max_send_le*sizeof(Particle));
	}
	memcpy(&part_send_le_buf[le_ind],&particles[n],sizeof(Particle));
	le_ind++;
      }
      if(particles[n].p[d] >= my_right[d]) {    /* to the right */
	if(ri_ind > max_send_ri) {
	  max_send_ri += PART_INCREMENT;
	  part_send_ri_buf = (Particle *)realloc(part_send_ri_buf,max_send_ri*sizeof(Particle));
	}
	memcpy(&part_send_ri_buf[ri_ind],&particles[n],sizeof(Particle));
	ri_ind++;
      }
    }
    /* send/recieve particles */
    /* First all processors with even pe grid position send their stuff 
       and then all the odd ones.  */
    //    if(pe_pos[d]%2==0) {
    //  MPI_Send(&le_ind, 1, MPI_INT, le_neighbors[d],0 , MPI_COMM_WORLD);
    //  MPI_Send(part_send_le_buf, le_ind, MPI_PARTICLE, le_neighbors[d]
    //	       ,0 , MPI_COMM_WORLD);
    //}
    //else {
    //  MPI_Recv(&max_recv_le, 1, MPI_INT, ri_neighbors[d],0 , MPI_COMM_WORLD);
    //  part_recv_le_buf = (Particle *)realloc(part_recv_le_buf,max_recv_le*sizeof(Particle));
    //  MPI_Recv(part_recv_le_buf, max_recv_le, MPI_PARTICLE, ri_neighbors[d]
    //       ,0 , MPI_COMM_WORLD);
    //}
  }
  
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
  
  free(send_cells);
  free(recv_cells);
  free(n_send_ghosts);
  free(n_recv_ghosts);
  
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


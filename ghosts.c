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
/** moves particle (ind) to the send buffer.
    subroutine of exchange_part(). */
void move_to_p_buf(int ind);
/** send particles in direction s_dir.
    subroutine of exchange_part(). */
void send_particles(int s_dir);
/** appends recieved particles to local particle array.
    subroutine of exchange_part(). */
void append_particles(void);

void send_ghosts(int s_dir);
void pack_ghost(Ghost *g_array, int g_ind, Particle *p_array, int p_ind);
void unpack_ghost(Particle *p_array, int p_ind, Ghost *g_array, int g_ind);

void send_posforce(int s_dir);

/******************************** variables for  exchange Particles */

/** Buffer for particles to send. */
int       n_p_send_buf;
int       max_p_send_buf;
Particle *p_send_buf;
/** Buffer for particles to recieve. */
int       n_p_recv_buf;
int       max_p_recv_buf;
Particle *p_recv_buf;

/******************************** variables for  exchange Ghosts */

/** maximal number of cells to send. */
int max_send_cells;
/** number of cells to send in direction X. */
int n_send_cells[6];
/** list of cell indices to send. */
int *send_cells;
/** number of cells to receive from direction X. */
int n_recv_cells[6];
/** list of cell indices to receive. */
int *recv_cells;
/** start indices for cells to send/recv in/from direction X. */ 
int cell_start[6];

/** Number of ghosts in each send cell. */ 
int *n_send_ghosts;
/** Number of ghosts in each recv cell. */ 
int *n_recv_ghosts;

/** Buffer for Ghosts to send. */
int   n_g_send_buf;
int   max_g_send_buf;
Ghost *g_send_buf;
/** Buffer for Ghosts to recieve. */
int   n_g_recv_buf;
int   max_g_recv_buf;
Ghost *g_recv_buf;

/** number of ghosts to send in direction X */
int ghost_send_size[6];
/** number of ghosts to recv from direction X */
int ghost_recv_size[6];

/** Buffer for forces/coordinates to send. */
double *send_buf;
int max_send_buf;
/** Buffer for forces/coordinates to recieve. */
double *recv_buf;
int max_recv_buf;

/*************************************** end variables */

void ghost_init()
{
  int i;
  /* ghost cell grid, cell grid */
  int gcg[3],cg[3];
  /* cell list sizes for the 3 direction, total number of sed cellsm. */
  int anz[3],n_total_cells;
  /* send/recv frames start indizes, end indizes */
  int lc[3],hc[3],done[3]={0,0,0};

  GHOST_TRACE(fprintf(stderr,"%d: ghost_init:\n",this_node));

  /* Init PE neighbors */
  calc_neighbors(this_node);

  /* Init exchange particles */
  max_p_send_buf = max_p_recv_buf = PART_INCREMENT;
  p_send_buf = (Particle *)malloc(max_p_send_buf*sizeof(Particle));
  p_recv_buf = (Particle *)malloc(max_p_recv_buf*sizeof(Particle));

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

  /* create send/recv cell index lists */
  send_cells = (int *)malloc(n_total_cells*sizeof(int));
  recv_cells = (int *)malloc(n_total_cells*sizeof(int));
  /* direction loop (sorry, it looks nasty, and it is!!!). */
  for(i=0;i<3;i++) {
    lc[(i+1)%3] = 1-done[(i+1)%3]; hc[(i+1)%3] = cg[(i+1)%3]+done[(i+1)%3];
    lc[(i+2)%3] = 1-done[(i+2)%3]; hc[(i+2)%3] = cg[(i+2)%3]+done[(i+2)%3];
    /* send to :   left, down, for */
    lc[(i+0)%3] = 1;               hc[(i+0)%3] = 1;
    n_send_cells[(2*i)] = sub_grid_indices(send_cells, cell_start[(2*i)], 
					     n_total_cells, lc, hc, gcg);
     /* recv from : right, up, back */
    lc[(i+0)%3] = 0;               hc[(i+0)%3] = 0;
    n_recv_cells[(2*i)+1] = sub_grid_indices(recv_cells, cell_start[(2*i)+1], 
					     n_total_cells, lc, hc, gcg);
    /* send to :   right, up, back */
    lc[(i+0)%3] = cg[(i+0)%3];     hc[(i+0)%3] = cg[(i+0)%3];
    n_send_cells[(2*i)+1] = sub_grid_indices(send_cells, cell_start[(2*i)+1], 
					 n_total_cells, lc, hc, gcg);
    /* recv from : left, down, for */
    lc[(i+0)%3] = cg[(i+0)%3]+1;   hc[(i+0)%3] = cg[(i+0)%3]+1;
    n_recv_cells[(2*i)] = sub_grid_indices(recv_cells, cell_start[(2*i)], 
					 n_total_cells, lc, hc, gcg);
    done[i] = 1;
  }
  
  /*
  for(i=0;i<6;i++) {
    int j;
    GHOST_TRACE(fprintf(stderr,"%d: dir %d: Send %d cells to node %d. Cells: {",
			this_node,i,n_send_cells[i],neighbors[i]));
    for(j=0;j<n_send_cells[i];j++)
      GHOST_TRACE(fprintf(stderr,"%d ",send_cells[j+cell_start[i]]));
    GHOST_TRACE(fprintf(stderr,"}\n"));
    GHOST_TRACE(fprintf(stderr,"%d: dir %d: Recv %d cells fr node %d. Cells: {",
     			this_node,i,n_recv_cells[i],neighbors[i]));
    for(j=0;j<n_recv_cells[i];j++)
      GHOST_TRACE(fprintf(stderr,"%d ",recv_cells[j+cell_start[i]]));
    GHOST_TRACE(fprintf(stderr,"}\n"));
  }
  */

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
  /* test 
  if(this_node ==0) {
    fprintf(stderr,"0: particles[0].p[0] = %f\n",particles[0].p[0]);
    particles[0].p[0] -= 5.0;
    fprintf(stderr,"0: new particles[0].p[0] = %f (myleft %f)\n",particles[0].p[0], my_left[0]);
  }
  end test */

  /* fold coordinates to primary simulation box */
  for(n=0;n<n_particles;n++) fold_particle(particles[n].p,particles[n].i);

  for(d=0; d<3; d++) {                           /* direction loop */  
    for(lr=0; lr<2; lr++) {
      dir = 2*d + lr;
      n_p_send_buf = n_p_recv_buf = 0;
      if(lr==0) {                                /* left */
	for(n=0;n<n_particles;n++)
	  if(particles[n].p[d] <   my_left[d] )  move_to_p_buf(n);
      }
      else {                                     /* right */
 	for(n=0;n<n_particles;n++)                  
	  if(particles[n].p[d] >=  my_right[d])  move_to_p_buf(n);
      }
      send_particles(dir);
      append_particles();
    }
    MPI_Barrier(MPI_COMM_WORLD);
  }
}

void exchange_ghost()
{
  int i,dir,c,n,m,c_ind;
  int old_max_send, old_max_recv;

  GHOST_TRACE(fprintf(stderr,"%d: exchange_ghost:\n",this_node));
  GHOST_TRACE(fprintf(stderr,"%d: start with n_p %d and n_g %d\n",this_node,n_particles,n_ghosts));

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
	pack_ghost(g_send_buf,n_g_send_buf,particles,cells[c_ind].particles[n]);
	n_g_send_buf++;
      }
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
	cells[c_ind].max_particles = n_recv_ghosts[c];
	cells[c_ind].particles = (int *)realloc(cells[c_ind].particles,
				 cells[c_ind].max_particles*sizeof(int));
      }
      for(i=0; i<n_recv_ghosts[c];i++) {
	unpack_ghost(particles,m,g_recv_buf,n);
	local_index[particles[m].identity] = m;
	cells[c_ind].particles[i] = m;
	m++; n++; 
      }
      cells[c_ind].n_particles = n_recv_ghosts[c];
    }
    if(dir%2==1) MPI_Barrier(MPI_COMM_WORLD);
  }                                               /* END direction loop */

  /* debug output */
  fprintf(stderr,"%d: max_send %d, max_recv %d (n_s,n_r) {",this_node,max_send_buf,max_recv_buf);
  for(i=0;i<6;i++) 
    fprintf(stderr,"(%d, %d),",ghost_send_size[i],ghost_recv_size[i]);
  fprintf(stderr,"}\n");
  GHOST_TRACE(fprintf(stderr,"%d: end with n_p %d and n_g %d\n",this_node,n_particles,n_ghosts));
  

  /* resize pos/force buffer if necessary */
  if(max_send_buf > old_max_send)
    send_buf = (double *)realloc(send_buf,3*max_send_buf*sizeof(double));
  if(max_recv_buf > old_max_recv)    
    recv_buf = (double *)realloc(recv_buf,3*max_recv_buf*sizeof(double));
}

void update_ghost_pos()
{
  int dir,c,c_ind,n,g;
  GHOST_TRACE(fprintf(stderr,"%d: update_ghost_pos:\n",this_node));
  for(dir=0; dir<6; dir++) {          /* direction loop forward */
    /* loop send cells -  copy positions to buffer*/
    g=0;
    for(c=0; c<n_send_cells[dir]; c++) {
      c_ind = send_cells[cell_start[dir]+c];
      for(n=0;n<cells[c_ind].n_particles;n++) {
	memcpy(&(send_buf[g]), particles[cells[c_ind].particles[n]].p, 3*sizeof(double));
	g += 3;
      }
    }
    GHOST_TRACE(fprintf(stderr,"%d: Dir %d send %d positions\n",this_node,g/3));
    /* send buffer */
    send_posforce(dir);
    /* loop recv cells - copy positions from buffer */
    g=0;
    for(c=0; c<n_recv_cells[dir]; c++) {
      c_ind = recv_cells[cell_start[dir]+c];
      for(n=0;n<cells[c_ind].n_particles;n++) {
	memcpy(particles[cells[c_ind].particles[n]].p, &(send_buf[g]), 3*sizeof(double));
	g += 3;
      }
    }
    if(dir%2==1) MPI_Barrier(MPI_COMM_WORLD);
  }
}

void collect_ghost_forces()
{
  int dir,c,c_ind,n,g,i;
  GHOST_TRACE(fprintf(stderr,"%d: collect_ghost_forces:\n",this_node));
  for(dir=5; dir>=0; dir--) {         /* direction loop backward */
    /* loop recv cells - copy forces to buffer*/
    g=0;
    for(c=0; c<n_recv_cells[dir]; c++) {
      c_ind = recv_cells[cell_start[dir]+c];
      for(n=0;n<cells[c_ind].n_particles;n++) {
	memcpy(&(send_buf[g]), particles[cells[c_ind].particles[n]].f, 3*sizeof(double));
	g += 3;
      }
    }
    GHOST_TRACE(fprintf(stderr,"%d: Dir %d send %d forces\n",this_node,g/3));
    /* send forces */
    send_posforce(dir);
    /* loop send cells - add buffer forces to local forces */
    g=0;
    for(c=0; c<n_send_cells[dir]; c++) {
      c_ind = send_cells[cell_start[dir]+c];
      for(n=0;n<cells[c_ind].n_particles;n++) {  
	for(i=0; i<3; i++) particles[cells[c_ind].particles[n]].f[i] += send_buf[g+i];
	g += 3;
      }
    }
    if(dir%2==1) MPI_Barrier(MPI_COMM_WORLD);
  }
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

  i=start;
  for(p0=lc[0];p0<=hc[0];p0++)
    for(p1=lc[1];p1<=hc[1];p1++)
      for(p2=lc[2];p2<=hc[2];p2++) {
	list[i] = get_linear_index(p0,p1,p2,gs[0],gs[1],gs[2]);
	i++;
      }

  return size;
}

void move_to_p_buf(int ind)
{
  /* check buffer size */
  if(n_p_send_buf >= max_p_send_buf) {
    max_p_send_buf += PART_INCREMENT;
    p_send_buf = (Particle *)realloc(p_send_buf,max_p_send_buf*sizeof(Particle));
  }
  /* copy particle from local array to send buffer */
  memcpy(&p_send_buf[n_p_send_buf],&particles[ind],sizeof(Particle));
  n_p_send_buf++;
  /* delete particle = 
     1: update the local_index list 
     2: move last particle (n_particles-1) to the free position (ind) */
  local_index[particles[ind].identity] = -1;
  local_index[particles[n_particles-1].identity] = ind;
  memcpy(&particles[ind],&particles[n_particles-1],sizeof(Particle));
  n_particles--;
}

void send_particles(int s_dir)
{
  int evenodd;
  int r_dir;
  MPI_Status status;

  /* calc recv direction (r_dir) from send direction (s_dir) */
  if(s_dir%2 == 0) r_dir = s_dir+1;
  else             r_dir = s_dir-1;
  /* two step communication: first all even positions than all odd */
  for(evenodd=0; evenodd<2;evenodd++) {
    if((pe_pos[s_dir/2]+evenodd)%2==0) {
      MPI_Send(&n_p_send_buf,1,MPI_INT,neighbors[s_dir],0,MPI_COMM_WORLD);
      MPI_Send(p_send_buf,n_p_send_buf*sizeof(Particle), MPI_BYTE, 
	       neighbors[s_dir],0,MPI_COMM_WORLD);
    }
    else {
      MPI_Recv(&n_p_recv_buf,1,MPI_INT,neighbors[r_dir],0,MPI_COMM_WORLD,&status);
      if(n_p_recv_buf >= max_p_recv_buf) {
	max_p_recv_buf = n_p_recv_buf;
	p_recv_buf = (Particle *)realloc(p_recv_buf,max_p_recv_buf*sizeof(Particle));
      }
      MPI_Recv(p_recv_buf,n_p_recv_buf*sizeof(Particle), MPI_BYTE,
	       neighbors[r_dir],0,MPI_COMM_WORLD,&status);
    }
  }
}

void append_particles(void)
{
  int n;
  /* check particle array size */
  if(n_particles + n_p_recv_buf >= max_particles) 
    realloc_particles(n_particles + n_p_recv_buf);
  /* copy particle data */
  memcpy(&particles[n_particles],p_send_buf,n_p_recv_buf*sizeof(Particle));
  /* update local_indes list */
  for(n=n_particles; n < n_particles+n_p_recv_buf; n++)
    local_index[particles[n].identity] = n;
  /* update particle number */
  n_particles += n_p_recv_buf;
}

void send_ghosts(int s_dir)
{
  int evenodd;
  int r_dir;
  MPI_Status status;

  /* calc recv direction (r_dir) from send direction (s_dir) */
  if(s_dir%2 == 0) r_dir = s_dir+1;
  else             r_dir = s_dir-1;
  /* two step communication: first all even positions than all odd */
  for(evenodd=0; evenodd<2;evenodd++) {
    if((pe_pos[s_dir/2]+evenodd)%2==0) {
      MPI_Send(n_send_ghosts, max_send_cells+1, MPI_INT,
	       neighbors[s_dir], 0, MPI_COMM_WORLD);
      MPI_Send(g_send_buf,n_send_ghosts[max_send_cells]*sizeof(Ghost), MPI_BYTE, 
	       neighbors[s_dir],0,MPI_COMM_WORLD);
    }
    else {
      MPI_Recv(n_recv_ghosts, max_send_cells+1, MPI_INT,
	       neighbors[r_dir],0,MPI_COMM_WORLD,&status);
      if(n_recv_ghosts[max_send_cells] >= max_g_recv_buf) {
	max_g_recv_buf = n_recv_ghosts[max_send_cells];
	g_recv_buf = (Ghost *)realloc(g_recv_buf,max_g_recv_buf*sizeof(Ghost));
      }
      MPI_Recv(g_recv_buf,n_recv_ghosts[max_send_cells]*sizeof(Ghost), MPI_BYTE,
	       neighbors[r_dir],0,MPI_COMM_WORLD,&status);
    }
  }
  /* store number of ghosts to send in direction s_dir */
  ghost_send_size[s_dir] = n_send_ghosts[max_send_cells];
  if(ghost_send_size[s_dir]>max_send_buf) max_send_buf = ghost_send_size[s_dir];
  /* store number of ghosts to recieve from direction r_dir */
  ghost_recv_size[r_dir] = n_recv_ghosts[max_send_cells];
  if(ghost_recv_size[r_dir]>max_recv_buf) max_recv_buf = ghost_recv_size[r_dir];
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

void send_posforce(int s_dir)
{
  int evenodd;
  int r_dir;
  MPI_Status status;

  /* calc recv direction (r_dir) from send direction (s_dir) */
  if(s_dir%2 == 0) r_dir = s_dir+1;
  else             r_dir = s_dir-1;
  /* two step communication: first all even positions than all odd */
  for(evenodd=0; evenodd<2;evenodd++) {
    if((pe_pos[s_dir/2]+evenodd)%2==0) {
      MPI_Send(send_buf,3*ghost_send_size[s_dir], MPI_DOUBLE, 
	       neighbors[s_dir],0,MPI_COMM_WORLD);
    }
    else {
      MPI_Recv(recv_buf,3*ghost_recv_size[r_dir], MPI_DOUBLE,
	       neighbors[r_dir],0,MPI_COMM_WORLD,&status);
    }
  }
}

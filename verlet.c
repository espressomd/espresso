/**************************************************/
/*******************  VERLET.C  *******************/
/**************************************************/

#include <stdio.h>
#include <stdlib.h>
#include "config.h"
#include "verlet.h"
#include "debug.h"
#include "cells.h"
#include "integrate.h"
#include "particle_data.h"
#include "communication.h"
#include "utils.h"

int   n_verletList;
int max_verletList;
int    *verletList;

int rebuild_verletlist = 1;

/** Granularity of the verlet list */
#define LIST_INCREMENT 100

/*******************  privat functions    *******************/
void resize_verlet_list();
double distance2(int p1, int p2);
void add_pair(int p1, int p2);

/*******************  exported functions  *******************/
void verlet_init()
{
  VERLET_TRACE(fprintf(stderr,"%d: verlet_init:\n",this_node));
  /* create empty verlet list of size LIST_INCREMENT */
  n_verletList   = 0;
  max_verletList = LIST_INCREMENT;
  verletList     = (int *)malloc(2*LIST_INCREMENT*sizeof(int));
}

void build_verlet_list()
{
  int i,j,nc;
  /* cell position */
  int m,n,o;
  /* cell indizes */
  int ci1,ci2;
  /* particle indizes */
  int p1,p2;
  /* pair distance square */
  double dist2;
  
  VERLET_TRACE(fprintf(stderr,"%d: build_verlet_list:\n",this_node));
  VERLET_TRACE(fprintf(stderr,"%d: Cell Grid: (%d, %d, %d)\n",
		       this_node,cell_grid[0],cell_grid[1],cell_grid[2]));
  VERLET_TRACE(fprintf(stderr,"%d: Max range %f:\n",this_node,sqrt(max_range2)));
 
  n_verletList   = 0; 

  /* define 'new old' coordinates*/
  for(i=0;i<n_particles+n_ghosts;i++)
    for(j=0;j<3;j++)
	particles[i].p_old[j] = particles[i].p[j];

  /* loop through all inner cells. */
  for(m=1; m<cell_grid[0]+1; m++)
    for(n=1; n<cell_grid[1]+1; n++)
      for(o=1; o<cell_grid[2]+1; o++) {
	ci1 = get_linear_index(m,n,o,ghost_cell_grid);
	/* interactions inside inner cells */
	for(i=0; i < cells[ci1].n_particles; i++) {
	  p1 = cells[ci1].particles[i];
	  for(j=i+1; j < cells[ci1].n_particles; j++) {
	    p2 = cells[ci1].particles[j];
	    dist2 = distance2(p1,p2);
	    VERLET_TRACE(fprintf(stderr,"%d: Pair (%d, %d) Dist %.2f",this_node,p1,p2,sqrt(dist2)));
	    if(dist2 <= max_range2) {
	      add_pair(p1,p2);
	      VERLET_TRACE(fprintf(stderr,"; added at %d",n_verletList));
	    }
	    VERLET_TRACE(fprintf(stderr,"\n"));
	  }
	}
	/* interactions with neighbour cells */
	for(nc=0; nc < cells[ci1].n_neighbours; nc++) {
	  ci2 = cells[ci1].neighbours[nc];
	  for(i=0; i < cells[ci1].n_particles; i++) {
	    p1 = cells[ci1].particles[i];
	    for(j=0; j < cells[ci2].n_particles; j++) {
	      p2 = cells[ci2].particles[j];
	      dist2 = distance2(p1,p2);
	      VERLET_TRACE(fprintf(stderr,"%d: Pair (%d, %d) Dist %.2f",this_node,p1,p2,sqrt(dist2)));
	      if(dist2 <= max_range2) {
		add_pair(p1,p2);
		VERLET_TRACE(fprintf(stderr,"; added at %d",n_verletList));
	      }
	      VERLET_TRACE(fprintf(stderr,"\n"));
	    }
	  }
	}
      }
  /* make the verlet list smaller again, if possible */
  resize_verlet_list();
  VERLET_TRACE(fprintf(stderr,"%d: n_verletList = %d \n",this_node,n_verletList));
  rebuild_verletlist = 0;
}

void verlet_exit()
{
  VERLET_TRACE(fprintf(stderr,"%d: verlet_exit:\n",this_node));
  free(verletList);
  n_verletList   = 0;
  max_verletList = 0; 
    
}

void resize_verlet_list()
{
  int diff;
  diff = max_verletList - n_verletList;
  if( diff > 2*LIST_INCREMENT ) {
    diff = (diff/LIST_INCREMENT)-1;
    max_verletList -= diff*LIST_INCREMENT;
    verletList = (int *)realloc(verletList, 2*max_verletList*sizeof(int));
  }
}

double distance2(int p1, int p2) 
{
  double dx,dy,dz;
 
  dx = particles[p1].p[0] - particles[p2].p[0];
  dy = particles[p1].p[1] - particles[p2].p[1];
  dz = particles[p1].p[2] - particles[p2].p[2];
  return (dx*dx + dy*dy + dz*dz);
}

void add_pair(int p1, int p2)
{
  /* check size of verlet List */
  if(n_verletList >= max_verletList) {
    max_verletList += LIST_INCREMENT;
    verletList = (int *)realloc(verletList, 2*max_verletList*sizeof(int));
  }
  /* add pair */
  verletList[2*n_verletList]   = p1;
  verletList[2*n_verletList+1] = p2;
  /* increase number of pairs */
  n_verletList++;
}

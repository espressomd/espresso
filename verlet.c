/** \file verlet.c   Verlet list.
 *  For more information see  \ref verlet.h "verlet.h"
 */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "verlet.h"
#include "debug.h"
#include "cells.h"
#include "integrate.h"
#include "particle_data.h"
#include "communication.h"
#include "utils.h"

/** Granularity of the verlet list */
#define LIST_INCREMENT 20

/*****************************************
 * Variables 
 *****************************************/

int rebuild_verletlist;

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Add a particle pair to a verlet pair list.
    Checks verlet pair list size and reallocates memory if necessary.
 *  \param p1 Pointer to paricle one.
 *  \param p2 Pointer to paricle two.
 *  \param vl Pointer to the verlet pair list.
 */
MDINLINE void add_pair(PairList *vl, Particle *p1, Particle *p2)
{
  /* check size of verlet List */
  if( (*vl).n+1 > (*vl).max ) {
    (*vl).max += LIST_INCREMENT;
    (*vl).pair = (Particle **)realloc((*vl).pair, 2*sizeof(Particle *)*(*vl).max);
  }
  /* add pair */
  (*vl).pair[(*vl).n * 2]       = p1;
  (*vl).pair[((*vl).n * 2) + 1] = p2;
  /* increase number of pairs */
  (*vl).n++;
}

/** Resizes a verlet pair list according to the actual content (*vl).n. 
    \param vl Pointer to the verlet pair list. */
void resize_verlet_list(PairList *vl);

/*@}*/

/*******************  exported functions  *******************/

void build_verlet_list()
{
  int c, nc, p1, p2;
  ParticleList *pl1, *pl2;
  PairList *vl;
  double dist2;

  VERLET_TRACE(fprintf(stderr,"%d: build_verlet_list:\n",this_node));

  /* loop cells */
  for(c=0; c<n_cells; c++) {
    pl1 = &(cells[c].pList);
    /* loop neighbor cells */
    for(nc=0; nc<cells[c].n_neighbors; nc++) {
      pl2 = cells[c].nList[nc].pList;
      vl  = &(cells[c].nList[nc].vList);
      /* loop cell particles */
      for(p1=0; p1<(*pl1).n; p1++) {
	/* set 'new old' coordinates */
	if(nc==0) memcpy((*pl1).part[p1].p_old, 
			 (*pl1).part[p1].r.p, 3*sizeof(double));
	/* loop neighbor cell particles */
	for(p2=0; p2<(*pl2).n; p2++) {
	  dist2 = distance2((*pl1).part[p1].r.p, (*pl2).part[p2].r.p);
	  if(dist2 <= max_range2) {
	    add_pair(vl, &((*pl1).part[p1]), &((*pl2).part[p2]));
	  }
	}
      }
      resize_verlet_list(vl);
    }
  }

#ifdef VERLET_DEBUG 
  {
    int sum,tot_sum=0;
    fprintf(stderr,"%d: Verlet list sizes: \n",this_node);
    for(c=0; c<n_cells; c++) { 
      sum=0;
      for(nc=0; nc<cells[c].n_neighbors; nc++) {
	sum += cells[c].nList[nc].vList.n;
      }
      fprintf(stderr,"%d: Cell %d: sum = %d \n",this_node,c,sum);
      tot_sum += sum;
    }
    fprintf(stderr,"%d: total number of interaction pairs: %d\n",this_node,tot_sum)
  }
#endif 

  rebuild_verletlist = 0;
}

/* preliminary routine ! */
void realloc_pairList(PairList *list, int size)
{
  int b_size;
  b_size = (LIST_INCREMENT*((size / LIST_INCREMENT)+1));
  if(b_size != (*list).max && b_size > (*list).n) {
    (*list).pair = (Particle **) realloc((*list).pair, 2*sizeof(Particle *)*b_size);
    (*list).max = b_size;
  }
}

/************************************************************/

void resize_verlet_list(PairList *vl)
{
  int diff;
  diff = (*vl).max - (*vl).n;
  if( diff > 2*LIST_INCREMENT ) {
    diff = (diff/LIST_INCREMENT)-1;
    (*vl).max -= diff*LIST_INCREMENT;
    (*vl).pair = (Particle **)realloc((*vl).pair, 2*sizeof(Particle *)*(*vl).max);
  }
}

/* Callback functions */
/************************************************************/

int rebuild_vlist_callback(Tcl_Interp *interp, void *_data)
{
  int data = *(int *)_data;
  if (data != 0 && data != 1) {
    Tcl_AppendResult(interp, "verletflag is an integrator flag and must be 0 or 1.", (char *) NULL);
    return (TCL_ERROR);
  }
  rebuild_verletlist = data;
  mpi_bcast_parameter(FIELD_VERLET);
  return (TCL_OK);
}




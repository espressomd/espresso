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
MDINLINE void add_pair(PairList *pl, Particle *p1, Particle *p2)
{
  /* check size of verlet List */
  if(pl->n >= pl->max) {
    pl->max += LIST_INCREMENT;
    pl->pair = (Particle **)realloc(pl->pair, 2*pl->max*sizeof(Particle *));
  }
  /* add pair */
  pl->pair[2*pl->n  ] = p1;
  pl->pair[2*pl->n+1] = p2;
  /* increase number of pairs */
  pl->n++;
}

/** Resizes a verlet pair list according to the actual content (*vl).n. 
    \param vl Pointer to the verlet pair list. */
void resize_verlet_list(PairList *pl);

/*@}*/

/*******************  exported functions  *******************/

void init_pairList(PairList *list)
{
  list->n       = 0;
  list->max     = 0;
  list->pair = NULL;
}


void build_verlet_lists()
{
  Cell *cell;
  PairList *pl;
  int i,j,nc;
  /* cell position */
  int m,n,o;
  /* particle lists */
  Particle *p1, *p2;
  int np1, np2;
  /* pair distance square */
  double dist2;

  VERLET_TRACE(fprintf(stderr,"%d: build_verlet_list:\n",this_node));

  INNER_CELLS_LOOP(m, n, o) {
    cell = CELL_PTR(m, n, o);
    p1   = cell->pList.part;
    np1  = cell->pList.n;

    for(nc=0; nc < cell->n_neighbors; nc++) {
      pl  = &cell->nList[nc].vList;
      p2  = cell->nList[nc].pList->part;
      np2 = cell->nList[nc].pList->n;
      for(i=0; i < np1; i++) {
	/* set 'new old' coordinates */
	if(nc == 0) memcpy(p1[i].p_old, 
			   p1[i].r.p, 3*sizeof(double));

	for(j = 0; j < np2; j++) {
	  dist2 = distance2(p1[i].r.p,p2[j].r.p);
	  if(dist2 <= max_range2) {
	    add_pair(pl, &p1[i], &p2[j]);
	  }
	}
      }
      /* make the verlet list smaller again, if possible */
      resize_verlet_list(pl);
    }
  }

#ifdef VERLET_DEBUG 
  {
    int c, sum,tot_sum=0;
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

/************************************************************/

void resize_verlet_list(PairList *pl)
{
  int diff;
  diff = pl->max - pl->n;
  if( diff > 2*LIST_INCREMENT ) {
    diff = (diff/LIST_INCREMENT)-1;
    pl->max -= diff*LIST_INCREMENT;
    pl->pair = (Particle **)realloc(pl->pair, 2*pl->max*sizeof(Particle *));
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

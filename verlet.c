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
#include "grid.h"

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
 *  \param pl Pointer to the verlet pair list.
 */
MDINLINE void add_pair(PairList *pl, Particle *p1, Particle *p2)
{
  /* check size of verlet List */
  if(pl->n+1 >= pl->max) {
    pl->max += LIST_INCREMENT;
    pl->pair = (Particle **)realloc(pl->pair, 2*pl->max*sizeof(Particle *));
  }
  /* add pair */
  pl->pair[(2*pl->n)  ] = p1;
  pl->pair[(2*pl->n)+1] = p2;
  /* increase number of pairs */
  pl->n++;
}

/** Resizes a verlet pair list according to the actual content (*vl).n. 
    \param pl Pointer to the verlet pair list. */
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
    
    /* interactions within the cell (neighbor cell 0)*/
    pl  = &cell->nList[0].vList;
    pl->n = 0;
    for(i=0; i < np1; i++) {
      memcpy(p1[i].p_old, p1[i].r.p, 3*sizeof(double));
      for(j = (i+1); j < np1; j++) {
	dist2 = distance2(p1[i].r.p,p1[j].r.p);
	if(dist2 <= max_range2) {
	  add_pair(pl, &p1[i], &p1[j]);
	  /* VERLET_TRACE(fprintf(stderr,"%d: cell(%d,%d,%d), nc=0, pair (%d-%d), dist=%f\n",
	     this_node,m,n,o,p1[i].r.identity, p1[j].r.identity,sqrt(dist2)));*/
	}
      }	
    }
    resize_verlet_list(pl);

    /* interactions with neighbor cells */
    for(nc=1; nc < cell->n_neighbors; nc++) {
      pl  = &cell->nList[nc].vList;
      pl->n = 0;
      p2  = cell->nList[nc].pList->part;
      np2 = cell->nList[nc].pList->n;
      for(i=0; i < np1; i++) {
	for(j = 0; j < np2; j++) {
	  dist2 = distance2(p1[i].r.p,p2[j].r.p);
	  if(dist2 <= max_range2) {
	    add_pair(pl, &p1[i], &p2[j]);
	    /* VERLET_TRACE(fprintf(stderr,"%d: cell(%d,%d,%d), nc=%d, pair (%d-%d), dist=%f\n",
	       this_node,m,n,o,nc,p1[i].r.identity, p2[j].r.identity,sqrt(dist2)));*/
	  }
	}	
      }
      resize_verlet_list(pl);
    }
  }

#ifdef VERLET_DEBUG 
  {
    int sum,tot_sum=0;
    int cind1,cind2;
    double estimate;

    estimate = 0.5*n_total_particles*(4.0/3.0*PI*pow(max_range,3.0))*(n_total_particles/(box_l[0]*box_l[1]*box_l[2]))/n_nodes;

    INNER_CELLS_LOOP(m, n, o) {
      cell = CELL_PTR(m, n, o);
      cind1 = get_linear_index(m,n,o,ghost_cell_grid);
      sum=0;
      for(nc=0; nc<cell->n_neighbors; nc++) {
	sum += cell->nList[nc].vList.n;
	cind2 = cell->nList[nc].cell_ind;
      }
      tot_sum += sum;
    }
    fprintf(stderr,"%d: total number of interaction pairs: %d (should be around %.1f)\n",this_node,tot_sum,estimate);
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

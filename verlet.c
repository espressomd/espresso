// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
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
#include "interaction_data.h"
#include "communication.h"
#include "utils.h"
#include "grid.h"
#include "forces.h"
#include "rotation.h"
#include "domain_decomposition.h"

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

void free_pairList(PairList *list)
{
  list->n       = 0;
  list->max     = 0;
  list->pair = (Particle **)realloc(list->pair, 0);
}

void build_verlet_lists()
{
  int c, np1, n, np2, i ,j, j_start;
  Cell *cell;
  IA_Neighbor neighbor;
  Particle *p1, *p2;
  PairList *pl;
  double dist2;

#ifdef VERLET_DEBUG 
  int estimate, sum=0;
  fprintf(stderr,"%d: build_verlet_list_and_force_calc:\n",this_node);
  /* estimate number of interactions: (0.5*n_part*ia_volume*density)/n_nodes */
  estimate = 0.5*n_total_particles*(4.0/3.0*PI*pow(max_range,3.0))*(n_total_particles/(box_l[0]*box_l[1]*box_l[2]))/n_nodes;
#endif
   
  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p1   = cell->part;
    np1  = cell->n;
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = dd.cell_inter[c].nList[n];
      p2  = neighbor.pList->part;
      np2 = neighbor.pList->n;
      /* init pair list */
      pl  = &neighbor.vList;
      pl->n = 0;
      /* Loop cell particles */
      for(i=0; i < np1; i++) {
	j_start = 0;
	/* Tasks within cell: store old position, avoid double counting */
	if(n == 0) {
	   memcpy(p1[i].l.p_old, p1[i].r.p, 3*sizeof(double));
	   j_start = i+1;
	}
	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {
	  dist2 = distance2(p1[i].r.p, p2[j].r.p);
	  if(dist2 <= max_range2) add_pair(pl, &p1[i], &p2[j]); 
	}
      }
      resize_verlet_list(pl);
      VERLET_TRACE(sum += pl->n);
    }
  }

    VERLET_TRACE(fprintf(stderr,"%d: total number of interaction pairs: %d (should be around %d)\n",this_node,sum,estimate));

  rebuild_verletlist = 0;
}

void calculate_verlet_ia()
{
  int c, np, n, i ,j;
  Cell *cell;
  Particle *p1, *p2, **pairs;
  double dist2, vec21[3];

  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p1   = cell->part;
    np  = cell->n;
    /* calculate bonded interactions (loop local particles) */
    for(i = 0; i < np; i++)  {
      add_bonded_force(&p1[i]);
#ifdef CONSTRAINTS
      add_constraints_forces(p1);
#endif
    }

    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      pairs = dd.cell_inter[c].nList[n].vList.pair;
      np    = dd.cell_inter[c].nList[n].vList.n;
      /* verlet list loop */
      for(i=0; i<2*np; i+=2) {
	p1 = pairs[i];                    /* pointer to particle 1 */
	p2 = pairs[i+1];                  /* pointer to particle 2 */
	dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);
	add_non_bonded_pair_force(p1, p2, vec21, sqrt(dist2), dist2);
      }
    }
  }
}

void build_verlet_lists_and_calc_verlet_ia()
{
  int c, np1, n, np2, i ,j, j_start;
  Cell *cell;
  IA_Neighbor neighbor;
  Particle *p1, *p2;
  PairList *pl;
  double dist2, vec21[3];
 
#ifdef VERLET_DEBUG 
  int estimate, sum=0;
  fprintf(stderr,"%d: build_verlet_list_and_calc_verlet_ia:\n",this_node);
  /* estimate number of interactions: (0.5*n_part*ia_volume*density)/n_nodes */
  estimate = 0.5*n_total_particles*(4.0/3.0*PI*pow(max_range,3.0))*(n_total_particles/(box_l[0]*box_l[1]*box_l[2]))/n_nodes;
#endif
 
  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p1   = cell->part;
    np1  = cell->n;
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = dd.cell_inter[c].nList[n];
      p2  = neighbor.pList->part;
      np2 = neighbor.pList->n;
      /* init pair list */
      pl  = &neighbor.vList;
      pl->n = 0;
      /* Loop cell particles */
      for(i=0; i < np1; i++) {
	j_start = 0;
	/* Tasks within cell: bonded forces, store old position, avoid double counting */
	if(n == 0) {
	  add_bonded_force(&p1[i]);
#ifdef CONSTRAINTS
	  add_constraints_forces(p1);
#endif
	  memcpy(p1[i].l.p_old, p1[i].r.p, 3*sizeof(double));
	  j_start = i+1;
	}
	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {
	  dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);
	  if(dist2 <= max_range2) {
	    ONEPART_TRACE(if(p1[i].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d,%d %d,%d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,i,n,j,sqrt(dist2)));
	    ONEPART_TRACE(if(p2[j].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d %d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,n,sqrt(dist2)));

	    add_pair(pl, &p1[i], &p2[j]); 
	    /* calc non bonded interactions */
	    add_non_bonded_pair_force(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
	  }
	}
      }
      resize_verlet_list(pl);
      VERLET_TRACE(sum += pl->n);
    }
  }

  VERLET_TRACE(fprintf(stderr,"%d: total number of interaction pairs: %d (should be around %d)\n",this_node,sum,estimate));
 
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
  mpi_bcast_parameter(FIELD_VERLETFLAG);
  return (TCL_OK);
}


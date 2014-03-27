/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
  This file is part of ESPResSo.
  
  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.
  
  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
  
  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>. 
*/
/** \file verlet.cpp   Verlet list.
 *  For more information see  \ref verlet.hpp "verlet.h"
 */
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "utils.hpp"
#include "verlet.hpp"
#include "cells.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "forces.hpp"
#include "energy.hpp"
#include "pressure.hpp"
#include "domain_decomposition.hpp"
#include "constraint.hpp"
#include "external_potential.hpp"

/** Granularity of the verlet list */
#define LIST_INCREMENT 20

/*****************************************
 * Variables 
 *****************************************/

/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Add a particle pair to a verlet pair list.
    Checks verlet pair list size and reallocates memory if necessary.
 *  \param p1 Pointer to particle one.
 *  \param p2 Pointer to particle two.
 *  \param pl Pointer to the verlet pair list.
 */
inline void add_pair(PairList *pl, Particle *p1, Particle *p2)
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
  IA_Neighbor *neighbor;
  Particle *p1, *p2;
  PairList *pl;
  double dist2;
#ifdef VERLET_DEBUG 
  double max_range_nonbonded2 = SQR(max_cut_nonbonded + skin);

  int estimate, sum=0;
  fprintf(stderr,"%d: build_verlet_list_and_force_calc:\n",this_node);
  /* estimate number of interactions: (0.5*n_part*ia_volume*density)/n_nodes */
  estimate = 0.5*n_part*(4.0/3.0*PI*pow(max_range_nonbonded2,1.5))*(n_part/(box_l[0]*box_l[1]*box_l[2]))/n_nodes;

  if (!dd.use_vList) { fprintf(stderr, "%d: build_verlet_lists, but use_vList == 0\n", this_node); errexit(); }
#endif
  
  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    VERLET_TRACE(fprintf(stderr,"%d: cell %d with %d neighbors\n",this_node,c, dd.cell_inter[c].n_neighbors));

    cell = local_cells.cell[c];
    p1   = cell->part;
    np1  = cell->n;
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = &dd.cell_inter[c].nList[n];
      p2  = neighbor->pList->part;
      np2 = neighbor->pList->n;
      /* init pair list */
      pl  = &neighbor->vList;
      pl->n = 0;

      /* no interaction set, Verlet list stays empty */
      if (max_cut_nonbonded == 0.0)
	continue;

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
#ifdef EXCLUSIONS
          if(do_nonbonded(&p1[i], &p2[j]))
#endif
	  {
	    dist2 = distance2(p1[i].r.p, p2[j].r.p);
	    if(dist2 <= SQR(get_ia_param(p1[i].p.type, p2[j].p.type)->max_cut + skin))
	      add_pair(pl, &p1[i], &p2[j]);
	  }
	}
      }
      resize_verlet_list(pl);
      VERLET_TRACE(fprintf(stderr,"%d: neighbor %d has %d particles\n",this_node,n,pl->n));
      VERLET_TRACE(sum += pl->n);
    }
  }

  rebuild_verletlist = 0;

  VERLET_TRACE(fprintf(stderr,"%d: total number of interaction pairs: %d (should be around %d)\n",this_node,sum,estimate));
}

void calculate_verlet_ia()
{
  int c, np, n, i;
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
      add_constraints_forces(&p1[i]);
#endif
      add_external_potential_forces(&p1[i]);
    }

    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      pairs = dd.cell_inter[c].nList[n].vList.pair;
      np    = dd.cell_inter[c].nList[n].vList.n;
      /* verlet list loop */
      for(i=0; i<2*np; i+=2) {
	p1 = pairs[i];                    /* pointer to particle 1 */
	p2 = pairs[i+1];                  /* pointer to particle 2 */
	dist2 = distance2vec(p1->r.p, p2->r.p, vec21);
	add_non_bonded_pair_force(p1, p2, vec21, sqrt(dist2), dist2);
      }
    }
  }
}

void build_verlet_lists_and_calc_verlet_ia()
{
  int c, np1, n, np2, i ,j, j_start;
  Cell *cell;
  IA_Neighbor *neighbor;
  Particle *p1, *p2;
  PairList *pl;
  double dist2, vec21[3];
 
#ifdef VERLET_DEBUG 
  int estimate, sum=0;
  double max_range_nonbonded2 = SQR(max_cut_nonbonded + skin);

  fprintf(stderr,"%d: build_verlet_list_and_calc_verlet_ia:\n",this_node);
  /* estimate number of interactions: (0.5*n_part*ia_volume*density)/n_nodes */
  estimate = 0.5*n_part*(4.0/3.0*PI*pow(max_range_nonbonded2,1.5))*(n_part/(box_l[0]*box_l[1]*box_l[2]))/n_nodes;

  if (!dd.use_vList) { fprintf(stderr, "%d: build_verlet_lists, but use_vList == 0\n", this_node); errexit(); }
#endif
 
  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    VERLET_TRACE(fprintf(stderr,"%d: cell %d with %d neighbors\n",this_node,c, dd.cell_inter[c].n_neighbors));

    cell = local_cells.cell[c];
    p1   = cell->part;
    np1  = cell->n;
    
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = &dd.cell_inter[c].nList[n];
      p2  = neighbor->pList->part;
      np2 = neighbor->pList->n;
      VERLET_TRACE(fprintf(stderr,"%d: neighbor %d contains %d parts\n",this_node,n,np2));
      /* init pair list */
      pl  = &neighbor->vList;
      pl->n = 0;
      /* Loop cell particles */
      for(i=0; i < np1; i++) {
	j_start = 0;
	/* Tasks within cell: bonded forces, store old position, avoid double counting */
	if(n == 0) {
	  add_bonded_force(&p1[i]);
#ifdef CONSTRAINTS
	  add_constraints_forces(&p1[i]);
#endif
    add_external_potential_forces(&p1[i]);
	  memcpy(p1[i].l.p_old, p1[i].r.p, 3*sizeof(double));
	  j_start = i+1;
	}
	
	/* no interaction set, no need for particle pairs */
	if (max_cut_nonbonded == 0.0)
	  continue;

	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {
#ifdef EXCLUSIONS
          if(do_nonbonded(&p1[i], &p2[j]))
#endif
	  {
	  dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);

	  VERLET_TRACE(fprintf(stderr,"%d: pair %d %d has distance %f\n",this_node,p1[i].p.identity,p2[j].p.identity,sqrt(dist2)));

	  if(dist2 <= SQR(get_ia_param(p1[i].p.type, p2[j].p.type)->max_cut + skin)) {
	    ONEPART_TRACE(if(p1[i].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d,%d %d,%d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,i,n,j,sqrt(dist2)));
	    ONEPART_TRACE(if(p2[j].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d %d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,n,sqrt(dist2)));

	    add_pair(pl, &p1[i], &p2[j]);
	    /* calc non bonded interactions */
	    add_non_bonded_pair_force(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
	  }
	 }
	}
      }
      resize_verlet_list(pl);
      VERLET_TRACE(fprintf(stderr,"%d: neighbor %d has %d pairs\n",this_node,n,pl->n));
      VERLET_TRACE(sum += pl->n);
    }
  }

  VERLET_TRACE(fprintf(stderr,"%d: total number of interaction pairs: %d (should be around %d)\n",this_node,sum,estimate));
 
  rebuild_verletlist = 0;
}

/************************************************************/

void calculate_verlet_energies()
{
  int c, np, n, i;
  Cell *cell;
  Particle *p1, *p2, **pairs;
  double dist2, vec21[3];

  VERLET_TRACE(fprintf(stderr,"%d: calculate verlet energies\n",this_node));

  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p1   = cell->part;
    np  = cell->n;
    /* calculate bonded interactions (loop local particles) */
    for(i = 0; i < np; i++)  {
      add_kinetic_energy(&p1[i]);
      add_bonded_energy(&p1[i]);
#ifdef CONSTRAINTS
      add_constraints_energy(&p1[i]);
#endif
      add_external_potential_energy(&p1[i]);
    }

    /* no interaction set */
    if (max_cut_nonbonded == 0.0)
      continue;

    VERLET_TRACE(fprintf(stderr,"%d: cell %d with %d neighbors\n",this_node,c, dd.cell_inter[c].n_neighbors));
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      pairs = dd.cell_inter[c].nList[n].vList.pair;
      np    = dd.cell_inter[c].nList[n].vList.n;
      VERLET_TRACE(fprintf(stderr,"%d: neighbor %d has %d particles\n",this_node,n,np));

      /* verlet list loop */
      for(i=0; i<2*np; i+=2) {
	p1 = pairs[i];                    /* pointer to particle 1 */
	p2 = pairs[i+1];                  /* pointer to particle 2 */
	dist2 = distance2vec(p1->r.p, p2->r.p, vec21);
	VERLET_TRACE(fprintf(stderr, "%d: %d <-> %d: dist2 dist2\n",this_node,p1->p.identity,p2->p.identity));
	add_non_bonded_pair_energy(p1, p2, vec21, sqrt(dist2), dist2);
      }
    }
  }
}

/************************************************************/

void calculate_verlet_virials(int v_comp)
{
  int c, np, n, i;
  Cell *cell;
  Particle *p1, *p2, **pairs;
  double dist2, vec21[3];

  VERLET_TRACE(fprintf(stderr,"%d: calculate verlet pressure\n",this_node));

  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p1   = cell->part;
    np  = cell->n;
    /* calculate bonded interactions (loop local particles) */
    for(i = 0; i < np; i++)  {
      add_kinetic_virials(&p1[i],v_comp);
      add_bonded_virials(&p1[i]);
#ifdef BOND_ANGLE_OLD
      add_three_body_bonded_stress(&p1[i]);
#endif
#ifdef BOND_ANGLE
      add_three_body_bonded_stress(&p1[i]);
#endif
    }

    /* no interaction set */
    if (max_cut_nonbonded == 0.0)
      continue;

    VERLET_TRACE(fprintf(stderr,"%d: cell %d with %d neighbors\n",this_node,c, dd.cell_inter[c].n_neighbors));
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      pairs = dd.cell_inter[c].nList[n].vList.pair;
      np    = dd.cell_inter[c].nList[n].vList.n;
      VERLET_TRACE(fprintf(stderr,"%d: neighbor %d has %d particles\n",this_node,n,np));

      /* verlet list loop */
      for(i=0; i<2*np; i+=2) {
	p1 = pairs[i];                    /* pointer to particle 1 */
	p2 = pairs[i+1];                  /* pointer to particle 2 */
	dist2 = distance2vec(p1->r.p, p2->r.p, vec21);
	add_non_bonded_pair_virials(p1, p2, vec21, sqrt(dist2), dist2);
      }
    }
  }
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


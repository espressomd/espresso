// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file iccp3m.c 

    <b>ICCP3M: A novel iterative scheme for  </b>
   This file contains function needed to compute induced charges on the dielectric
   interface for any geometry. Method utilizes idea of ICC (Boda et al PRE 69,046702, 2004) 
   and here we are using p3m to compute forces for determining induced charges by an iterative scheme
   which is proposed by Sandeep Tyagi, for matrix equation Ah=c (refer to original paper) rather then
   converting whole A matrix. Induced charge grids are point charges which is defined as first set of particles
   (id has to be specified in the options of iccp3m command)

   For more information about the iccp3m algorithm
   The corresponding header file is \ref iccp3m.h "iccp3m.h".

*/

#include <stdlib.h>
#include <stdio.h>
/* #define _GNU_SOURCE  */
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>
#include "iccp3m.h"
#include "tcl.h"
#include "parser.h"
#include "p3m.h"
#include "verlet.h"
#include "utils.h"
#include "cells.h"
#include "particle_data.h"
#include "domain_decomposition.h"
#include "verlet.h"
#include "forces.h"
#include "maggs.h"
#include "config.h"
#include "global.h"
#include "elc.h"

#ifdef ELP3M

iccp3m_struct iccp3m_cfg;

/* p3m functions that is used in iccp3m to avoid short range forces computation! */
void init_forces_iccp3m();
void calc_long_range_forces_iccp3m();
MDINLINE void add_pair_iccp3m(PairList *pl, Particle *p1, Particle *p2);
void resize_verlet_list_iccp3m(PairList *pl);
MDINLINE void init_local_particle_force_iccp3m(Particle *part);
MDINLINE void init_ghost_force_iccp3m(Particle *part);

Cell *local;
CellPList me_do_ghosts;

/* other iccp3m functions */
int parse_normal(Tcl_Interp *interp,int normal_args, char *string);
int parse_area(Tcl_Interp *interp,int normal_args, char *string);
int imod(int x,int y);
void iccp3m_revive_forces();
void iccp3m_store_forces();

/* Granularity of the verlet list */
#define LIST_INCREMENT 20

int iccp3m(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int last_ind_id,num_iteration,normal_args,normal_check=0,area_args,area_check,update;
  char buffer[TCL_DOUBLE_SPACE];
  double e1,e2,convergence,relax;
    iccp3m_cfg.selection=1;  /* now it is selected */
    /* printf("Warning: iccp3m command update option has no effect for now (reserved option for future )!\n"); */
   if(argc != 10) { 
  Tcl_AppendResult(interp, "Wrong # of args! Usage: iccp3m <last_ind_id> <e1> <e2> <num_iteration> <convergence> <relaxation> <area> <normal_components> <update>", (char *)NULL); 
      return (TCL_ERROR); 
     }
     if(!ARG_IS_I(1, last_ind_id)) {
       Tcl_ResetResult(interp);
       Tcl_AppendResult(interp, "Last induced id must be integer (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
       } else if(last_ind_id < 2) { 
       Tcl_ResetResult(interp);
       Tcl_AppendResult(interp, "Last induced id can not be smaller then 2 (got: ", argv[1],")!", (char *)NULL); return (TCL_ERROR);
      }
     if(!ARG_IS_D(2, e1)) {
       Tcl_ResetResult(interp);
       Tcl_AppendResult(interp, "Dielectric constant e1(inner) must be double(got: ", argv[2],")!", (char *)NULL); return (TCL_ERROR);
      }
     if(!ARG_IS_D(3, e2)) {
       Tcl_ResetResult(interp);
       Tcl_AppendResult(interp, "Dielectric constant e2(outer) must be double (got: ", argv[3],")!", (char *)NULL); return (TCL_ERROR);
      }
     if(!ARG_IS_I(4, num_iteration)) {
       Tcl_ResetResult(interp);
       Tcl_AppendResult(interp, "Number of maximum iterations must be integer (got: ", argv[4],")!", (char *)NULL); return (TCL_ERROR);
      }
     if(!ARG_IS_D(5, convergence)) {
       Tcl_ResetResult(interp);
       Tcl_AppendResult(interp, "Convergence criterion must be double (got: ", argv[5],")!", (char *)NULL); return (TCL_ERROR);
       } else if( (convergence < 1e-10) || (convergence > 1e-1) ) { 
         Tcl_ResetResult(interp);
         Tcl_AppendResult(interp, "Convergence criterion can not be smaller then 1e-10 and greater then 1e-2(got: ", argv[5],")!", (char *)NULL); return (TCL_ERROR);
      }
       if(!ARG_IS_D(6, relax)) {
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "Relaxation parameter must be double (got: ", argv[6],")!", (char *)NULL); return (TCL_ERROR);
        }

       if(!ARG_IS_I(9, update)) {
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "iccp3m Update must be integer (got: ", argv[9],")!", (char *)NULL); return (TCL_ERROR);
          } else if(update < 1) { 
           Tcl_ResetResult(interp);
           Tcl_AppendResult(interp, "iccp3m Update not be smaller then 1 (got: ", argv[9],")!", (char *)NULL); return (TCL_ERROR);
        }

  iccp3m_cfg.last_ind_id=last_ind_id;      /* Assign needed options */
  iccp3m_cfg.num_iteration=num_iteration; /* maximum number of iterations to go */
  iccp3m_cfg.e1=e1;
  iccp3m_cfg.e2=e2;
  iccp3m_cfg.convergence=convergence;
  iccp3m_cfg.relax=relax;
  iccp3m_cfg.update=update;

   normal_args=(iccp3m_cfg.last_ind_id+1)*3;
   /* Now get Normal Vectors components consecutively */
     normal_check=parse_normal(interp,normal_args,argv[8]);
         if(normal_check == 1) { 
           Tcl_ResetResult(interp);
           Tcl_AppendResult(interp, "Error in number of normal vectors", argv[8],")!", (char *)NULL); return (TCL_ERROR);
          }
     area_args=(iccp3m_cfg.last_ind_id+1);
     /* Now get area of the boundary elements */
     area_check=parse_area(interp,area_args,argv[7]);
          if(area_check == 1) {
           Tcl_ResetResult(interp);
           Tcl_AppendResult(interp, "Error in number of area vectors", argv[7],")!", (char *)NULL); return (TCL_ERROR);
          }
   Tcl_PrintDouble(interp,iccp3m_iteration(),buffer); 
   Tcl_AppendResult(interp, buffer, (char *) NULL);
  return TCL_OK;
}
int iccp3m_iteration() {
   double fdot,hold,hnew,del_eps,diff=0.0,difftemp=0.0,qold;
   Cell *cell;
   int c,i,np;
   Particle *part;
   int j;

 /* allocate force arrays for getting for interface particles */
  iccp3m_cfg.fx=malloc(sizeof(double)*n_total_particles);
  iccp3m_cfg.fy=malloc(sizeof(double)*n_total_particles);
  iccp3m_cfg.fz=malloc(sizeof(double)*n_total_particles);

    iccp3m_store_forces(); /* Not to over-write previous forces */
  /* Now iterate for getting charges */
      if((iccp3m_cfg.e1 == 0) && (iccp3m_cfg.e2 ==0)) {
          del_eps=0.0;
         } else {
          del_eps=(iccp3m_cfg.e1-iccp3m_cfg.e2)/(iccp3m_cfg.e2+iccp3m_cfg.e1);
        }
      iccp3m_cfg.citeration=0;
  for(j=0;j<iccp3m_cfg.num_iteration;j++) {
         force_calc_iccp3m(); /* Calculate electrostatic forces excluding source source interaction and short range contributions*/
      diff=0; 
    for(c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      part = cell->part;
      np   = cell->n;
      for(i=0 ; i < np; i++) {
          if(part[i].p.identity <= iccp3m_cfg.last_ind_id) {   /* interface grids */
               fdot=part[i].f.f[0]*iccp3m_cfg.nvectorx[part[i].p.identity]+
                part[i].f.f[1]*iccp3m_cfg.nvectory[part[i].p.identity]+part[i].f.f[2]*iccp3m_cfg.nvectorz[part[i].p.identity];
                hold=part[i].p.q/iccp3m_cfg.areas[part[i].p.identity];
                qold=part[i].p.q;
		hnew=hold;
              if(fabs(hold) > 1e-7) {
                fdot=fdot/(6.283185307); /* divide fdot with 2pi */
                hnew=hold + iccp3m_cfg.relax*((del_eps*fdot)/qold-hold);
              }
                difftemp=fabs(hold-hnew);
                if(difftemp > diff) { diff=difftemp; } /* Take the largest error for convergence */
                part[i].p.q=hnew*iccp3m_cfg.areas[part[i].p.identity];
                  if(fabs(part[i].p.q) > 100) { /* this is kind a arbitraru measure but, does a good job spotting divergence !*/
                    char *errtxt = runtime_error(128 + 2*TCL_DOUBLE_SPACE);
      	            ERROR_SPRINTF(errtxt, "{error occured 990 : too big charge assignment in iccp3m! q >100 , normal vectors or computed forces might be wrong or too big! assigned charge= %f } \n",part[i].p.q);
                     break;
                   }
           }
       }  /* cell particles */
    } /* local cells */
      /* printf(" iccp3m iteration j= %d convergence_cre = %f \r",j,diff);*/
      iccp3m_cfg.citeration++;
      if(diff < iccp3m_cfg.convergence) { 
	/* printf("ICCP3M converged step=%d \n",j); */
	break;
      }
  } /* iteration */
         free(iccp3m_cfg.areas);
    iccp3m_revive_forces(); /* revive originally computed forces by Espresso */
         free(iccp3m_cfg.fx);
         free(iccp3m_cfg.fy);
         free(iccp3m_cfg.fz); 
         free(iccp3m_cfg.nvectorx);
         free(iccp3m_cfg.nvectory);
         free(iccp3m_cfg.nvectorz);
      /* printf("ICCP3M_ITERATION iccp3m iteration finished \n");*/
  return iccp3m_cfg.citeration++;
}

void force_calc_iccp3m() {

#ifdef DIPOLES
   convert_quat_to_dip_all();
#endif
  init_forces_iccp3m();
  switch (cell_structure.type) {
  case CELL_STRUCTURE_LAYERED:
    layered_calculate_ia_iccp3m();
    break;
  case CELL_STRUCTURE_DOMDEC:
    if(dd.use_vList) {
      if (rebuild_verletlist) {
        build_verlet_lists_and_calc_verlet_ia_iccp3m();
       } else  {
        calculate_verlet_ia_iccp3m();
      }
    }
    else
      calc_link_cell_iccp3m();
    break;
  case CELL_STRUCTURE_NSQUARE:
    nsq_calculate_ia_iccp3m();
  }

  calc_long_range_forces_iccp3m();
#ifdef COMFORCE
 calc_comforce();
#endif

/* this must be the last force to be calculated (Mehmet)*/
#ifdef COMFIXED
 calc_comfixed();
#endif

}

/** nonbonded and bonded force calculation using the verlet list */
void layered_calculate_ia_iccp3m()
{
  int c, i, j;
  Cell  *celll, *cellb;
  int      npl,    npb;
  Particle *pl,    *pb, *p1;
  double dist2, d[3];

  CELL_TRACE(fprintf(stderr, "%d: rebuild_v=%d\n", this_node, rebuild_verletlist));

  for (c = 1; c <= n_layers; c++) {
    celll = &cells[c];
    pl    = celll->part;
    npl   = celll->n;

    cellb = &cells[c-1];
    pb    = cellb->part;
    npb   = cellb->n;

    for(i = 0; i < npl; i++) {
      p1 = &pl[i];

      if (rebuild_verletlist)
        memcpy(p1->l.p_old, p1->r.p, 3*sizeof(double));

      /* cell itself and bonded / constraints */
      for(j = i+1; j < npl; j++) {
        layered_get_mi_vector(d, p1->r.p, pl[j].r.p);
        dist2 = sqrlen(d);
#ifdef EXCLUSIONS
        if (do_nonbonded(p1, &pl[j])) {
#endif
          /* avoid source-source computation */
         if(!(p1->p.identity > iccp3m_cfg.last_ind_id && pl[j].p.identity >iccp3m_cfg.last_ind_id)) {
           add_non_bonded_pair_force_iccp3m(p1, &pl[j], d, sqrt(dist2), dist2);
          }
#ifdef EXCLUSIONS  
         }
#endif
      }

      /* bottom neighbor */
      for(j = 0; j < npb; j++) {
        layered_get_mi_vector(d, p1->r.p, pb[j].r.p);
        dist2 = sqrlen(d);
#ifdef EXCLUSIONS
        if (do_nonbonded(p1, &pl[j])) {
#endif
             /* avoid source-source computation */
         if(!(p1->p.identity > iccp3m_cfg.last_ind_id && pb[j].p.identity >iccp3m_cfg.last_ind_id)) {
           add_non_bonded_pair_force_iccp3m(p1, &pb[j], d, sqrt(dist2), dist2);
          }
#ifdef EXCLUSIONS
         }
#endif
      }
    }
  }
  rebuild_verletlist = 0;
}

void build_verlet_lists_and_calc_verlet_ia_iccp3m()
{
  int c, np1, n, np2, i ,j, j_start=0;
  Cell *cell;
  IA_Neighbor *neighbor;
  Particle *p1, *p2;
  PairList *pl;
  double dist2, vec21[3];
 
#ifdef VERLET_DEBUG 
  int estimate, sum=0;
  fprintf(stderr,"%d: build_verlet_list_and_calc_verlet_ia:\n",this_node);
  /* estimate number of interactions: (0.5*n_part*ia_volume*density)/n_nodes */
  estimate = 0.5*n_total_particles*(4.0/3.0*PI*pow(max_range_non_bonded,3.0))*(n_total_particles/(box_l[0]*box_l[1]*box_l[2]))/n_nodes;

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
	  memcpy(p1[i].l.p_old, p1[i].r.p, 3*sizeof(double));
	  j_start = i+1;
	}
	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {
#ifdef EXCLUSIONS
          if(do_nonbonded(&p1[i], &p2[j]))
#endif
	  {
	  dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);

	  VERLET_TRACE(fprintf(stderr,"%d: pair %d %d has distance %f\n",this_node,p1[i].p.identity,p2[j].p.identity,sqrt(dist2)));
	  if(dist2 <= max_range_non_bonded2) {
	    ONEPART_TRACE(if(p1[i].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d,%d %d,%d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,i,n,j,sqrt(dist2)));
	    ONEPART_TRACE(if(p2[j].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d %d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,n,sqrt(dist2)));

	    add_pair_iccp3m(pl, &p1[i], &p2[j]);
	    /* calc non bonded interactions */ 
           if(!(p1[i].p.identity > iccp3m_cfg.last_ind_id && p2[j].p.identity >iccp3m_cfg.last_ind_id)) {
	       add_non_bonded_pair_force_iccp3m(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
                } else {
                   printf("Ever here \n");
             }
             /* avoid source-source computation */ 
	  }
	 }
	}
      }
      resize_verlet_list_iccp3m(pl);
      VERLET_TRACE(fprintf(stderr,"%d: neighbor %d has %d pairs\n",this_node,n,pl->n));
      VERLET_TRACE(sum += pl->n);
    }
  }

  VERLET_TRACE(fprintf(stderr,"%d: total number of interaction pairs: %d (should be around %d)\n",this_node,sum,estimate));
 
  rebuild_verletlist = 0;
}

void calculate_verlet_ia_iccp3m()
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
       if(!(p1->p.identity > iccp3m_cfg.last_ind_id && p2->p.identity >iccp3m_cfg.last_ind_id))  /* avoid source-source computation */
	  add_non_bonded_pair_force_iccp3m(p1, p2, vec21, sqrt(dist2), dist2);
      }
    }
  }
}

void calc_link_cell_iccp3m()
{
  int c, np1, n, np2, i ,j, j_start=0;
  Cell *cell;
  IA_Neighbor *neighbor;
  Particle *p1, *p2;
  double dist2, vec21[3];
 
  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {

    cell = local_cells.cell[c];
    p1   = cell->part;
    np1  = cell->n;
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      neighbor = &dd.cell_inter[c].nList[n];
      p2  = neighbor->pList->part;
      np2 = neighbor->pList->n;
      /* Loop cell particles */
      for(i=0; i < np1; i++) {
        j_start = 0;
        /* Tasks within cell: bonded forces */
        if(n == 0) {
          add_bonded_force(&p1[i]);
#ifdef CONSTRAINTS
          add_constraints_forces(&p1[i]);
#endif
          j_start = i+1;
        }
	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {
	    {
	      dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);
	      if(dist2 <= max_range_non_bonded2) {
		/* calc non bonded interactions */
                if(!(p1[i].p.identity > iccp3m_cfg.last_ind_id && p2[j].p.identity >iccp3m_cfg.last_ind_id)) 
		add_non_bonded_pair_force_iccp3m(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
                 /* avoid source-source computation */
	      }
	    }
	}
      }
    }
  }
}

void nsq_calculate_ia_iccp3m()
{
  Particle *partl, *partg;
  Particle *pt1, *pt2;
  int p, p2, npl, npg, c;
  double d[3], dist2, dist;

  npl   = local->n;
  partl = local->part;

  /* calculate bonded interactions and non bonded node-node */
  for (p = 0; p < npl; p++) {
    pt1 = &partl[p];
    /* other particles, same node */
    for (p2 = p + 1; p2 < npl; p2++) {
      pt2 = &partl[p2];
      get_mi_vector(d, pt1->r.p, pt2->r.p);
      dist2 = sqrlen(d);
      dist = sqrt(dist2);
          /* avoid source-source computation */
         if(!(pt1->p.identity > iccp3m_cfg.last_ind_id && pt2->p.identity >iccp3m_cfg.last_ind_id))
	 add_non_bonded_pair_force_iccp3m(pt1, pt2, d, dist, dist2);
    }

    /* calculate with my ghosts */
    for (c = 0; c < me_do_ghosts.n; c++) {
      npg   = me_do_ghosts.cell[c]->n;
      partg = me_do_ghosts.cell[c]->part;

      for (p2 = 0; p2 < npg; p2++) {
	pt2 = &partg[p2];
	get_mi_vector(d, pt1->r.p, pt2->r.p);
	dist2 = sqrlen(d);
	dist = sqrt(dist2);
           /* avoid source-source computation */
        if(!(pt1->p.identity > iccp3m_cfg.last_ind_id && pt2->p.identity >iccp3m_cfg.last_ind_id))
	add_non_bonded_pair_force_iccp3m(pt1, pt2, d, dist, dist2);
      }
    }
  }
}


void init_forces_iccp3m()
{
  Cell *cell;
  Particle *p;
  int np, c, i;

  /* The force initialization depends on the used thermostat and the
     thermodynamic ensemble */

#ifdef NPT
  /* reset virial part of instantaneous pressure */
  if(integ_switch == INTEG_METHOD_NPT_ISO)
    nptiso.p_vir[0] = nptiso.p_vir[1] = nptiso.p_vir[2] = 0.0;
#endif


  /* initialize forces with langevin thermostat forces
     or zero depending on the thermostat
     set torque to zero for all and rescale quaternions
  */
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++)
      init_local_particle_force_iccp3m(&p[i]);
  }

  /* initialize ghost forces with zero
     set torque to zero for all and rescale quaternions
  */
  for (c = 0; c < ghost_cells.n; c++) {
    cell = ghost_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++)
      init_ghost_force_iccp3m(&p[i]);
  }
   
}

void calc_long_range_forces_iccp3m()
{
#ifdef ELECTROSTATICS
  /* calculate k-space part of electrostatic interaction. */
  switch (coulomb.method) {
  case COULOMB_P3M:
      P3M_charge_assign();
 #ifdef NPT
     if(integ_switch == INTEG_METHOD_NPT_ISO) 
       nptiso.p_vir[0] += P3M_calc_kspace_forces(1,1);
#endif
      P3M_calc_kspace_forces(1,0);
    break;
  }

#endif
}
/** \name Privat Functions */
/************************************************************/
/*@{*/

/** Add a particle pair to a verlet pair list.
    Checks verlet pair list size and reallocates memory if necessary.
 *  \param p1 Pointer to paricle one.
 *  \param p2 Pointer to paricle two.
 *  \param pl Pointer to the verlet pair list.
 */
MDINLINE void add_pair_iccp3m(PairList *pl, Particle *p1, Particle *p2)
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

void resize_verlet_list_iccp3m(PairList *pl)
{
  int diff;
  diff = pl->max - pl->n;
  if( diff > 2*LIST_INCREMENT ) {
    diff = (diff/LIST_INCREMENT)-1;
    pl->max -= diff*LIST_INCREMENT;
    pl->pair = (Particle **)realloc(pl->pair, 2*pl->max*sizeof(Particle *));
  }
}

/** initialize the forces for a real particle */
MDINLINE void init_local_particle_force_iccp3m(Particle *part)
{

    part->f.f[0] = 0.0; /* no need to friction_thermo_langevin function */
    part->f.f[1] = 0.0;
    part->f.f[2] = 0.0;

}

/** initialize the forces for a ghost particle */
MDINLINE void init_ghost_force_iccp3m(Particle *part)
{
  part->f.f[0] = 0.0;
  part->f.f[1] = 0.0;
  part->f.f[2] = 0.0;
}

/* integer mod*/
int imod(int x,int y) {
  double p,q,m;
  int z;
   p=x; q=y;
   m=fmod(p,q);
   z=m;
  return z;
}

int parse_area(Tcl_Interp *interp,int normal_args, char *string) {

   int size,i;
   const char delimiters[] = " ";
   char *token,*cp;
   float temp;

   size= normal_args;
   iccp3m_cfg.areas=malloc((size+1)*sizeof(double));
   cp = strdup(string);                /* Make writable copy.  */
   token = strtok (cp, delimiters);
     sscanf(token,"%f",&temp);
     iccp3m_cfg.areas[0]=temp;

    for(i=1;i<size;i++) {
      token = strtok (NULL, delimiters);
          if(token == NULL){
             return(1);
            }
       sscanf(token,"%f",&temp);
      iccp3m_cfg.areas[i]=temp;
    }
  return(0);
}

int parse_normal(Tcl_Interp *interp,int normal_args, char *string) {

   int size,i,k;
   double *numbers;
   const char delimiters[] = " ";
   char *token,*cp;
   float temp;

   size= normal_args;
   numbers=malloc((size+1)*sizeof(double));
   cp = strdup(string);                /* Make writable copy.  */
   token = strtok (cp, delimiters);
    sscanf(token,"%f",&temp);
     numbers[0]=temp;
     
    for(i=1;i<size;i++) {
      token = strtok (NULL, delimiters);
          if(token == NULL){
             return(1);
            }
       sscanf(token,"%f",&temp);
      numbers[i]=temp;
    }
    iccp3m_cfg.nvectorx=malloc(sizeof(double)*(iccp3m_cfg.last_ind_id+1));
    iccp3m_cfg.nvectory=malloc(sizeof(double)*(iccp3m_cfg.last_ind_id+1));
    iccp3m_cfg.nvectorz=malloc(sizeof(double)*(iccp3m_cfg.last_ind_id+1));

   /* assign normal vectors */
     k=0;
    for(i=0;i<size;i++) {
       if( (i == 0) || (imod(i,3) == 0) ) {iccp3m_cfg.nvectorx[k]=numbers[i]; } /* assign components */
       if( (i == 1) || (imod(i,3) == 1) ) {iccp3m_cfg.nvectory[k]=numbers[i]; }
       if( (i == 2) || (imod(i,3) == 2) ) {iccp3m_cfg.nvectorz[k]=numbers[i];  k++; } 
     }

   free(numbers);
  return(0);
}
void reset_forces() {
     Cell *cell;    
     int c,i,np;
     Particle *part;
   for(c = 0; c < local_cells.n; c++) {
       cell = local_cells.cell[c];
       part = cell->part;
       np   = cell->n;
     for(i=0 ; i < np; i++) {
      part[i].f.f[0]=0.0;  part[i].f.f[1]=0.0;  part[i].f.f[2]=0.0;
     }
   }
}
void iccp3m_revive_forces() {
  /* restore forces that are computed before in Espresso integrate_vv function*/
   Cell *cell;
   int c,i,np;
   Particle *part;
   for(c = 0; c < local_cells.n; c++) {
     cell = local_cells.cell[c];
     part = cell->part;
     np   = cell->n;
    for(i=0 ; i < np; i++) {
	part[i].f.f[0]=iccp3m_cfg.fx[part[i].p.identity];
	part[i].f.f[1]=iccp3m_cfg.fy[part[i].p.identity];
	part[i].f.f[2]=iccp3m_cfg.fz[part[i].p.identity];
     }
   }
}
void iccp3m_store_forces() {
  /* store forces that are computed before in Espresso integrate_vv function */
  /* iccp3m will re-compute electrostatic-forces on boundary particles */
   Cell *cell;
   int c,i,np;
   Particle *part;
   for(c = 0; c < local_cells.n; c++) {
     cell = local_cells.cell[c];
     part = cell->part;
     np   = cell->n;
    for(i=0 ; i < np; i++) {
       iccp3m_cfg.fx[part[i].p.identity]=part[i].f.f[0];
       iccp3m_cfg.fy[part[i].p.identity]=part[i].f.f[1];
       iccp3m_cfg.fz[part[i].p.identity]=part[i].f.f[2];
     }
   }
}

int inter_parse_iccp3m(Tcl_Interp * interp, int argc, char ** argv) {
   char buffer[TCL_DOUBLE_SPACE];
   Tcl_PrintDouble(interp, iccp3m_cfg.citeration, buffer);
   Tcl_AppendResult(interp, buffer,(char *) NULL);
  return TCL_OK;
}
int gettime(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
    char buffer[TCL_DOUBLE_SPACE];
    clock_t cputime;
    double timing;
    cputime= clock();
    timing=((double) cputime) / CLOCKS_PER_SEC;
    Tcl_PrintDouble(interp, timing, buffer);
    Tcl_AppendResult(interp, buffer,(char *) NULL);
    return TCL_OK;
}

#endif

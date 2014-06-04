/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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

/** \file iccp3m.cpp
    Detailed Information about the method is included in the corresponding header file \ref iccp3m.hpp.

 */
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <cstddef>
#include <cmath>
#include <ctime>

#include "iccp3m.hpp"
#include "p3m.hpp"
#include "elc.hpp"
#include "mmm2d.hpp"
#include "mmm1d.hpp"

#include "communication.hpp"

#include "utils.hpp"
#include "verlet.hpp"
#include "cells.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "domain_decomposition.hpp"
#include "verlet.hpp"
#include "forces.hpp"
#include "config.hpp"
#include "global.hpp"

#ifdef ELECTROSTATICS

iccp3m_struct iccp3m_cfg;

int iccp3m_initialized = 0;
/* functions that are used in icc* to compute the electric field acting on the induced charges,
 * excluding forces other than the electrostatic ones */
void init_forces_iccp3m();
void calc_long_range_forces_iccp3m();
inline void add_pair_iccp3m(PairList *pl, Particle *p1, Particle *p2);
void resize_verlet_list_iccp3m(PairList *pl);
inline void init_local_particle_force_iccp3m(Particle *part);
inline void init_ghost_force_iccp3m(Particle *part);
extern void on_particle_change();

Cell *local_icc;
CellPList me_do_ghosts_icc;

/* other icc*  functions */
int imod(int x,int y);
void iccp3m_revive_forces();
void iccp3m_store_forces();

/** Granularity of the verlet list */
#define LIST_INCREMENT 20

void iccp3m_init(void){
   iccp3m_cfg.set_flag=0;
   iccp3m_cfg.areas = NULL;
   iccp3m_cfg.ein = NULL;
   iccp3m_cfg.nvectorx = NULL;
   iccp3m_cfg.nvectory = NULL;
   iccp3m_cfg.nvectorz = NULL;
   iccp3m_cfg.extx = 0;
   iccp3m_cfg.exty = 0;
   iccp3m_cfg.extz = 0;
   iccp3m_cfg.first_id = 0;
   iccp3m_cfg.num_iteration=30;
   iccp3m_cfg.convergence=1e-2;
   iccp3m_cfg.relax=0.7;
   iccp3m_cfg.eout=1;
   iccp3m_cfg.citeration=0;


}



int bcast_iccp3m_cfg(void){
  int i;


  MPI_Bcast((int*)&iccp3m_cfg.n_ic, 1, MPI_INT, 0, comm_cart); 

  /* allocates Memory on slave nodes 
   * Master node allocates the memory when parsing tcl arguments
   * */
  if (this_node != 0) {
    iccp3m_cfg.areas      = (double*) realloc (iccp3m_cfg.areas     ,(iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.ein        = (double*) realloc (iccp3m_cfg.ein       ,(iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.nvectorx   = (double*) realloc (iccp3m_cfg.nvectorx  ,(iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.nvectory   = (double*) realloc (iccp3m_cfg.nvectory  ,(iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.nvectorz   = (double*) realloc (iccp3m_cfg.nvectorz  ,(iccp3m_cfg.n_ic) * sizeof(double));
    iccp3m_cfg.sigma      = (double*) realloc (iccp3m_cfg.sigma     ,(iccp3m_cfg.n_ic) * sizeof(double));
  }

  MPI_Bcast((int*)&iccp3m_cfg.num_iteration, 1, MPI_INT, 0, comm_cart); 
  MPI_Bcast((int*)&iccp3m_cfg.first_id, 1, MPI_INT, 0, comm_cart); 
  MPI_Bcast((double*)&iccp3m_cfg.convergence, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast((double*)&iccp3m_cfg.eout, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast((double*)&iccp3m_cfg.relax, 1, MPI_DOUBLE, 0, comm_cart);
  
  /* broadcast the vectors element by element. This is slow
   * but safe and only performed at the beginning of each simulation*/
  for ( i = 0; i < iccp3m_cfg.n_ic; i++) {
    MPI_Bcast((double*)&iccp3m_cfg.areas[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double*)&iccp3m_cfg.ein[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double*)&iccp3m_cfg.nvectorx[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double*)&iccp3m_cfg.nvectory[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double*)&iccp3m_cfg.nvectorz[i], 1, MPI_DOUBLE, 0, comm_cart);
    MPI_Bcast((double*)&iccp3m_cfg.sigma[i], 1, MPI_DOUBLE, 0, comm_cart);
  }
  MPI_Bcast((double*)&iccp3m_cfg.extx, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast((double*)&iccp3m_cfg.exty, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast((double*)&iccp3m_cfg.extz, 1, MPI_DOUBLE, 0, comm_cart);

  MPI_Bcast(&iccp3m_cfg.citeration, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&iccp3m_cfg.set_flag, 1, MPI_DOUBLE, 0, comm_cart);

  return 0 ;
    
}

int iccp3m_iteration() {
   double fdot,hold,hnew,hmax,del_eps,diff=0.0,difftemp=0.0, ex, ey, ez, l_b;
   Cell *cell;
   int c,np;
   Particle *part;
   int i, j,id;
   char* errtxt;
   double globalmax;
   double f1, f2 = 0;

   iccp3m_sanity_check();
   
   l_b = coulomb.bjerrum;
   if((iccp3m_cfg.eout <= 0)) {
     errtxt = runtime_error(128);
     ERROR_SPRINTF(errtxt, "ICCP3M: nonpositive dielectric constant is not allowed. Put a decent tcl error here\n");
   }
   
   
   iccp3m_cfg.citeration=0;
   for(j=0;j<iccp3m_cfg.num_iteration;j++) {
     hmax=0.;
     force_calc_iccp3m(); /* Calculate electrostatic forces (SR+LR) excluding source source interaction*/
     diff=0;
     for(c = 0; c < local_cells.n; c++) {
       cell = local_cells.cell[c];
       part = cell->part;
       np   = cell->n;
       for(i=0 ; i < np; i++) {
         if( part[i].p.identity < iccp3m_cfg.n_ic+iccp3m_cfg.first_id && part[i].p.identity >= iccp3m_cfg.first_id ) {
           id = part[i].p.identity - iccp3m_cfg.first_id;
           /* the dielectric-related prefactor: */                     
           del_eps = (iccp3m_cfg.ein[id]-iccp3m_cfg.eout)/(iccp3m_cfg.ein[id] + iccp3m_cfg.eout)/6.283185307;
           /* calculate the electric field at the certain position */
           ex=part[i].f.f[0]/part[i].p.q;
           ey=part[i].f.f[1]/part[i].p.q;
           ez=part[i].f.f[2]/part[i].p.q;
           
           /* let's add the contribution coming from the external field */
           ex += iccp3m_cfg.extx; 
           ey += iccp3m_cfg.exty; 
           ez += iccp3m_cfg.extz;
           
           if (ex == 0 && ey == 0 && ez == 0) {
             errtxt = runtime_error(128);
             ERROR_SPRINTF(errtxt, "ICCP3M found zero electric field on a charge. This must never happen");
           }
           /* the dot product   */
           fdot = ex*iccp3m_cfg.nvectorx[id]+
             ey*iccp3m_cfg.nvectory[id]+
             ez*iccp3m_cfg.nvectorz[id];
           
           /* recalculate the old charge density */
                      hold=part[i].p.q/iccp3m_cfg.areas[id];
          /* determine if it is higher than the previously highest charge density */            
                      if(fabs(hold)>hmax)hmax=fabs(hold); 
                      f1 =  (+del_eps*fdot/l_b);
//                      double f2 = (1- 0.5*(iccp3m_cfg.ein[id]-iccp3m_cfg.eout)/(iccp3m_cfg.eout + iccp3m_cfg.ein[id] ))*(iccp3m_cfg.sigma[id]);
                      if (iccp3m_cfg.sigma!=0) {
                        f2 = (2*iccp3m_cfg.eout)/(iccp3m_cfg.eout + iccp3m_cfg.ein[id] )*(iccp3m_cfg.sigma[id]);
                      } 

                      hnew=(1.-iccp3m_cfg.relax)*hold + (iccp3m_cfg.relax)*(f1 + f2);
                      difftemp=fabs( 1*(hnew - hold)/(hmax + fabs(hnew+hold)) ); /* relative variation: never use 
                                                                              an estimator which can be negative
                                                                              here */
                      if(difftemp > diff && part[i].p.q > 1e-5)
                      {
//                          if (fabs(difftemp - 1./(1./iccp3m_cfg.relax - 1.)) > 1e-10) 
                        diff=difftemp;  /* Take the largest error to check for convergence */
                      }
                      part[i].p.q = hnew * iccp3m_cfg.areas[id];
         /* check if the charge now is more than 1e6, to determine if ICC still leads to reasonable results */
         /* this is kind a arbitrary measure but, does a good job spotting divergence !*/
                      if(fabs(part[i].p.q) > 1e6) { 
                        errtxt = runtime_error(128);
            	          ERROR_SPRINTF(errtxt, "{error occured 990 : too big charge assignment in iccp3m! q >1e6 , \
                               assigned charge= %f } \n",part[i].p.q);
                        diff = 1e90; /* A very high value is used as error code */
                        break;
                      }
                 }
             }  /* cell particles */
           // printf("cell %d w %d particles over (node %d)\n",c,np,this_node); fflush(stdout);
       } /* local cells */
       iccp3m_cfg.citeration++;
       MPI_Allreduce(&diff, &globalmax, 1,MPI_DOUBLE, MPI_MAX, comm_cart);

       if (globalmax < iccp3m_cfg.convergence) 
         break; 
       if ( diff > 1e89 ) /* Error happened */
         return iccp3m_cfg.citeration++;

  } /* iteration */
  on_particle_change();

  return iccp3m_cfg.citeration;
}

void force_calc_iccp3m() {
/* The following ist mostly copied from forces.cpp */

/*  I don´t see the point of this part until there are electrical dipoles in Espresso, BTW, it generates a warning .. JJCP
#ifdef DIPOLES
   convert_quat_to_dip_all();
#endif

*/

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

      /* cell itself. No bonded / constraints considered in ICCP3M */
      for(j = i+1; j < npl; j++) {
        layered_get_mi_vector(d, p1->r.p, pl[j].r.p);
        dist2 = sqrlen(d);
#ifdef EXCLUSIONS
        if (do_nonbonded(p1, &pl[j])) {
#endif
          /* avoid source-source computation */
           add_non_bonded_pair_force_iccp3m(p1, &pl[j], d, sqrt(dist2), dist2);
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
           add_non_bonded_pair_force_iccp3m(p1, &pb[j], d, sqrt(dist2), dist2);
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
  estimate = 0.5*n_part*(4.0/3.0*PI*pow(max_cut_nonbonded,3.0))*(n_part/(box_l[0]*box_l[1]*box_l[2]))/n_nodes;

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
	/* Tasks within cell: (no bonded forces) store old position, avoid double counting */
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
	  dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);

	  VERLET_TRACE(fprintf(stderr,"%d: pair %d %d has distance %f\n",this_node,p1[i].p.identity,p2[j].p.identity,sqrt(dist2)));
	  if(dist2 <= SQR(get_ia_param(p1[i].p.type, p2[j].p.type)->max_cut + skin)) {
	    ONEPART_TRACE(if(p1[i].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d,%d %d,%d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,i,n,j,sqrt(dist2)));
	    ONEPART_TRACE(if(p2[j].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d %d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,n,sqrt(dist2)));

	    add_pair_iccp3m(pl, &p1[i], &p2[j]);
	    /* calc non bonded interactions */ 
	       add_non_bonded_pair_force_iccp3m(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
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
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      pairs = dd.cell_inter[c].nList[n].vList.pair;
      np    = dd.cell_inter[c].nList[n].vList.n;
      /* verlet list loop */
      for(i=0; i<2*np; i+=2) {
	      p1 = pairs[i];                    /* pointer to particle 1 */
	      p2 = pairs[i+1];                  /* pointer to particle 2 */
	      dist2 = distance2vec(p1->r.p, p2->r.p, vec21); 
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
	  if (rebuild_verletlist)
	    memcpy(p1[i].l.p_old, p1[i].r.p, 3*sizeof(double));

          j_start = i+1;
        }
	/* Loop neighbor cell particles */
	for(j = j_start; j < np2; j++) {
	    {
	      dist2 = distance2vec(p1[i].r.p, p2[j].r.p, vec21);
	      if(dist2 <= SQR(get_ia_param(p1[i].p.type, p2[j].p.type)->max_cut + skin)) {
		/* calc non bonded interactions */
		add_non_bonded_pair_force_iccp3m(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
	      }
	    }
	}
      }
    }
  }
  rebuild_verletlist = 0;
}

void nsq_calculate_ia_iccp3m()
{
  Particle *partl, *partg;
  Particle *pt1, *pt2;
  int p, p2, npl, npg, c;
  double d[3], dist2, dist;

  npl   = local_icc->n;
  partl = local_icc->part;

  /* calculate bonded interactions and non bonded node-node */
  for (p = 0; p < npl; p++) {
    pt1 = &partl[p];

    if (rebuild_verletlist)
      memcpy(pt1->l.p_old, pt1->r.p, 3*sizeof(double));

    /* other particles, same node */
    for (p2 = p + 1; p2 < npl; p2++) {
      pt2 = &partl[p2];
      get_mi_vector(d, pt1->r.p, pt2->r.p);
      dist2 = sqrlen(d);
      dist = sqrt(dist2);
	    add_non_bonded_pair_force_iccp3m(pt1, pt2, d, dist, dist2);
    }

    /* calculate with my ghosts */
    for (c = 0; c < me_do_ghosts_icc.n; c++) {
      npg   = me_do_ghosts_icc.cell[c]->n;
      partg = me_do_ghosts_icc.cell[c]->part;

      for (p2 = 0; p2 < npg; p2++) {
	      pt2 = &partg[p2];
	      get_mi_vector(d, pt1->r.p, pt2->r.p);
	      dist2 = sqrlen(d);
	      dist = sqrt(dist2);
	      add_non_bonded_pair_force_iccp3m(pt1, pt2, d, dist, dist2);
      }
    }
  }
  rebuild_verletlist = 0;
}


void init_forces_iccp3m()
{
  /* copied from forces.cpp */
  Cell *cell;
  Particle *p;
  int np, c, i;

  /* The force initialization depends on the used thermostat and the
     thermodynamic ensemble */

#ifdef NPT
  char* errtxt;
  /* reset virial part of instantaneous pressure */
  if(integ_switch == INTEG_METHOD_NPT_ISO){
      errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{ICCP3M cannot be used with pressure coupling} ");
  }
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
  char* errtxt;
  /* calculate k-space part of electrostatic interaction. */
  if (!(coulomb.method == COULOMB_ELC_P3M || 
        coulomb.method == COULOMB_P3M     || 
        coulomb.method == COULOMB_MMM2D   ||
        coulomb.method == COULOMB_MMM1D)  ) { 
                errtxt = runtime_error(128);
                ERROR_SPRINTF(errtxt, "{ICCP3M implemented only for MMM1D,MMM2D,ELC or P3M ");
     }
  switch (coulomb.method) {
#ifdef P3M
    case COULOMB_ELC_P3M:
       if (elc_params.dielectric_contrast_on) {
                errtxt = runtime_error(128);
                ERROR_SPRINTF(errtxt, "{ICCP3M conflicts with ELC dielectric constrast");
       }
       p3m_charge_assign();
       p3m_calc_kspace_forces(1,0);
       ELC_add_force(); 
    break;

    case COULOMB_P3M:
       p3m_charge_assign();
       p3m_calc_kspace_forces(1,0);
    break;
#endif
    case COULOMB_MMM2D:
       MMM2D_add_far_force();
       MMM2D_dielectric_layers_force_contribution();
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
inline void add_pair_iccp3m(PairList *pl, Particle *p1, Particle *p2)
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
inline void init_local_particle_force_iccp3m(Particle *part)
{
    part->f.f[0] = 0.0; /* no need to friction_thermo_langevin function */
    part->f.f[1] = 0.0;
    part->f.f[2] = 0.0;

#ifdef ROTATION
    /* set torque to zero */
    part->f.torque[0] = 0;
    part->f.torque[1] = 0;
    part->f.torque[2] = 0;
#endif
}

/** initialize the forces for a ghost particle */
inline void init_ghost_force_iccp3m(Particle *part)
{
  part->f.f[0] = 0.0;
  part->f.f[1] = 0.0;
  part->f.f[2] = 0.0;

#ifdef ROTATION
  /* set torque to zero */
  part->f.torque[0] = 0;
  part->f.torque[1] = 0;
  part->f.torque[2] = 0;
#endif
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

int iccp3m_sanity_check()
{
  switch (coulomb.method) {
#ifdef P3M
    case COULOMB_ELC_P3M: {
      if (elc_params.dielectric_contrast_on) {
	char *errtxt = runtime_error(128);
	ERROR_SPRINTF(errtxt, "ICCP3M conflicts with ELC dielectric constrast");
	return 1;
      }
      break;
    }
#endif
    case COULOMB_DH: {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "ICCP3M does not work with Debye-Hueckel iccp3m.h");
      return 1;
    }
    case COULOMB_RF: {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "ICCP3M does not work with COULOMB_RF iccp3m.h");
      return 1;
    }
  }
  
#ifdef NPT
  if(integ_switch == INTEG_METHOD_NPT_ISO) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "ICCP3M does not work in the NPT ensemble");
    return 1;
  }
#endif

  return 0;
}

#endif


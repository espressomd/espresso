/*
  Copyright (C) 2010 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany
  
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

/** \file iccp3m.c
    Detailed Information about the method is included in the corresponding header file \ref iccp3m.h.

 */
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <stddef.h>
#include <math.h>
#include <time.h>

#include "iccp3m.h"
#include "p3m.h"
#include "elc.h"
#include "mmm2d.h"
#include "mmm1d.h"

#include "communication.h"

#include "utils.h"
#include "tcl.h"
#include "parser.h"
#include "verlet.h"
#include "cells.h"
#include "particle_data.h"
#include "domain_decomposition.h"
#include "verlet.h"
#include "forces.h"
#include "config.h"
#include "global.h"

#ifdef ELECTROSTATICS

enum { ICCP3M_AREA , ICCP3M_EPSILON, ICCP3M_NORMAL, ICCP3M_EXTFIELD } ;
iccp3m_struct iccp3m_cfg;

int iccp3m_initialized = 0;
/* functions that are used in icc* to compute the electric field acting on the induced charges,
 * excluding forces other than the electrostatic ones */
void init_forces_iccp3m();
void calc_long_range_forces_iccp3m();
MDINLINE void add_pair_iccp3m(PairList *pl, Particle *p1, Particle *p2);
void resize_verlet_list_iccp3m(PairList *pl);
MDINLINE void init_local_particle_force_iccp3m(Particle *part);
MDINLINE void init_ghost_force_iccp3m(Particle *part);
extern void on_particle_change();

Cell *local_icc;
CellPList me_do_ghosts_icc;

/* other icc*  functions */
static int tclcommand_iccp3m_parse_params(Tcl_Interp *interp,int normal_args, char *string, int flag);
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
   iccp3m_cfg.extx = NULL;
   iccp3m_cfg.exty = NULL;
   iccp3m_cfg.extz = NULL;
   iccp3m_cfg.fx = NULL;
   iccp3m_cfg.fy = NULL;
   iccp3m_cfg.fz = NULL;
}

/** Parses the ICCP3M command.
 */
int tclcommand_iccp3m(ClientData data, Tcl_Interp *interp, int argc, char **argv) {
  int last_ind_id,num_iteration,normal_args,area_args;
  char buffer[TCL_DOUBLE_SPACE];
  double e1,convergence,relax;

  Tcl_AppendResult(interp, "The ICCP3M algorithm is still experimental. Function can not be guaranteed, therefore it is still disabled.\n", (char *)NULL);
  return (TCL_ERROR);

  if(iccp3m_initialized==0){
      iccp3m_init();
      iccp3m_initialized=1;
  }

  if(argc != 9 && argc != 2 && argc != 10) { 
         Tcl_AppendResult(interp, "Wrong # of args! Usage: iccp3m { iterate | <last_ind_id> <e1> <num_iteration> <convergence> <relaxation> <area> <normal_components> <e_in/e_out>  [<ext_field>] }", (char *)NULL); 
         return (TCL_ERROR); 
   }
   if (argc == 2 ){
      if(ARG_IS_S(1,"iterate")) { 
           if (iccp3m_cfg.set_flag==0) {
                 Tcl_AppendResult(interp, "iccp3m parameters not set!", (char *)NULL);
                 return (TCL_ERROR);
           }
           else{ 
              Tcl_PrintDouble(interp,mpi_iccp3m_iteration(0),buffer); 
              Tcl_AppendResult(interp, buffer, (char *) NULL);
              return TCL_OK;
	   }
      }
   }
   else {
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
       if(!ARG_IS_I(3, num_iteration)) {
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "Number of maximum iterations must be integer (got: ", argv[4],")!", (char *)NULL); return (TCL_ERROR);
       }
       if(!ARG_IS_D(4, convergence)) {
          Tcl_ResetResult(interp);
          Tcl_AppendResult(interp, "Convergence criterion must be double (got: ", argv[5],")!", (char *)NULL); return (TCL_ERROR);
          } else if( (convergence < 1e-10) || (convergence > 1e-1) ) { 
            Tcl_ResetResult(interp);
            Tcl_AppendResult(interp, "Convergence criterion can not be smaller then 1e-10 and greater then 1e-2(got: ", argv[5],")!", (char *)NULL); return (TCL_ERROR);
         }
          if(!ARG_IS_D(5, relax)) {
             Tcl_ResetResult(interp);
             Tcl_AppendResult(interp, "Relaxation parameter must be double (got: ", argv[6],")!", (char *)NULL); return (TCL_ERROR);
       }
   
       iccp3m_cfg.last_ind_id = last_ind_id;      /* Assign needed options */
       iccp3m_cfg.num_iteration = num_iteration; /* maximum number of iterations to go */
       iccp3m_cfg.eout = e1;
       iccp3m_cfg.convergence = convergence;
       iccp3m_cfg.relax = relax;
       iccp3m_cfg.update = 0;
       iccp3m_cfg.set_flag = 1;
       
       normal_args = (iccp3m_cfg.last_ind_id+1)*3;
       /* Now get Normal Vectors components consecutively */
       
       if( tclcommand_iccp3m_parse_params(interp,normal_args,argv[7],ICCP3M_NORMAL) == 1) { 
              Tcl_ResetResult(interp);
              Tcl_AppendResult(interp, "ICCP3M: Error in following normal vectors\n", argv[7],"\nICCP3M: Error in previous normal vectors\n", (char *)NULL); 
              return (TCL_ERROR);
       }      

       area_args=(iccp3m_cfg.last_ind_id+1);
       /* Now get area of the boundary elements */
       
       if ( tclcommand_iccp3m_parse_params(interp,area_args,argv[6],ICCP3M_AREA) == 1 ){
             Tcl_ResetResult(interp);
             Tcl_AppendResult(interp, "ICCP3M: Error in following areas\n", argv[6],"\nICCP3M: Error in previous areas\n", (char *)NULL); return (TCL_ERROR);
       }
       /* Now get the ration ein/eout for each element. 
          It's the user's duty to make sure that only disconnected 
          regions have different ratios */
      if ( tclcommand_iccp3m_parse_params(interp,area_args,argv[8],ICCP3M_EPSILON) == 1 ) {
             Tcl_ResetResult(interp);
             Tcl_AppendResult(interp, "ICCP3M: Error in following dielectric constants\n", argv[8],"\nICCP3M:  Error in previous dielectric constants\n", (char *)NULL); return (TCL_ERROR);
       } 

       if( argc == 10 ) {
         if( tclcommand_iccp3m_parse_params(interp,normal_args,argv[9],ICCP3M_EXTFIELD) == 1) { 
              Tcl_ResetResult(interp);
              Tcl_AppendResult(interp, "ICCP3M: Error in following external field vectors\n", argv[9],"\nICCP3M: Error in previous external field vectors\n", (char *)NULL); return (TCL_ERROR);
         }      
       }
       else {
         printf("allocating zeroes for external field \n");
         iccp3m_cfg.extx = (double*)calloc((last_ind_id +1), sizeof(double));
         iccp3m_cfg.exty = (double*)calloc((last_ind_id +1), sizeof(double));
         iccp3m_cfg.extz = (double*)calloc((last_ind_id +1), sizeof(double));
       }
      
       mpi_iccp3m_init(0);
       Tcl_PrintDouble(interp,mpi_iccp3m_iteration(0),buffer); 
       Tcl_AppendResult(interp, buffer, (char *) NULL);
       return TCL_OK;
   } /* else (argc==10) */
   return TCL_OK;
}


int bcast_iccp3m_cfg(void){
  int i;


  MPI_Bcast((int*)&iccp3m_cfg.last_ind_id, 1, MPI_INT, 0, MPI_COMM_WORLD); 

  /* allocates Memory on slave nodes 
   * Master node allocates the memory when parsing tcl arguments
   * */
  if (this_node != 0) {
    iccp3m_cfg.areas      = (double*) realloc (iccp3m_cfg.areas     ,(iccp3m_cfg.last_ind_id+1) * sizeof(double));
    iccp3m_cfg.ein        = (double*) realloc (iccp3m_cfg.ein       ,(iccp3m_cfg.last_ind_id+1) * sizeof(double));
    iccp3m_cfg.nvectorx   = (double*) realloc (iccp3m_cfg.nvectorx  ,(iccp3m_cfg.last_ind_id+1) * sizeof(double));
    iccp3m_cfg.nvectory   = (double*) realloc (iccp3m_cfg.nvectory  ,(iccp3m_cfg.last_ind_id+1) * sizeof(double));
    iccp3m_cfg.nvectorz   = (double*) realloc (iccp3m_cfg.nvectorz  ,(iccp3m_cfg.last_ind_id+1) * sizeof(double));
    iccp3m_cfg.extx       = (double*) realloc (iccp3m_cfg.extx      ,(iccp3m_cfg.last_ind_id+1) * sizeof(double));
    iccp3m_cfg.exty       = (double*) realloc (iccp3m_cfg.exty      ,(iccp3m_cfg.last_ind_id+1) * sizeof(double));
    iccp3m_cfg.extz       = (double*) realloc (iccp3m_cfg.extz      ,(iccp3m_cfg.last_ind_id+1) * sizeof(double));
  }


  MPI_Bcast((int*)&iccp3m_cfg.num_iteration, 1, MPI_INT, 0, MPI_COMM_WORLD); 
  MPI_Bcast((double*)&iccp3m_cfg.convergence, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast((double*)&iccp3m_cfg.eout, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast((double*)&iccp3m_cfg.relax, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast((int*)&iccp3m_cfg.update, 1, MPI_INT, 0, MPI_COMM_WORLD);
  
  /* broadcast the vectors element by element. This is slow
   * but safe and only performed at the beginning of each simulation*/
  for ( i = 0; i < iccp3m_cfg.last_ind_id +1; i++) {
    MPI_Bcast((double*)&iccp3m_cfg.areas[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast((double*)&iccp3m_cfg.ein[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast((double*)&iccp3m_cfg.nvectorx[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast((double*)&iccp3m_cfg.nvectory[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast((double*)&iccp3m_cfg.nvectorz[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast((double*)&iccp3m_cfg.extx[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast((double*)&iccp3m_cfg.exty[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast((double*)&iccp3m_cfg.extz[i], 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  }

  MPI_Bcast(&iccp3m_cfg.citeration, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(&iccp3m_cfg.set_flag, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

  printf("node %d: no iterations: %d\n", this_node, iccp3m_cfg.num_iteration);
  return 0 ;
    
}

int iccp3m_iteration() {
   double fdot,hold,hnew,hmax,del_eps,diff=0.0,difftemp=0.0,qold, ex, ey, ez, l_b;
   Cell *cell;
   int c,np;
   Particle *part;
   int i, j,id;
   char* errtxt;
   double globalmax;

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
                id = part[i].p.identity ;
                if( id <= iccp3m_cfg.last_ind_id) {
           /* the dielectric-related prefactor: */                     
                      del_eps = (iccp3m_cfg.ein[id]-iccp3m_cfg.eout)/(iccp3m_cfg.ein[id] + iccp3m_cfg.eout)/6.283185307;
           /* calculate the electric field at the certain position */
                      ex=part[i].f.f[0]/part[i].p.q;
                      ey=part[i].f.f[1]/part[i].p.q;
                      ez=part[i].f.f[2]/part[i].p.q;

          /* let's add the contribution coming from the external field */
                      ex += iccp3m_cfg.extx[id]; 
                      ey += iccp3m_cfg.exty[id]; 
                      ez += iccp3m_cfg.extz[id];
                      
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
                      qold=part[i].p.q;
          /* determine if it is higher than the previously highest charge density */            
                      if(hold>fabs(hmax))hmax=fabs(hold); 

                      hnew=(1.-iccp3m_cfg.relax)*hold + (iccp3m_cfg.relax)*del_eps*fdot/l_b;
                      difftemp=fabs( 2.*(hnew - hold)/(hold+hnew) ); /* relative variation: never use 
                                                                              an estimator which can be negative
                                                                              here */
//                      printf("(%d) difftemp = %f diff = %f\n",this_node,difftemp,diff);
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
       MPI_Allreduce(&diff, &globalmax, 1,MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

       if (globalmax < iccp3m_cfg.convergence) 
         break; 
       if ( diff > 1e89 ) /* Error happened */
         return iccp3m_cfg.citeration++;

  } /* iteration */

  on_particle_change();
  return iccp3m_cfg.citeration++;
}

void force_calc_iccp3m() {
/* The following ist mostly copied from forces.c */

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
	  if(dist2 <= max_range_non_bonded2) {
	    ONEPART_TRACE(if(p1[i].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d,%d %d,%d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,i,n,j,sqrt(dist2)));
	    ONEPART_TRACE(if(p2[j].p.identity==check_id) fprintf(stderr,"%d: OPT: Verlet Pair %d %d (Cells %d %d dist %f)\n",this_node,p1[i].p.identity,p2[j].p.identity,c,n,sqrt(dist2)));

	    add_pair_iccp3m(pl, &p1[i], &p2[j]);
	    /* calc non bonded interactions */ 
           if(!(p1[i].p.identity > iccp3m_cfg.last_ind_id && p2[j].p.identity >iccp3m_cfg.last_ind_id)) {
	       add_non_bonded_pair_force_iccp3m(&(p1[i]), &(p2[j]), vec21, sqrt(dist2), dist2);
                } 
/*           else {
             printf("Ever here \n");
             }*/
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

  npl   = local_icc->n;
  partl = local_icc->part;

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
    for (c = 0; c < me_do_ghosts_icc.n; c++) {
      npg   = me_do_ghosts_icc.cell[c]->n;
      partg = me_do_ghosts_icc.cell[c]->part;

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
  /* copied from forces.c */
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
#ifdef ELC_P3M
    case COULOMB_ELC_P3M:
       if (elc_params.dielectric_contrast_on) {
                errtxt = runtime_error(128);
                ERROR_SPRINTF(errtxt, "{ICCP3M conflicts with ELC dielectric constrast");
       }
       P3M_charge_assign();
       P3M_calc_kspace_forces(1,0);
       ELC_add_force(); 
    break;
#endif

#ifdef ELP3M
    case COULOMB_P3M:
         P3M_charge_assign();

#ifdef NPT
        if(integ_switch == INTEG_METHOD_NPT_ISO) exit(0);
#endif
        P3M_calc_kspace_forces_for_charges(1,0);
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

static int tclcommand_iccp3m_parse_params(Tcl_Interp *interp,int normal_args, char *string, int flag) {
  /* This function parses a vector give as C-String */
  /* It currently is not too elegant. But works. */
  int size,i,k=0;
  int scan_succes;
  static double *numbers=NULL;
  const char delimiters[] = " ";
  char *token,*cp;
  float temp;

  size= normal_args;
  numbers = malloc((size)*sizeof(double));

  cp = strdup(string);                /* Make writable copy.  */
  token = strtok (cp, delimiters);
  scan_succes = sscanf(token,"%f",&temp);
  if (scan_succes < 1) 
    return 1;

  numbers[0]=temp;
    
  for(i=1;i<size;i++) {
    token = strtok (NULL, delimiters);
    if(token == NULL)
      return 1;
    scan_succes = sscanf(token,"%lf",&(numbers[i]));
    if (scan_succes < 1 ) 
      return 1;
  }

  switch(flag) {
    case ICCP3M_AREA: 
      iccp3m_cfg.areas = (double*) realloc(iccp3m_cfg.areas, (size+1)*sizeof(double)); 
      for( i = 0 ; i < size ; i++ )  
        iccp3m_cfg.areas[i]=numbers[i];
      break;
    case ICCP3M_EPSILON:
      iccp3m_cfg.ein = (double*) realloc(iccp3m_cfg.ein,(size+1)*sizeof(double));
      for( i = 0 ; i < size; i++)  
        iccp3m_cfg.ein[i]=numbers[i];
    break;
    case ICCP3M_NORMAL:
      iccp3m_cfg.nvectorx = (double*) realloc(iccp3m_cfg.nvectorx,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      iccp3m_cfg.nvectory = (double*) realloc(iccp3m_cfg.nvectory,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      iccp3m_cfg.nvectorz = (double*) realloc(iccp3m_cfg.nvectorz,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      for(i=0;i<size;i++) {
        if( i%3 == 0 ) { iccp3m_cfg.nvectorx[k] = numbers[i]; } 
        if( i%3 == 1 ) { iccp3m_cfg.nvectory[k] = numbers[i]; }
        if( i%3 == 2 ) { iccp3m_cfg.nvectorz[k] = numbers[i];  k++; } 
       }
    break;

    case ICCP3M_EXTFIELD:
      iccp3m_cfg.extx = (double*) realloc(iccp3m_cfg.extx,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      iccp3m_cfg.exty = (double*) realloc(iccp3m_cfg.exty,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      iccp3m_cfg.extz = (double*) realloc(iccp3m_cfg.extz,sizeof(double)*(iccp3m_cfg.last_ind_id+1));
      for(i=0;i<size;i++) {
        if( i%3 == 0 ) { iccp3m_cfg.extx[k] = numbers[i]; } 
        if( i%3 == 1 ) { iccp3m_cfg.exty[k] = numbers[i]; }
        if( i%3 == 2 ) { iccp3m_cfg.extz[k] = numbers[i];  k++; } 
      }
    break;
  }

  free(numbers);
  return (0);
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

#endif


// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
/** \file statistics.c
    This is the place for analysis (so far...).
    Implementation of statistics.h
*/
#include <tcl.h>
#include <stdlib.h>
#include <string.h>
#include "utils.h"
#include "statistics.h"
#include "statistics_chain.h"
#include "statistics_molecule.h"
#include "statistics_cluster.h"
#include "statistics_fluid.h"
#include "energy.h"
#include "modes.h"
#include "pressure.h"
#include "communication.h"
#include "grid.h"
#include "parser.h"
#include "particle_data.h"
#include "interaction_data.h"
#include "domain_decomposition.h"
#include "verlet.h"
#include "lb.h"
#include "virtual_sites.h"
#include "initialize.h"

/** Previous particle configurations (needed for offline analysis and correlation analysis in \ref #analyze) */
double **configs = NULL; int n_configs = 0; int n_part_conf = 0;

/** Variables for measuring the compressibility from volume fluctuations.
    Will be used by \ref parse_Vkappa exclusively. */
typedef struct {
  /** sum of all the considered volumes resp. volumes squared so far */
  double Vk1; double Vk2; 
  /** amount of considered volumes so far */
  double avk;
} Vkappa_struct;
Vkappa_struct Vkappa = { 0.0, 0.0, 0.0 };

/****************************************************************************************
 *                                 helper functions
 ****************************************************************************************/

double min_distance2(double pos1[3], double pos2[3])
{
  double diff[3];
  get_mi_vector(diff, pos1, pos2);
  return sqrlen(diff);
}

/** Parses a reference point, i.e. either three doubles representing a
    position, or a particle id. */
static int get_reference_point(Tcl_Interp *interp, int *argc, char ***argv,
			       double pos[3], int *pid)
{
  *pid = -1;

  if (*argc >= 3 &&
      (Tcl_GetDouble(interp, (*argv)[0], &pos[0]) != TCL_ERROR) &&
      (Tcl_GetDouble(interp, (*argv)[1], &pos[1]) != TCL_ERROR) &&
      (Tcl_GetDouble(interp, (*argv)[2], &pos[2]) != TCL_ERROR)) {
    /* Found three doubles representing a position. */
    (*argc) -= 3;
    (*argv) += 3;
    return TCL_OK;
  }
  /* else */
  if (*argc >= 1 && 
      (Tcl_GetInt(interp, (*argv)[0], pid) != TCL_ERROR)) {
    /* Found a single integer representing a particle id. */
    Particle ref;
    
    if (get_particle_data(*pid, &ref) != TCL_OK) {
      Tcl_AppendResult(interp, "reference particle does not exist", (char *)NULL);
      return TCL_ERROR;
    }
    pos[0] = ref.r.p[0];
    pos[1] = ref.r.p[1];
    pos[2] = ref.r.p[2];
    
    (*argc)--;
    (*argv)++;

    free_particle(&ref);
    return TCL_OK;
  }
  /* else */
  return TCL_ERROR;
}

/****************************************************************************************
 *                                 basic observables calculation
 ****************************************************************************************/

double mindist(IntList *set1, IntList *set2)
{
  double mindist, pt[3];
  int i, j, in_set;

  mindist = SQR(box_l[0] + box_l[1] + box_l[2]);

  updatePartCfg(WITHOUT_BONDS);
  for (j=0; j<n_total_particles-1; j++) {
    pt[0] = partCfg[j].r.p[0];
    pt[1] = partCfg[j].r.p[1];
    pt[2] = partCfg[j].r.p[2];
    /* check which sets particle j belongs to
       bit 0: set1, bit1: set2
    */
    in_set = 0;
    if (!set1 || intlist_contains(set1, partCfg[j].p.type))
      in_set = 1;
    if (!set2 || intlist_contains(set2, partCfg[j].p.type))
      in_set |= 2;
    if (in_set == 0)
      continue;

    for (i=j+1; i<n_total_particles; i++)
      /* accept a pair if particle j is in set1 and particle i in set2 or vice versa. */
      if (((in_set & 1) && (!set2 || intlist_contains(set2, partCfg[i].p.type))) ||
	  ((in_set & 2) && (!set1 || intlist_contains(set1, partCfg[i].p.type))))
	mindist = dmin(mindist, min_distance2(pt, partCfg[i].r.p));
  }
  mindist = sqrt(mindist);
  return mindist;
}

void merge_aggregate_lists(int *head_list, int *agg_id_list, int p1molid, int p2molid, int *link_list)
{
    int target1, target2, head_p1;
    /* merge list containing p2molid into list containing p1molid*/
    target1=head_list[agg_id_list[p2molid]];
    head_list[agg_id_list[p2molid]]=-2;
    head_p1=head_list[agg_id_list[p1molid]];
    head_list[agg_id_list[p1molid]]=target1;
    agg_id_list[target1]=agg_id_list[p1molid];
    target2=link_list[target1];
    while(target2 != -1) {
	target1=target2;
	target2=link_list[target1];
	agg_id_list[target1]=agg_id_list[p1molid];
    }
    agg_id_list[target1]=agg_id_list[p1molid];
    link_list[target1]=head_p1;
}

int aggregation(double dist_criteria2, int min_contact, int s_mol_id, int f_mol_id, int *head_list, int *link_list, int *agg_id_list, int *agg_num, int *agg_size, int *agg_max, int *agg_min, int *agg_avg, int *agg_std, int charge)
{
  int c, np, n, i;
  Particle *p1, *p2, **pairs;
  double dist2;
  int target1;
  int p1molid, p2molid;
  int *contact_num, ind;

  if (min_contact > 1) {
    contact_num = (int *) malloc(n_molecules*n_molecules *sizeof(int));
    for (i = 0; i < n_molecules *n_molecules; i++) contact_num[i]=0;
  } else {
    contact_num = (int *) 0; /* Just to keep the compiler happy */
  }

  on_observable_calc();
  build_verlet_lists();

  for (i = s_mol_id; i <= f_mol_id; i++) {
    head_list[i]=i;
    link_list[i]=-1;
    agg_id_list[i]=i;
    agg_size[i]=0;
  }
  
  /* Loop local cells */
  for (c = 0; c < local_cells.n; c++) {
    /* Loop cell neighbors */
    for (n = 0; n < dd.cell_inter[c].n_neighbors; n++) {
      pairs = dd.cell_inter[c].nList[n].vList.pair;
      np    = dd.cell_inter[c].nList[n].vList.n;
      /* verlet list loop */
      for(i=0; i<2*np; i+=2) {
	p1 = pairs[i];                    /* pointer to particle 1 */
	p2 = pairs[i+1];                  /* pointer to particle 2 */
	p1molid = p1->p.mol_id;
	p2molid = p2->p.mol_id;
	if (((p1molid <= f_mol_id) && (p1molid >= s_mol_id)) && ((p2molid <= f_mol_id) && (p2molid >= s_mol_id))) {
	  if (agg_id_list[p1molid] != agg_id_list[p2molid]) {
	    dist2 = min_distance2(p1->r.p, p2->r.p);

#ifdef ELECTROSTATICS
	    if (charge && (p1->p.q * p2->p.q >= 0)) {continue;}
#endif
	    if (dist2 < dist_criteria2) {
	      if ( p1molid > p2molid ) { ind=p1molid*n_molecules + p2molid;} 
	      else { ind=p2molid*n_molecules +p1molid;}
	      if (min_contact > 1) {
		contact_num[ind] ++;
		if (contact_num[ind] >= min_contact) {
		    merge_aggregate_lists( head_list, agg_id_list, p1molid, p2molid, link_list);				    
		}
	      } else {
		  merge_aggregate_lists( head_list, agg_id_list, p1molid, p2molid, link_list);				    
	      }
	    }
	  }
	}
      }
    }
  }
  
  /* count number of aggregates 
     find aggregate size
     find max and find min size, and std */
  for (i = s_mol_id ; i <= f_mol_id ; i++) {
    if (head_list[i] != -2) {
      (*agg_num)++;
      agg_size[*agg_num -1]++;
      target1= head_list[i];
      while( link_list[target1] != -1) {
	target1= link_list[target1];
	agg_size[*agg_num -1]++;
      }
    }
  }
  
  for (i = 0 ; i < *agg_num; i++) {
    *agg_avg += agg_size[i];
    *agg_std += agg_size[i] * agg_size[i];
    if (*agg_min > agg_size[i]) { *agg_min = agg_size[i]; }
    if (*agg_max < agg_size[i]) { *agg_max = agg_size[i]; }
  }
  
  return 0;
}

/** Calculate momentum of all particles in the local domain
 * @param result Result for this processor (Output)
 */
void predict_momentum_particles(double *result)
{
  Cell *cell;
  Particle *p;
  int i, c, np;

  double momentum[3] = { 0.0, 0.0, 0.0 };

  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    np = cell->n;
    p  = cell->part;

    for(i=0; i < np; i++) {
         momentum[0] +=  p[i].m.v[0] + p[i].f.f[0];
         momentum[1] +=  p[i].m.v[1] + p[i].f.f[1];
         momentum[2] +=  p[i].m.v[2] + p[i].f.f[2];
    }
  }

  momentum[0] /= time_step;
  momentum[1] /= time_step;
  momentum[2] /= time_step;

  MPI_Reduce(momentum, result, 3, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
}

/** Calculate total momentum of the system (particles & LB fluid)
 * @param momentum Rsult for this processor (Output)
 */
void momentum_calc(double *momentum) 
{
    double momentum_fluid[3] = { 0., 0., 0. };
    double momentum_particles[3] = { 0., 0., 0. };

    mpi_gather_stats(4, momentum_particles, NULL, NULL, NULL);
    mpi_gather_stats(6, momentum_fluid, NULL, NULL, NULL);

    momentum[0] = momentum_fluid[0] + momentum_particles[0];
    momentum[1] = momentum_fluid[1] + momentum_particles[1];
    momentum[2] = momentum_fluid[2] + momentum_particles[2];

}

void centermass(int type, double *com)
{
  int i, j;
  double M = 0.0;
  com[0]=com[1]=com[2]=0.;
   	
  updatePartCfg(WITHOUT_BONDS);
  for (j=0; j<n_total_particles; j++) {
    if ((partCfg[j].p.type == type) || (type == -1)) {
      for (i=0; i<3; i++) {
      	com[i] += partCfg[j].r.p[i]*PMASS(partCfg[j]);
      }
      M += PMASS(partCfg[j]);
    }
  }
  
  for (i=0; i<3; i++) {
    com[i] /= M;
  }
  return;
}

void centermass_vel(int type, double *com)
{
  /*center of mass velocity scaled with time_step*/
  int i, j;
  int count = 0;
  com[0]=com[1]=com[2]=0.;

  updatePartCfg(WITHOUT_BONDS);
  for (j=0; j<n_total_particles; j++) {
    if (type == partCfg[j].p.type) {
      for (i=0; i<3; i++) {
      	com[i] += partCfg[j].m.v[i];
      }
      count++;
    }
  }

  for (i=0; i<3; i++) {
    com[i] /= count;
  }
  return;
}

void angularmomentum(int type, double *com)
{
  int i, j;
  double tmp[3];
  double pre_factor;
  com[0]=com[1]=com[2]=0.;

  updatePartCfg(WITHOUT_BONDS);
  for (j=0; j<n_total_particles; j++) 
  {
    if (type == partCfg[j].p.type) 
    {
      vector_product(partCfg[j].r.p,partCfg[j].m.v,tmp);
      pre_factor=PMASS(partCfg[j]);
      for (i=0; i<3; i++) {
        com[i] += tmp[i]*pre_factor;
      }
    }
  }
  return;
}

void  momentofinertiamatrix(int type, double *MofImatrix)
{
  int i,j,count;
  double p1[3],com[3],massi;

  count=0;
  updatePartCfg(WITHOUT_BONDS);
  for(i=0;i<9;i++) MofImatrix[i]=0.;
  centermass(type, com);
  for (j=0; j<n_total_particles; j++) {
    if (type == partCfg[j].p.type) {
      count ++;
      for (i=0; i<3; i++) {
      	p1[i] = partCfg[j].r.p[i] - com[i];
      }
      massi= PMASS(partCfg[j]);
      MofImatrix[0] += massi * (p1[1] * p1[1] + p1[2] * p1[2]) ; 
      MofImatrix[4] += massi * (p1[0] * p1[0] + p1[2] * p1[2]);
      MofImatrix[8] += massi * (p1[0] * p1[0] + p1[1] * p1[1]);
      MofImatrix[1] -= massi * (p1[0] * p1[1]);
      MofImatrix[2] -= massi * (p1[0] * p1[2]); 
      MofImatrix[5] -= massi * (p1[1] * p1[2]);
    }
  }
  /* use symmetry */
  MofImatrix[3] = MofImatrix[1]; 
  MofImatrix[6] = MofImatrix[2]; 
  MofImatrix[7] = MofImatrix[5];
  return;
}


void nbhood(double pt[3], double r, IntList *il, int planedims[3] )
{
  double d[3],dsize;
  int i,j;
  double r2;

  r2 = r*r;

  init_intlist(il);
 
  updatePartCfg(WITHOUT_BONDS);

  for (i = 0; i<n_total_particles; i++) {
    if ( (planedims[0] + planedims[1] + planedims[2]) == 3 ) {
      get_mi_vector(d, pt, partCfg[i].r.p);
    } else {
      /* Calculate the in plane distance */
      dsize = 0.0;
      for ( j= 0 ; j < 3 ; j++ ) {
	d[j] = planedims[j]*(partCfg[i].r.p[j]-pt[j]);
      }
    }

    if (sqrlen(d) < r2) {
      realloc_intlist(il, il->n + 1);
      il->e[il->n] = partCfg[i].p.identity;
      il->n++;
    }
  }
}

double distto(double p[3], int pid)
{
  int i;
  double d[3];
  double mindist;

  /* larger than possible */
  mindist=SQR(box_l[0] + box_l[1] + box_l[2]);
  for (i=0; i<n_total_particles; i++) {
    if (pid != partCfg[i].p.identity) {
      get_mi_vector(d, p, partCfg[i].r.p);
      mindist = dmin(mindist, sqrlen(d));
    }
  }
  return sqrt(mindist);
}

void calc_cell_gpb(double xi_m, double Rc, double ro, double gacc, int maxtry, double *result) {
  double LOG,xi_min, RM, gamma,g1,g2,gmid=0,dg,ig, f,fmid, rtb;
  int i;
  LOG    = log(Rc/ro);
  xi_min = LOG/(1+LOG);
  if(maxtry < 1) maxtry = 1;

  /* determine which of the regimes we are in: */
  if(xi_m > 1) {
    ig = 1.0;
    g1 = PI / LOG;
    g2 = PI / (LOG + xi_m/(xi_m-1.0)); }
  else if(xi_m == 1) {
    ig = 1.0;
    g1 = (PI/2.0) / LOG;
    g2 = (PI/2.0) / (LOG + 1.0); }
  else if(xi_m == xi_min) {
    ig = 1.0;
    g1 = g2 = 0.0; }
  else if(xi_m > xi_min) {
    ig = 1.0;
    g1 = (PI/2.0) / LOG;
    g2 = sqrt(3.0*(LOG-xi_m/(1.0-xi_m))/(1-pow((1.0-xi_m),-3.0))); }
  else if(xi_m > 0.0) {
    ig = -1.0;
    g1 = 1-xi_m;
    g2 = xi_m*(6.0-(3.0-xi_m)*xi_m)/(3.0*LOG); }
  else if(xi_m == 0.0) {
    ig = -1.0;
    g1 = g2 = 1-xi_m; }
  else {
    result[2]=-5.0; return;
  }

  /* decide which method to use (if any): */
  if(xi_m == xi_min) {
    gamma = 0.0;
    RM    = 0.0; }
  else if(xi_m == 0.0) {
    gamma = 1-xi_m;
    RM    = -1.0; }
  else if(ig == 1.0) {
    /* determine gamma via a bisection-search: */
    f    = atan(1.0/g1) + atan( (xi_m-1.0) / g1 ) - g1 * LOG;
    fmid = atan(1.0/g2) + atan( (xi_m-1.0) / g2 ) - g2 * LOG;
    if (f*fmid >= 0.0) {
      /* failed to bracket function value with intial guess - abort: */
      result[0]=f; result[1]=fmid; result[2]=-3.0; return;  }

    /* orient search such that the positive part of the function lies to the right of the zero */
    rtb = f < 0.0 ? (dg=g2-g1,g1) : (dg=g1-g2,g2);
    for (i = 1; i <= maxtry; i++) {
      gmid = rtb + (dg *= 0.5);
      fmid = atan(1.0/gmid) + atan((xi_m-1.0)/gmid) - gmid*LOG;
      if (fmid <= 0.0) rtb = gmid;
      if (fabs(dg) < gacc || fmid == 0.0) break;
    }

    if (fabs(dg) > gacc) {
      /* too many iterations without success - abort: */
      result[0]=gmid; result[1]=dg; result[2]=-2.0; return;  }

    /* So, these are the values for gamma and Manning-radius: */
    gamma = gmid;
    RM    = Rc * exp( -(1.0/gamma) * atan(1.0/gamma) ); }
  else if(ig == -1.0) {
    /* determine -i*gamma: */
    f = -1.0*(atanh(g2) + atanh(g2/(xi_m-1))) - g2*LOG;

    /* modified orient search, this time starting from the upper bound backwards: */
    if (f < 0.0) {
      rtb = g1;  dg = g1-g2; }
    else {
      fprintf(stderr,"WARNING: Lower boundary is actually larger than l.h.s, flipping!\n");
      rtb = g1;  dg = g1;    }
    for (i = 1; i <= maxtry; i++) {
      gmid = rtb - (dg *= 0.5);
      fmid = -1.0*(atanh(gmid) + atanh(gmid/(xi_m-1))) - gmid*LOG;
      if (fmid >= 0.0) rtb = gmid;
      if (fabs(dg) < gacc || fmid == 0.0) break;
    }
    
    if (fabs(dg) > gacc) {
      /* too many iterations without success - abort: */
      result[0]=gmid; result[1]=dg; result[2]=-2.0; return;  }

    /* So, these are the values for -i*gamma and Manning-radius: */
    gamma = gmid;
    RM    = Rc * exp( atan(1.0/gamma)/gamma ); }
  else {
    result[2]=-5.0; return;
  }

  result[0]=gamma;
  result[1]=RM;
  result[2]=ig;
  return;
}

void calc_part_distribution(int *p1_types, int n_p1, int *p2_types, int n_p2, 
			    double r_min, double r_max, int r_bins, int log_flag, 
			    double *low, double *dist)
{
  int i,j,t1,t2,ind,cnt=0;
  double inv_bin_width=0.0;
  double min_dist,min_dist2=0.0,start_dist2,act_dist2;

  start_dist2 = SQR(box_l[0] + box_l[1] + box_l[2]);
  /* bin preparation */
  *low = 0.0;
  for(i=0;i<r_bins;i++) dist[i] = 0.0;
  if(log_flag == 1) inv_bin_width = (double)r_bins/(log(r_max)-log(r_min));
  else              inv_bin_width = (double)r_bins / (r_max-r_min);

  /* particle loop: p1_types*/
  for(i=0; i<n_total_particles; i++) {
    for(t1=0; t1<n_p1; t1++) {
      if(partCfg[i].p.type == p1_types[t1]) {
	min_dist2 = start_dist2;
	/* particle loop: p2_types*/
	for(j=0; j<n_total_particles; j++) {
	  if(j != i) {
	    for(t2=0; t2<n_p2; t2++) {
	      if(partCfg[j].p.type == p2_types[t2]) {
		act_dist2 =  min_distance2(partCfg[i].r.p, partCfg[j].r.p);
		if(act_dist2 < min_dist2) { min_dist2 = act_dist2; }
	      }
	    }
	  }
	}
	min_dist = sqrt(min_dist2);
	if(min_dist <= r_max) {
	  if(min_dist >= r_min) {
	    /* calculate bin index */
	    if(log_flag == 1) ind = (int) ((log(min_dist) - log(r_min))*inv_bin_width);
	    else              ind = (int) ((min_dist - r_min)*inv_bin_width);
	    if(ind >= 0 && ind < r_bins) {
	      dist[ind] += 1.0;
	    }
	  }
	  else {
	    *low += 1.0;
	  }
	}
	cnt++;    
      }
    }
  }
  
  /* normalization */
  *low /= (double)cnt;
  for(i=0;i<r_bins;i++) dist[i] /= (double)cnt;
}

void calc_vel_distr(Tcl_Interp *interp, int type,int bins,double given_max)
{
   int i,j,p_count,dist_count,ind;
   double min,max,bin_width,inv_bin_width,com[3],vel;
   long distribution[bins];
   char buffer[2*TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE+256];

   max=given_max*time_step;
   min=-given_max*time_step;
   p_count=0;
   for(i=0; i<bins; i++) {distribution[i]=0;}

   centermass_vel(type,com);

   for (i=0;i<n_total_particles;i++)
   {
      if (partCfg[i].p.type==type)
      {
        p_count++;
        for (j=0;j<3;j++)
        {
           vel=partCfg[i].m.v[j] - com[j];
           if (min > vel ){min = vel;}
           if (max < vel ){max = vel;}
        }
      }
   }

   if (p_count==0) {return;}

   if ( (-min) > max ) {
      max = -min;
   }
   else{
      min = -max;
   }

   bin_width     = (max-min) / (double)bins;
   inv_bin_width = 1.0 / bin_width;
   dist_count=0;
   for (i=0;i<n_total_particles;i++)
   {
     if (partCfg[i].p.type==type)
     {
        for (j=0;j<3;j++)
        {
          vel= partCfg[i].m.v[j] - com[j];
          ind = (int) ( (vel - min)*inv_bin_width );
          distribution[ind]++;
          dist_count++;
        }
     }
   }
   
   vel=min + bin_width/2.0;
   Tcl_AppendResult(interp, " {\n", (char *)NULL);
   for(i=0; i<bins; i++) {
      sprintf(buffer,"%f %f",vel/time_step,distribution[i]/(double)dist_count);
      Tcl_AppendResult(interp, "{ ",buffer," }\n", (char *)NULL);
      vel += bin_width;
   }
   Tcl_AppendResult(interp, "}\n", (char *)NULL);
}

void calc_rdf(int *p1_types, int n_p1, int *p2_types, int n_p2, 
	      double r_min, double r_max, int r_bins, double *rdf)
{
  int i,j,t1,t2,ind,cnt=0;
  int mixed_flag=0,start;
  double inv_bin_width=0.0,bin_width=0.0, dist;
  double volume, bin_volume, r_in, r_out;

  if(n_p1 == n_p2) {
    for(i=0;i<n_p1;i++) 
      if( p1_types[i] != p2_types[i] ) mixed_flag=1;
  }
  else mixed_flag=1;

  bin_width     = (r_max-r_min) / (double)r_bins;
  inv_bin_width = 1.0 / bin_width;
  for(i=0;i<r_bins;i++) rdf[i] = 0.0;
  /* particle loop: p1_types*/
  for(i=0; i<n_total_particles; i++) {
    for(t1=0; t1<n_p1; t1++) {
      if(partCfg[i].p.type == p1_types[t1]) {
	/* distinguish mixed and identical rdf's */
	if(mixed_flag == 1) start = 0;
	else                start = (i+1);
	/* particle loop: p2_types*/
	for(j=start; j<n_total_particles; j++) {
	  for(t2=0; t2<n_p2; t2++) {
	    if(partCfg[j].p.type == p2_types[t2]) {
	      dist = min_distance(partCfg[i].r.p, partCfg[j].r.p);
	      if(dist > r_min && dist < r_max) {
		ind = (int) ( (dist - r_min)*inv_bin_width );
		rdf[ind]++;
	      }
	      cnt++;
	    }
	  }
	}
      }
    }
  }

  /* normalization */
  volume = box_l[0]*box_l[1]*box_l[2];
  for(i=0; i<r_bins; i++) {
    r_in       = i*bin_width + r_min; 
    r_out      = r_in + bin_width;
    bin_volume = (4.0/3.0) * PI * ((r_out*r_out*r_out) - (r_in*r_in*r_in));
    rdf[i] *= volume / (bin_volume * cnt);
  }
}

void calc_rdf_av(int *p1_types, int n_p1, int *p2_types, int n_p2,
		 double r_min, double r_max, int r_bins, double *rdf, int n_conf)
{
  int i,j,k,l,t1,t2,ind,cnt=0,cnt_conf=1;
  int mixed_flag=0,start;
  double inv_bin_width=0.0,bin_width=0.0, dist;
  double volume, bin_volume, r_in, r_out;
  double *rdf_tmp, p1[3],p2[3];

  rdf_tmp = malloc(r_bins*sizeof(double));

  if(n_p1 == n_p2) {
    for(i=0;i<n_p1;i++)
      if( p1_types[i] != p2_types[i] ) mixed_flag=1;
  }
  else mixed_flag=1;
  bin_width     = (r_max-r_min) / (double)r_bins;
  inv_bin_width = 1.0 / bin_width;
  volume = box_l[0]*box_l[1]*box_l[2];
  for(l=0;l<r_bins;l++) rdf_tmp[l]=rdf[l] = 0.0;

  while(cnt_conf<=n_conf) {
    for(l=0;l<r_bins;l++) rdf_tmp[l]=0.0;
    cnt=0;
    k=n_configs-cnt_conf;
    for(i=0; i<n_total_particles; i++) {
      for(t1=0; t1<n_p1; t1++) {
	if(partCfg[i].p.type == p1_types[t1]) {
	  // distinguish mixed and identical rdf's
	  if(mixed_flag == 1) start = 0;
	  else                start = (i+1);
	  //particle loop: p2_types
	  for(j=start; j<n_total_particles; j++) {
	    for(t2=0; t2<n_p2; t2++) {
	      if(partCfg[j].p.type == p2_types[t2]) {
		p1[0]=configs[k][3*i  ];p1[1]=configs[k][3*i+1];p1[2]=configs[k][3*i+2];
		p2[0]=configs[k][3*j  ];p2[1]=configs[k][3*j+1];p2[2]=configs[k][3*j+2];
		dist =min_distance(p1, p2);
		if(dist > r_min && dist < r_max) {
		  ind = (int) ( (dist - r_min)*inv_bin_width );
		  rdf_tmp[ind]++;
		}
		cnt++;
	      }
	    }
	  }
	}
      }
    }
    // normalization
  
    for(i=0; i<r_bins; i++) {
      r_in       = i*bin_width + r_min;
      r_out      = r_in + bin_width;
      bin_volume = (4.0/3.0) * PI * ((r_out*r_out*r_out) - (r_in*r_in*r_in));
      rdf[i] += rdf_tmp[i]*volume / (bin_volume * cnt);
    }

    cnt_conf++;
  } //cnt_conf loop
  for(i=0; i<r_bins; i++) {
    rdf[i] /= (cnt_conf-1);
  }
  free(rdf_tmp);

}

void calc_rdf_intermol_av(int *p1_types, int n_p1, int *p2_types, int n_p2,
			  double r_min, double r_max, int r_bins, double *rdf, int n_conf)
{
  int i,j,k,l,t1,t2,ind,cnt=0,cnt_conf=1;
  int mixed_flag=0,start;
  double inv_bin_width=0.0,bin_width=0.0, dist;
  double volume, bin_volume, r_in, r_out;
  double *rdf_tmp, p1[3],p2[3];

  rdf_tmp = malloc(r_bins*sizeof(double));

  if(n_p1 == n_p2) {
    for(i=0;i<n_p1;i++)
      if( p1_types[i] != p2_types[i] ) mixed_flag=1;
  }
  else mixed_flag=1;
  bin_width     = (r_max-r_min) / (double)r_bins;
  inv_bin_width = 1.0 / bin_width;
  volume = box_l[0]*box_l[1]*box_l[2];
  for(l=0;l<r_bins;l++) rdf_tmp[l]=rdf[l] = 0.0;

  while(cnt_conf<=n_conf) {
    for(l=0;l<r_bins;l++) rdf_tmp[l]=0.0;
    cnt=0;
    k=n_configs-cnt_conf;
    for(i=0; i<n_total_particles; i++) {
      for(t1=0; t1<n_p1; t1++) {
	if(partCfg[i].p.type == p1_types[t1]) {
	  // distinguish mixed and identical rdf's
	  if(mixed_flag == 1) start = 0;
	  else                start = (i+1);
	  //particle loop: p2_types
	  for(j=start; j<n_total_particles; j++) {
	    for(t2=0; t2<n_p2; t2++) {
	      if(partCfg[j].p.type == p2_types[t2]) {
		/*see if particles i and j belong to different molecules*/
		if(partCfg[i].p.mol_id!=partCfg[j].p.mol_id) {
		  p1[0]=configs[k][3*i  ];p1[1]=configs[k][3*i+1];p1[2]=configs[k][3*i+2];
		  p2[0]=configs[k][3*j  ];p2[1]=configs[k][3*j+1];p2[2]=configs[k][3*j+2];
		  dist =min_distance(p1, p2);
		  if(dist > r_min && dist < r_max) {
		    ind = (int) ( (dist - r_min)*inv_bin_width );
		    rdf_tmp[ind]++;
		  }
		  cnt++;
		}
	      }
	    }
	  }
	}
      }
    }
    // normalization

    for(i=0; i<r_bins; i++) {
      r_in       = i*bin_width + r_min;
      r_out      = r_in + bin_width;
      bin_volume = (4.0/3.0) * PI * ((r_out*r_out*r_out) - (r_in*r_in*r_in));
      rdf[i] += rdf_tmp[i]*volume / (bin_volume * cnt);
    }

    cnt_conf++;
  } //cnt_conf loop
  for(i=0; i<r_bins; i++) {
    rdf[i] /= (cnt_conf-1);
  }
  free(rdf_tmp);

}

/*addes this line*/
void calc_rdf_adress(int *p1_types, int n_p1, int *p2_types, int n_p2,
			   double x_min, double x_max, double r_min, double r_max, int r_bins, double *rdf, int n_conf)
{
  int i,j,k,l,t1,t2,ind,cnt=0,cnt_conf=1;
  int mixed_flag=0,start;
  double inv_bin_width=0.0,bin_width=0.0, dist;
  double volume, bin_volume, r_in, r_out;
  double *rdf_tmp, p1[3],p2[3];

  rdf_tmp = malloc(r_bins*sizeof(double));

  if(n_p1 == n_p2) {
    for(i=0;i<n_p1;i++)
      if( p1_types[i] != p2_types[i] ) mixed_flag=1;
  }
  else mixed_flag=1;
  bin_width     = (r_max-r_min) / (double)r_bins;
  inv_bin_width = 1.0 / bin_width;
  volume = (x_max-x_min)*box_l[1]*box_l[2];
  //volume = box_l[0]*box_l[1]*box_l[2];
  for(l=0;l<r_bins;l++) rdf_tmp[l]=rdf[l] = 0.0;

  while(cnt_conf<=n_conf) {
    for(l=0;l<r_bins;l++) rdf_tmp[l]=0.0;
    cnt=0;
    k=n_configs-cnt_conf;
    for(i=0; i<n_total_particles; i++) {
      for(t1=0; t1<n_p1; t1++) {
	if(partCfg[i].p.type == p1_types[t1]) {
          //accepts particles i's contained in explicit region 
          if(configs[k][3*i]<x_max && configs[k][3*i]>=x_min) {
	  // distinguish mixed and identical rdf's
	  if(mixed_flag == 1) start = 0;
	  else                start = (i+1);
	  //particle loop: p2_types
	  for(j=start; j<n_total_particles; j++) {
	    for(t2=0; t2<n_p2; t2++) {
	      if(partCfg[j].p.type == p2_types[t2]) {
                //accepts particles j's contained in explicit region 
                if(configs[k][3*j]<x_max && configs[k][3*j]>=x_min) {
		/*see if particles i and j belong to different molecules*/
		if(partCfg[i].p.mol_id!=partCfg[j].p.mol_id) {
		  p1[0]=configs[k][3*i  ];p1[1]=configs[k][3*i+1];p1[2]=configs[k][3*i+2];
		  p2[0]=configs[k][3*j  ];p2[1]=configs[k][3*j+1];p2[2]=configs[k][3*j+2];
		  dist =min_distance(p1, p2);
		  if(dist > r_min && dist < r_max) {
		    ind = (int) ( (dist - r_min)*inv_bin_width );
		    rdf_tmp[ind]++;
		  }
		  cnt++;
		}
	      }
             }
	    }
	  }
	}
        }
      }
    }
    // normalization

    for(i=0; i<r_bins; i++) {
      r_in       = i*bin_width + r_min;
      r_out      = r_in + bin_width;
      bin_volume = (4.0/3.0) * PI * ((r_out*r_out*r_out) - (r_in*r_in*r_in));
      rdf[i] += rdf_tmp[i]*volume / (bin_volume * cnt);
    }

    cnt_conf++;
  } //cnt_conf loop
  for(i=0; i<r_bins; i++) {
    rdf[i] /= (cnt_conf-1);
  }
  free(rdf_tmp);

}
/*Up to here*/

void calc_structurefactor(int type, int order, double **_ff) {
  int i, j, k, n, qi, p, order2;
  double qr, twoPI_L, C_sum, S_sum, *ff=NULL;
  
  order2 = order*order;
  *_ff = ff = realloc(ff,2*order2*sizeof(double));
  twoPI_L = 2*PI/box_l[0];
  
  if ((type < 0) || (type > n_particle_types)) { fprintf(stderr,"WARNING: Type %i does not exist!",type); fflush(NULL); errexit(); }
  else if (order < 1) { fprintf(stderr,"WARNING: parameter \"order\" has to be a whole positive number"); fflush(NULL); errexit(); }
  else {
    for(qi=0; qi<2*order2; qi++) {
      ff[qi] = 0.0;
    }
    for(i=0; i<=order; i++) {
      for(j=-order; j<=order; j++) {
        for(k=-order; k<=order; k++) {
	  n = i*i + j*j + k*k;
	  if ((n<=order2) && (n>=1)) {
	    C_sum = S_sum = 0.0;
	    for(p=0; p<n_total_particles; p++) {
	      if (partCfg[p].p.type == type) {
		qr = twoPI_L * ( i*partCfg[p].r.p[0] + j*partCfg[p].r.p[1] + k*partCfg[p].r.p[2] );
		C_sum+= cos(qr);
		S_sum+= sin(qr);
	      }
	    }
	    ff[2*n-2]+= C_sum*C_sum + S_sum*S_sum;
	    ff[2*n-1]++;
	  }
	}
      }
    }
    n = 0;
    for(p=0; p<n_total_particles; p++) {
      if (partCfg[p].p.type == type) n++;
    }
    for(qi=0; qi<order2; qi++) 
      if (ff[2*qi+1]!=0) ff[2*qi]/= n*ff[2*qi+1];
  }
}

//calculates average density profile in dir direction over last n_conf configurations
void density_profile_av(int n_conf, int n_bin, double density, int dir, double *rho_ave, int type)
{
  int i,j,k,m,n;
  double r;
  double r_bin;
  double pos[3];
  int  image_box[3];
  
  //calculation over last n_conf configurations  
  
  //bin width
  r_bin = box_l[dir]/(double)(n_bin);
  
  for (i=0; i<n_bin; i++)
    rho_ave[i]=0;
  
  k=n_configs-n_conf;
  
  while(k<n_configs) {
    r = 0;
    j = 0;
    while (r < box_l[dir]) { 
      n = 0;
      for(i=0; i<n_total_particles; i++) {
	//com particles
	if(partCfg[i].p.type == type) {
	  for(m=0; m<3; m++) {
	    pos[m] = configs[k][3*i+m];
	    image_box[m] = 0;
	  }
	  fold_coordinate(pos, image_box, dir);
	  if (pos[dir] <= r+r_bin && pos[dir] > r)
	    n++;
	}
      }
      
      rho_ave[j] += (double)(n)/(box_l[1]*box_l[2]*r_bin)/density;
      j++;
      r += r_bin;
    }     
    k++;
  } //k loop
  
  // normalization
  for (i=0; i<n_bin; i++)
    rho_ave[i]/=n_conf;
}

void calc_diffusion_profile(int dir, double xmin, double xmax, int nbins, int n_part, int n_conf, int time, int type, double *bins) 
{
  int i,t, count,index;
  double tcount=0;
  double xpos;
  double tpos[3];
  int img_box[3] = {0,0,0};
  //double delta_x = (box_l[0])/((double) nbins);
  
  /* create and initialize the array of bins */
  
  // double *bins;
  
  int *label;
  label = malloc(n_part*sizeof(int));
  
  /* calculation over last n_conf configurations */
  t=n_configs-n_conf;
  
  while (t<n_configs-time) {
    /* check initial condition */
    count = 0;
    
    for (i=0;i<n_total_particles;i++) {
      if(partCfg[i].p.type == type) {
	tpos[0] = configs[t][3*i];
	tpos[1] = configs[t][3*i+1];
	tpos[2] = configs[t][3*i+2];
	fold_coordinate(tpos, img_box,dir);
	xpos = tpos[dir];
	if(xpos > xmin && xpos < xmax) {
	  label[count] = i;
	}
	else label[count] = -1;
	count ++;
      }
    }
    
    /* check at time 'time' */
    for (i=0;i<n_part;i++) {
      if (label[i]>0) {
	tpos[0] = configs[t+time][3*label[i]];
	tpos[1] = configs[t+time][3*label[i]+1];
	tpos[2] = configs[t+time][3*label[i]+2];
	fold_coordinate(tpos, img_box,dir);
	xpos = tpos[dir];
	
	index = (int)(xpos/box_l[dir]*nbins);
	bins[index]++;
      }
    }
    t++;
    tcount++;
  }
  
  /* normalization */
  for (i=0;i<nbins;i++) {
    bins[i]=bins[i]/(tcount);
  }
  free(label);
}


int calc_radial_density_map (int xbins,int ybins,int thetabins,double xrange,double yrange, double axis[3], double center[3], IntList *beadids, DoubleList *density_map, DoubleList *density_profile) {
  int i,j,t;
  int pi,bi;
  int nbeadtypes;
  int beadcount;
  double vectprod[3];
  double pvector[3];
  double xdist,ydist,rdist,xav,yav,theta;
  double xbinwidth,ybinwidth,binvolume;
  double thetabinwidth;
  double *thetaradii;
  int *thetacounts;
  int xindex,yindex,tindex;
  xbinwidth = xrange/(double)(xbins);
  ybinwidth = yrange/(double)(ybins);

  nbeadtypes = beadids->n;
  /* Update particles */
  updatePartCfg(WITHOUT_BONDS);

  /*Make sure particles are folded  */
  for (i = 0 ; i < n_total_particles ; i++) {
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,0);
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,1);
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,2);
  }

  beadcount = 0;
  xav = 0.0;
  yav = 0.0;
  for ( pi = 0 ; pi < n_total_particles ; pi++ ) {
    for ( bi = 0 ; bi < nbeadtypes ; bi++ ) {
      if ( beadids->e[bi] == partCfg[pi].p.type ) {


	/* Find the vector from the point to the center */
	vecsub(center,partCfg[pi].r.p,pvector);

	/* Work out x and y coordinates with respect to rotation axis */
	
	/* Find the minimum distance of the point from the axis */
	vector_product(axis,pvector,vectprod);
	xdist = sqrt(sqrlen(vectprod)/sqrlen(axis));

	/* Find the projection of the vector from the point to the center
	   onto the axis vector */
	ydist = scalar(axis,pvector)/sqrt(sqrlen(axis));
	
    
	/* Work out relevant indices for x and y */
	xindex = (int)(floor(xdist/xbinwidth));
	yindex = (int)(floor((ydist+yrange*0.5)/ybinwidth));
	/*
	printf("x %d y %d \n",xindex,yindex);
	printf("p %f %f %f \n",partCfg[pi].r.p[0],partCfg[pi].r.p[1],partCfg[pi].r.p[2]);
	printf("pvec %f %f %f \n",pvector[0],pvector[1],pvector[2]);
	printf("axis %f %f %f \n",axis[0],axis[1],axis[2]);
	printf("dists %f %f \n",xdist,ydist);
	fflush(stdout);
	*/
	/* Check array bounds */
	if ( (xindex < xbins && xindex > 0) && (yindex < ybins && yindex > 0) ) {
	  density_map[bi].e[ybins*xindex+yindex] += 1;
	  xav += xdist;
	  yav += ydist;
	  beadcount += 1;
	} else {
	  //	    fprintf(stderr,"ERROR: outside array bounds in calc_radial_density_map"); fflush(NULL); errexit(); 
	}
      }

    }
  }


  /* Now turn counts into densities for the density map */
  for ( bi = 0 ; bi < nbeadtypes ; bi++ ) {
    for ( i = 0 ; i < xbins ; i++ ) {
      /* All bins are cylinders and therefore constant in yindex */
      binvolume = PI*(2*i*xbinwidth + xbinwidth*xbinwidth)*yrange;
      for ( j = 0 ; j < ybins ; j++ ) {
	density_map[bi].e[ybins*i+j] /= binvolume;
      }
    }
  }


  /* if required calculate the theta density profile */
  if ( thetabins > 0 ) {
    /* Convert the center to an output of the density center */
    xav = xav/(double)(beadcount);
    yav = yav/(double)(beadcount);
    thetabinwidth = 2*PI/(double)(thetabins);
    thetaradii = malloc(thetabins*nbeadtypes*sizeof(double));
    thetacounts = malloc(thetabins*nbeadtypes*sizeof(int));
    for ( bi = 0 ; bi < nbeadtypes ; bi++ ) {
      for ( t = 0 ; t < thetabins ; t++ ) {
	thetaradii[bi*thetabins+t] = 0.0;
	thetacounts[bi*thetabins+t] = 0.0;
      }
    }
    /* Maybe there is a nicer way to do this but now I will just repeat the loop over all particles */
      for ( pi = 0 ; pi < n_total_particles ; pi++ ) {
	for ( bi = 0 ; bi < nbeadtypes ; bi++ ) {
	  if ( beadids->e[bi] == partCfg[pi].p.type ) {
	    vecsub(center,partCfg[pi].r.p,pvector);
	    vector_product(axis,pvector,vectprod);
	    xdist = sqrt(sqrlen(vectprod)/sqrlen(axis));
	    ydist = scalar(axis,pvector)/sqrt(sqrlen(axis));
	    /* Center the coordinates */

	    xdist = xdist - xav;
	    ydist = ydist - yav;
	    rdist = sqrt(xdist*xdist+ydist*ydist);
	    if ( ydist >= 0 ) {
	      theta = acos(xdist/rdist);
	    } else {
	      theta = 2*PI-acos(xdist/rdist);
	    }
	    tindex = (int)(floor(theta/thetabinwidth));
	    thetaradii[bi*thetabins+tindex] += xdist + xav;
	    thetacounts[bi*thetabins+tindex] += 1;
	    if ( tindex >= thetabins ) {
	      fprintf(stderr,"ERROR: outside density_profile array bounds in calc_radial_density_map"); fflush(NULL); errexit(); 
	    } else {
	      density_profile[bi].e[tindex] += 1;
	    }
	  }	  
	}
      }



      /* normalize the theta densities*/
      for ( bi = 0 ; bi < nbeadtypes ; bi++ ) {
	for ( t = 0 ; t < thetabins ; t++ ) {
	  rdist = thetaradii[bi*thetabins+t]/(double)(thetacounts[bi*thetabins+t]);
	  density_profile[bi].e[t] /= rdist*rdist;
	}
      }
       


      free(thetaradii);
      free(thetacounts);

  }
  






  //  printf("done \n");
  return TCL_OK;
}

double calc_vanhove(int ptype, double rmin, double rmax, int rbins, int tmax, double *msd, double **vanhove) 
{
  int i, c1, c3, c3_max, ind, np=0;
  double p1[3],p2[3],dist;
  double bin_width, inv_bin_width;
  IntList p;

  /* create particle list */
  init_intlist(&p);
  for(i=0; i<n_total_particles; i++) { if( partCfg[i].p.type == ptype ) { np ++; } }
  if(np==0) { return 0; }
  alloc_intlist(&p,np);
  for(i=0; i<n_total_particles; i++) { if( partCfg[i].p.type == ptype ) { p.e[p.n]=i; p.n++; } }

  /* preparation */
  bin_width     = (rmax-rmin) / (double)rbins;
  inv_bin_width = 1.0 / bin_width;
 
  /* calculate msd and store distribution in vanhove */
  for(c1=0; c1<n_configs; c1++) { 
    c3_max=(c1+tmax+1)>n_configs ? n_configs : c1+tmax+1;
    for(c3=(c1+1); c3<c3_max; c3++) { 
      for(i=0; i<p.n; i++) {
	p1[0]=configs[c1][3*p.e[i] ]; p1[1]=configs[c1][3*p.e[i]+1]; p1[2]=configs[c1][3*p.e[i]+2];
	p2[0]=configs[c3][3*p.e[i]  ]; p2[1]=configs[c3][3*p.e[i]+1]; p2[2]=configs[c3][3*p.e[i]+2];
	dist = distance(p1, p2);
	if(dist > rmin && dist < rmax) {
	  ind = (int) ( (dist - rmin)*inv_bin_width );
	  vanhove[(c3-c1-1)][ind]++;
	}
	msd[(c3-c1-1)] += dist*dist;
      }
    }
  }

  /* normalize */
  for(c1=0; c1<(tmax); c1++) { 
    for(i=0; i<rbins; i++) { vanhove[c1][i] /= (double) (n_configs-c1-1)*p.n; }
    msd[c1] /= (double) (n_configs-c1-1)*p.n;
  }

  realloc_intlist(&p,0);
  return np;
}

/****************************************************************************************
 *                                 config storage functions
 ****************************************************************************************/

void analyze_append() {
  int i;
  n_part_conf = n_total_particles;
  configs = realloc(configs,(n_configs+1)*sizeof(double *));
  configs[n_configs] = (double *) malloc(3*n_part_conf*sizeof(double));
  for(i=0; i<n_part_conf; i++) {
    configs[n_configs][3*i]   = partCfg[i].r.p[0];
    configs[n_configs][3*i+1] = partCfg[i].r.p[1];
    configs[n_configs][3*i+2] = partCfg[i].r.p[2];
  }
  n_configs++;
}

void analyze_push() {
  int i;
  n_part_conf = n_total_particles;
  free(configs[0]);
  for(i=0; i<n_configs-1; i++) {
    configs[i]=configs[i+1];
  }
  configs[n_configs-1] = (double *) malloc(3*n_part_conf*sizeof(double));
  for(i=0; i<n_part_conf; i++) {
    configs[n_configs-1][3*i]   = partCfg[i].r.p[0];
    configs[n_configs-1][3*i+1] = partCfg[i].r.p[1];
    configs[n_configs-1][3*i+2] = partCfg[i].r.p[2];
  }
}

void analyze_replace(int ind) {
  int i;
  n_part_conf = n_total_particles;
  for(i=0; i<n_part_conf; i++) {
    configs[ind][3*i]   = partCfg[i].r.p[0];
    configs[ind][3*i+1] = partCfg[i].r.p[1];
    configs[ind][3*i+2] = partCfg[i].r.p[2];
  }
}

void analyze_remove(int ind) {
  int i;
  free(configs[ind]);
  for(i=ind; i<n_configs-1; i++) {
    configs[i]=configs[i+1];
  }
  n_configs--;
  configs = realloc(configs,n_configs*sizeof(double *));
  if (n_configs == 0) n_part_conf = 0;
}

void analyze_configs(double *tmp_config, int count) {
  int i;
  n_part_conf = count;
  configs = realloc(configs,(n_configs+1)*sizeof(double *));
  configs[n_configs] = (double *) malloc(3*n_part_conf*sizeof(double));
  for(i=0; i<n_part_conf; i++) {
    configs[n_configs][3*i]   = tmp_config[3*i];
    configs[n_configs][3*i+1] = tmp_config[3*i+1];
    configs[n_configs][3*i+2] = tmp_config[3*i+2];
  }
  n_configs++;
}

void analyze_activate(int ind) {
  int i;
  double pos[3];
  n_part_conf = n_total_particles;

  for(i=0; i<n_part_conf; i++) {
    pos[0] = configs[ind][3*i];
    pos[1] = configs[ind][3*i+1];
    pos[2] = configs[ind][3*i+2];
    if (place_particle(i, pos)==TCL_ERROR) {
      char *errtxt = runtime_error(128 + TCL_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt, "{057 failed upon replacing particle %d in Espresso} ", i); 
    }
  }
}


/****************************************************************************************
 *                                 Observables handling
 ****************************************************************************************/

void obsstat_realloc_and_clear(Observable_stat *stat, int n_pre, int n_bonded, int n_non_bonded,
			       int n_coulomb, int n_dipolar, int c_size)
{
  int i, total = c_size*(n_pre + n_bonded_ia + n_non_bonded + n_coulomb + n_dipolar);

  realloc_doublelist(&(stat->data), stat->data.n = total);
  stat->chunk_size = c_size;
  stat->n_coulomb    = n_coulomb;
  stat->n_dipolar    = n_dipolar;
  stat->n_non_bonded = n_non_bonded;
  stat->bonded     = stat->data.e + c_size*n_pre;
  stat->non_bonded = stat->bonded + c_size*n_bonded_ia;
  stat->coulomb    = stat->non_bonded + c_size*n_non_bonded;
  stat->dipolar    = stat->coulomb    + c_size*n_coulomb;

  for(i = 0; i < total; i++)
    stat->data.e[i] = 0.0;
}

void obsstat_realloc_and_clear_non_bonded(Observable_stat_non_bonded *stat_nb, int n_nonbonded, int c_size)
{
  int i, total = c_size*(n_nonbonded + n_nonbonded);

  realloc_doublelist(&(stat_nb->data_nb), stat_nb->data_nb.n = total);
  stat_nb->chunk_size_nb = c_size;
  stat_nb->n_nonbonded = n_nonbonded;
  stat_nb->non_bonded_intra = stat_nb->data_nb.e;
  stat_nb->non_bonded_inter = stat_nb->non_bonded_intra + c_size*n_nonbonded;
  
  for(i = 0; i < total; i++)
    stat_nb->data_nb.e[i] = 0.0;
}

void invalidate_obs()
{
  total_energy.init_status = 0;
  total_pressure.init_status = 0;
}


/****************************************************************************************
 *                                 basic observables parsing
 ****************************************************************************************/

static int parse_get_folded_positions(Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[10 + 3*TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i,change ;
  double shift[3];
  float *coord;

  enum flag { NONE , FOLD_MOLS};

  int flag;

  shift[0] = shift[1] = shift[2] = 0.0;
  flag = NONE;

  change = 0;
  shift[0] = shift[1] = shift[2] = 0.0;

  STAT_TRACE(fprintf(stderr,"%d,parsing get_folded_positions \n",this_node));
  while (argc > 0)
    {
      if ( ARG0_IS_S("-molecule") ) {
	flag = FOLD_MOLS;
	change = 1;
      }

      if ( ARG0_IS_S("shift") ) {
	if ( !ARG_IS_D(1,shift[0]) || !ARG_IS_D(2,shift[1]) || !ARG_IS_D(3,shift[2]) ) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp,"usage: analyze get_folded_positions [-molecule] [shift <xshift> <yshift> <zshift>]", (char *)NULL);
	  return (TCL_ERROR);
	}
	change = 4;
      }
      argc -= change;
      argv += change;
      STAT_TRACE(fprintf(stderr,"%d,argc = %d \n",this_node, argc));
    }

  freePartCfg();
  updatePartCfg(WITH_BONDS);
  if (!sortPartCfg()) {
    char *errtxt = runtime_error(128);
    ERROR_SPRINTF(errtxt, "{058 could not sort partCfg, particles have to start at 0 and have consecutive identities} ");
    return TCL_ERROR;
  }
  coord = malloc(n_total_particles*3*sizeof(float));
  /* Construct the array coord*/
  for (i = 0; i < n_total_particles; i++) {
    int dummy[3] = {0,0,0};
    double tmpCoord[3];
    tmpCoord[0] = partCfg[i].r.p[0];
    tmpCoord[1] = partCfg[i].r.p[1];
    tmpCoord[2] = partCfg[i].r.p[2];
    if (flag == NONE)  {   // perform folding by particle
      fold_position(tmpCoord, dummy);
    }    
    coord[i*3    ] = (float)(tmpCoord[0]);
    coord[i*3 + 1] = (float)(tmpCoord[1]);
    coord[i*3 + 2] = (float)(tmpCoord[2]);
  }


  // Use information from the analyse set command to fold chain molecules
  if ( flag == FOLD_MOLS ) {
    if( analyze_fold_molecules(coord, shift) != TCL_OK ){
      Tcl_AppendResult(interp, "could not fold chains: \"analyze set chains <chain_start> <n_chains> <chain_length>\" must be used first",
		       (char *) NULL);
      return (TCL_ERROR);;   
    }
  }

  //  Tcl_AppendResult(interp, "{ ", (char *)NULL);
  for ( i = 0 ; i < n_total_particles ; i++) {
    sprintf(buffer, " { %d %f %f %f } ", partCfg[i].p.identity , coord[i*3] , coord[i*3+1] , coord[i*3+2] );
    Tcl_AppendResult(interp, buffer , (char *)NULL);
  }
  //  Tcl_AppendResult(interp, "} ", (char *)NULL);

  return TCL_OK;

}

#ifdef MODES

static int parse_get_lipid_orients(Tcl_Interp *interp, int argc, char **argv)
{
  char buffer[TCL_DOUBLE_SPACE + TCL_INTEGER_SPACE];
  int i,change ;
  IntList l_orient;
  init_intlist(&l_orient);

  change = 0;

  STAT_TRACE(fprintf(stderr,"%d,parsing get_lipid_orients \n",this_node));
  /* Check that the grid has been set */
  if ( xdir + ydir + zdir == -3 ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"get_lipid_orients requires a grid. You must set grid first", (char *)NULL);
    return (TCL_ERROR);
  }
  
  realloc_intlist(&l_orient, n_molecules);
  get_lipid_orients(&l_orient);
  

  Tcl_AppendResult(interp, "{ Lipid_orientations } { ", (char *)NULL);
  for ( i = 0 ; i < n_molecules ; i++) {
    sprintf(buffer, "%d ", l_orient.e[i]);
    Tcl_AppendResult(interp, buffer , (char *)NULL);
  }
  Tcl_AppendResult(interp, "} ", (char *)NULL);

  realloc_intlist(&l_orient,0);

  return TCL_OK;

}

static int parse_modes2d(Tcl_Interp *interp, int argc, char **argv)
{
  STAT_TRACE(fprintf(stderr,"%d,parsing modes2d height grid \n",this_node);)
    /* 'analyze modes2d' */
    char buffer[TCL_DOUBLE_SPACE];
  int i,j,change ;
  fftw_complex* result_ht;
  fftw_complex* result_th;

  change = 0;

  /* Check that the grid has been set */
  if ( xdir + ydir + zdir == -3 ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"modes2d requires a grid. You must set grid first", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if (n_total_particles <= 2) {
    Tcl_AppendResult(interp, "(not enough particles for mode analysis)",
		     (char *)NULL);
    return (TCL_OK);
  }
  
  result_ht = malloc((mode_grid_3d[ydir]/2+1)*(mode_grid_3d[xdir])*sizeof(fftw_complex));
  result_th = malloc((mode_grid_3d[ydir]/2+1)*(mode_grid_3d[xdir])*sizeof(fftw_complex));

  if (!modes2d(result_th, 0) || !modes2d(result_ht,1)) {
    fprintf(stderr,"%d,mode analysis failed \n",this_node);
    return TCL_ERROR;
  }
  else {    STAT_TRACE(fprintf(stderr,"%d,mode analysis done \n",this_node));}
  

  /* Output is of the form :
     ht_RE ht_IM th_RE th_IM
  */
  Tcl_AppendResult(interp, "{ Modes } { ", (char *)NULL);
  for ( i = 0 ; i < mode_grid_3d[xdir] ; i++) {
    Tcl_AppendResult(interp, " { ", (char *)NULL);
    for ( j = 0 ; j < mode_grid_3d[ydir]/2 + 1 ; j++) {
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      Tcl_PrintDouble(interp,FFTW_REAL(result_ht[j+i*(mode_grid_3d[ydir]/2+1)]),buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " ", (char *)NULL);
      Tcl_PrintDouble(interp,FFTW_IMAG(result_ht[j+i*(mode_grid_3d[ydir]/2+1)]),buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " ", (char *)NULL);
      Tcl_PrintDouble(interp,FFTW_REAL(result_th[j+i*(mode_grid_3d[ydir]/2+1)]),buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " ", (char *)NULL);
      Tcl_PrintDouble(interp,FFTW_IMAG(result_th[j+i*(mode_grid_3d[ydir]/2+1)]),buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " } ", (char *)NULL);
    }
    Tcl_AppendResult(interp, " } ", (char *)NULL);
  }


  Tcl_AppendResult(interp, " } ", (char *)NULL);

  free(result_ht);
  free(result_th);

  return TCL_OK;

}

static int parse_bilayer_set(Tcl_Interp *interp, int argc, char **argv)
{
  int change;
  change = 0;
  while ( argc > 0 ) {
    if ( ARG0_IS_S("grid") ) { 
      if ( !ARG_IS_I(1,mode_grid_3d[0]) || !ARG_IS_I(2,mode_grid_3d[1]) || !ARG_IS_I(3,mode_grid_3d[2]) ) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp,"usage: analyze bilayer_set [grid <xdim> <ydim> <zdim>] [straycut <stray_cut_off>]", (char *)NULL);
	return (TCL_ERROR);
      }
      STAT_TRACE(fprintf(stderr,"%d,set has args %d,%d,%d \n",this_node,mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]));
      change = 4;
      /* Update global parameters */
      map_to_2dgrid();
      mode_grid_changed = 1;
      argc -= change;
      argv += change;      
    }
    
    if ( ARG0_IS_S("straycut") ) { 
      if ( !ARG_IS_D(1,stray_cut_off) ) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp,"usage: analyze bilayer_set [grid <xdim> <ydim> <zdim>] [straycut <stray_cut_off>]", (char *)NULL);
	return (TCL_ERROR);
      }
      STAT_TRACE(fprintf(stderr,"%d,setgrid has args %d,%d,%d \n",this_node,mode_grid_3d[0],mode_grid_3d[1],mode_grid_3d[2]));
      change = 2;
      argc -= change;
      argv += change;
    }
    
  }
  return TCL_OK;
}


static int parse_radial_density_map(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze radial density profile ' */
  char buffer[TCL_DOUBLE_SPACE+256];
  int i,j,k;
  IntList beadtypes;
  double rotation_axis[3];
  double rotation_center[3];
  int xbins,ybins,thetabins,dobyangle;
  double xrange,yrange;
  DoubleList *density_profile = NULL;
  DoubleList *density_map = NULL;
  dobyangle = 0;
  thetabins = 0;
  init_intlist(&beadtypes);
  alloc_intlist(&beadtypes,1);



  if ( !ARG_IS_I(0,xbins) || !ARG_IS_I(1,ybins) || !ARG_IS_D(2,xrange) || !ARG_IS_D(3,yrange) ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"first 4 args usage:  analyze radial_density_map <xbins> <ybins> <xrange> <yrange> {<axisofrotation>} {<centerofrotation>} { { <beadtypelist> } [thetabins]", (char *)NULL);
    return (TCL_ERROR);
  } else {
    argc -= 4;
    argv += 4;
  }

  if ( !ARG_IS_D(0,rotation_axis[0]) || !ARG_IS_D(1,rotation_axis[1]) || !ARG_IS_D(2,rotation_axis[2]) ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"rotation_axis usage:  analyze radial_density_map <xbins> <ybins> <xrange> <yrange> {<axisofrotation>} {<centerofrotation>} { { <beadtypelist> } [thetabins]", (char *)NULL);
    return (TCL_ERROR);
  } else {
    argc -= 3;
    argv += 3;
  }
  
  if ( !ARG_IS_D(0,rotation_center[0] ) || !ARG_IS_D(1,rotation_center[1] ) || !ARG_IS_D(2,rotation_center[2] ) ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"rotation_center usage: analyze radial_density_map <xbins> <ybins> <xrange> <yrange> {<axisofrotation>} {<centerofrotation>} { { <beadtypelist> } [thetabins]", (char *)NULL);
    return (TCL_ERROR);
  } else {
    argc -= 3;
    argv += 3;
  }

  if ( !ARG_IS_INTLIST(0,beadtypes)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp," beadtypes usage: analyze radial_density_map <xbins> <ybins> <xrange> <yrange> {<axisofrotation>} {<centerofrotation>} { { <beadtypelist> } [thetabins] ", (char *)NULL);
    return (TCL_ERROR);
  } else {
    argc -= 1;
    argv += 1;
  }

  /* Now check for the optional thetabins argument */
  if ( argc > 0 ) {
    if ( !ARG_IS_I(0,thetabins)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp,"thetabins usage:  analyze radial_density_map <xbins> <ybins> <xrange> <yrange> {<axisofrotation>} {<centerofrotation>} { { <beadtypelist> } [thetabins] ", (char *)NULL);
      return (TCL_ERROR);
    } else {
      dobyangle = 1;
      argc -= 1;
      argv += 1;
    }
  }

  
  printf("nbtypes %d \n",beadtypes.max);
  printf("bt %d \n",beadtypes.e[0]);
  printf("bins %d %d %d %f %f \n",xbins,ybins,thetabins,xrange,yrange);
  printf("rotationaxis %f %f %f \n",rotation_axis[0],rotation_axis[1],rotation_axis[2]);
  printf("center %f %f %f \n",rotation_center[0],rotation_center[1],rotation_center[2]);

  /* allocate memory for the profile if necessary */
  if (thetabins > 0 ) {
    density_profile = malloc(beadtypes.max*sizeof(DoubleList));
    if (density_profile) {
      for ( i = 0 ; i < beadtypes.max ; i++ ) {
	init_doublelist(&density_profile[i]);
	alloc_doublelist(&density_profile[i],thetabins);
	for ( j = 0 ; j < thetabins ; j++ ) {
	  density_profile[i].e[j] = 0.0;
	}
      }
    } else {
	Tcl_AppendResult(interp,"could not allocate memory for density_profile", (char *)NULL);
	return (TCL_ERROR);
    }
  }
  /* Allocate a doublelist of bins for each beadtype so that we
     can keep track of beads separately */
  density_map = malloc(beadtypes.max*sizeof(DoubleList));
  if ( density_map ) {
  /* Initialize all the subprofiles in density profile */
    for ( i = 0 ; i < beadtypes.max ; i++ ) {
      //      printf("initializing for beadtype %d and i: %d \n",beadtypes.e[i],i);
      init_doublelist(&density_map[i]);
      alloc_doublelist(&density_map[i],xbins*ybins);
      if (!density_map[i].e) { 
	Tcl_AppendResult(interp,"could not allocate memory for density_map", (char *)NULL);
	return (TCL_ERROR);
      } else {
	for ( j = 0 ; j < xbins*ybins ; j++ ) {
	  // printf("j: %d max: %d",j,density_map[i].max);
	  density_map[i].e[j] = 0.0;
	}

      }
      //      printf("done initializing at i: %d j: %d\n",i,j);
    }
  } else {
    Tcl_AppendResult(interp,"could not allocate memory for density_map", (char *)NULL);
    return (TCL_ERROR);
  }

  //  printf("about to calculate profile \n");
  if(calc_radial_density_map(xbins,ybins,thetabins,xrange,yrange,rotation_axis,rotation_center,&beadtypes, density_map ,density_profile) != TCL_OK) {
    Tcl_AppendResult(interp, "Error calculating radial density profile ", (char *)NULL);
    return TCL_ERROR;
  }

  //  printf("done calculating profile \n");

  Tcl_AppendResult(interp, "{ RadialDensityMap  { ", (char *)NULL);
  for ( i = 0 ; i < beadtypes.n ; i++) {
    Tcl_AppendResult(interp, " { ", (char *)NULL);
    for ( j = 0 ; j < xbins ; j++) {
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      for ( k = 0 ; k < ybins ; k++) {
	Tcl_PrintDouble(interp,density_map[i].e[j*ybins+k],buffer);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
	Tcl_AppendResult(interp, " ", (char *)NULL);
      }
      Tcl_AppendResult(interp, " } ", (char *)NULL);
    }
    Tcl_AppendResult(interp, " } ", (char *)NULL);
  }
  Tcl_AppendResult(interp, " } ", (char *)NULL);    
  Tcl_AppendResult(interp, " } ", (char *)NULL); 
  if ( thetabins > 0 ) {
    Tcl_AppendResult(interp, "{ RadialDensityThetaProfile  { ", (char *)NULL);
    for ( i = 0 ; i < beadtypes.n ; i++) {
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      for ( j = 0 ; j < thetabins ; j++ ) {
	Tcl_PrintDouble(interp,density_profile[i].e[j],buffer);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
	Tcl_AppendResult(interp, " ", (char *)NULL);
      }
      Tcl_AppendResult(interp, " } ", (char *)NULL);
    }
    Tcl_AppendResult(interp, " } ", (char *)NULL);
    Tcl_AppendResult(interp, " } ", (char *)NULL);

    /* Dealloc the density profile*/
    for ( i = 0 ; i < beadtypes.max ; i++ ) {
      realloc_doublelist(&density_profile[i],0);
    }
    free(density_profile);
  }


  /* Dealloc all subprofiles in density map */
  for ( i = 0 ; i < beadtypes.max ; i++ ) {
    realloc_doublelist(&density_map[i],0);
  }
  free(density_map);

  realloc_intlist(&beadtypes,0);
  return TCL_OK;
}


static int parse_bilayer_density_profile(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze bilayer density profile ' */
  char buffer[TCL_DOUBLE_SPACE];
  int change,nbins,i,j;
  IntList beadtypes;
  double hrange;
  int usegrid;
  double sradius, cx,cy,cz;
  DoubleList* density_profile = NULL;
  int errorValue;
  double scenter[3] = { cx, cy ,cz };
  /* Defaults */
  sradius = 0.0;
  cx = cy = cz = 0.0;

  init_intlist(&beadtypes);
  usegrid = 1;
  if ( !ARG0_IS_D(hrange) ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"usage: analyze bilayer_density_profile <hrange> <nbins> <beadtype1> .. <beadtypeN> [withsphere <radius> <cx> <cy> <cz> ] [nogrid]", (char *)NULL);
    return (TCL_ERROR);
  } else {
    argc -= 1;
    argv += 1;
  }

  if ( !ARG0_IS_I(nbins) ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"usage: analyze bilayer_density_profile <hrange> <nbins> <beadtype1> .. <beadtypeN> [withsphere <radius> <cx> <cy> <cz> ] [nogrid]", (char *)NULL);
    return (TCL_ERROR);
  } else {
    argc -= 1;
    argv += 1;
  }

  if ( !ARG0_IS_INTLIST(beadtypes) ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"usage: analyze bilayer_density_profile <hrange> <nbins> <beadtype1> .. <beadtypeN> [withsphere <radius> <cx> <cy> <cz> ] [nogrid]", (char *)NULL);
    return (TCL_ERROR);
  } else {
    change = 1;
    argc -= change;
    argv += change;
  }


  while ( argc > 0 ) {

    if ( ARG0_IS_S("withsphere") ) {
      if ( !ARG_IS_D(1,sradius) || !ARG_IS_D(2,cx) || !ARG_IS_D(3,cy) || !ARG_IS_D(4,cz) ) {
	Tcl_ResetResult(interp);
	Tcl_AppendResult(interp,"usage: analyze bilayer_density_profile <hrange> <nbins> <beadtype1> .. <beadtypeN> [withsphere <radius> <cx> <cy> <cz> ] [nogrid]", (char *)NULL);
	return (TCL_ERROR);
      }
      change = 5;
    }

    if ( ARG0_IS_S("nogrid") ) { 
      usegrid = 0;
      change = 1;
    }
    argc -= change;
    argv += change;
  }


  /* If we said we wanted to use a grid then check to make sure one is set */
  if ( xdir + ydir + zdir == -3 && usegrid ) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp,"analyze bilayer_density_profile requires a grid. You must set grid first", (char *)NULL);
    return (TCL_ERROR);
  }



  /* Allocate a two doublelists of bins for each beadtype so that we
     can keep track of beads in upper or lower lipids */
  density_profile = (DoubleList *)malloc(beadtypes.max*2*sizeof(DoubleList));
  if ( density_profile ) {
  /* Initialize all the subprofiles in density profile */
    for ( i = 0 ; i < beadtypes.max*2 ; i++ ) {
      init_doublelist(&density_profile[i]);
      alloc_doublelist(&density_profile[i],nbins);
      for ( j = 0 ; j < nbins ; j++ ) {
	density_profile[i].e[j] = 0.0;
      }
    }
  } else {
    Tcl_AppendResult(interp,"could not allocate memory for density_profile", (char *)NULL);
    return (TCL_ERROR);
  }

  if ( sradius > 0 ) {
    errorValue = bilayer_density_profile_sphere(&beadtypes, hrange, density_profile,sradius,scenter );
  } else {
    errorValue = bilayer_density_profile(&beadtypes, hrange, density_profile, usegrid );
  }

  if (!errorValue) {
    Tcl_AppendResult(interp, "Error calculating bilayer density profile ", (char *)NULL);
    return TCL_ERROR;
  }

  Tcl_AppendResult(interp, "{ DensityProfile } { ", (char *)NULL);
  for ( i = 0 ; i < beadtypes.n*2 ; i++) {
    Tcl_AppendResult(interp, " { ", (char *)NULL);
    for ( j = 0 ; j < nbins ; j++) {
      Tcl_PrintDouble(interp,density_profile[i].e[j],buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " ", (char *)NULL);
    }
    Tcl_AppendResult(interp, " } ", (char *)NULL);
  }
  Tcl_AppendResult(interp, " } ", (char *)NULL);    

  /* Initialize all the subprofiles in density profile */
  for ( i = 0 ; i < beadtypes.max*2 ; i++ ) {
    realloc_doublelist(&density_profile[i],0);
  }
  free(density_profile);
  realloc_intlist(&beadtypes,0);
  return TCL_OK;
}

static int parse_lipid_orient_order(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze lipid_orient_order ' */
  double result;
  double* stored_dirs;
  char buffer[TCL_DOUBLE_SPACE];
  int i , j ;
  result = 0;

  if (n_total_particles <= 1) {
    Tcl_AppendResult(interp, "(not enough particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  stored_dirs = malloc(sizeof(double)*n_molecules*3);
  /* Do the calculation */
  if ( orient_order(&result,stored_dirs) != TCL_OK ) {
    Tcl_AppendResult(interp, "Error calculating orientational order ", (char *)NULL);
    return TCL_ERROR;
  }

  if ( argc == 0 ) {
    /* If no options are specified then only give the average
       orientational order */
      Tcl_PrintDouble(interp, result, buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      return TCL_OK;
  } else {
    /* If the -all option is specified then print everything */
    if ( ARG0_IS_S("all") ) {
      Tcl_PrintDouble(interp, result, buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      Tcl_AppendResult(interp, " { ", (char *)NULL);
      for ( i = 0 ; i < n_molecules ; i++ ) {
	Tcl_AppendResult(interp, " { ", (char *)NULL);
	for ( j = 0 ; j < 3 ; j++ ) {
	  Tcl_PrintDouble(interp,stored_dirs[i*3+j],buffer);
	  Tcl_AppendResult(interp, buffer, (char *)NULL);
	  Tcl_AppendResult(interp, " ", (char *)NULL);
	}
	Tcl_AppendResult(interp, "} ", (char *)NULL);
      } 
      Tcl_AppendResult(interp, "} ", (char *)NULL);
    }
    return TCL_OK;
  }
  
  Tcl_AppendResult(interp, "Error calculating orientational order ", (char *)NULL);
  free(stored_dirs);
  return TCL_ERROR;
}

#endif

static int parse_aggregation(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze centermass [<type>]' */
  char buffer[256 + 3*TCL_INTEGER_SPACE + 2*TCL_DOUBLE_SPACE];
  int i, target1;
  int *agg_id_list;
  double dist_criteria, dist_criteria2;
  int charge_criteria, min_contact, *head_list, *link_list;
  int agg_num =0, *agg_size, agg_min= n_molecules, agg_max = 0,  agg_std = 0, agg_avg = 0; 
  float fagg_avg;
  int s_mol_id, f_mol_id;

  agg_id_list = (int *) malloc(n_molecules *sizeof(int));
  head_list =  (int *) malloc(n_molecules *sizeof(int));
  link_list = (int *) malloc(n_molecules *sizeof(int));
  agg_size = (int *) malloc(n_molecules *sizeof(int));

  /* parse arguments */
  if (argc < 3) {
    Tcl_AppendResult(interp, "usage: analyze aggregation <dist_criteria> <start mol_id> <finish mol_id> [<min_contact>] [<charge_criteria>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG_IS_D(0,dist_criteria)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze aggregation <dist_criteria> <start mol_id> <finish mol_id> [<min_contact>] [<charge_criteria>]", (char *)NULL);
    return (TCL_ERROR);
  }
  dist_criteria2 = dist_criteria * dist_criteria;

  if (!ARG_IS_I(1,s_mol_id)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze aggregation <dist_criteria> <start mol_id> <finish mol_id> [<min_contact>] [<charge_criteria>]", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG_IS_I(2,f_mol_id)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze aggregation <dist_criteria> <start mol_id> <finish mol_id> [<min_contact>] [<charge_criteria>]", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if (n_nodes > 1) {
    Tcl_AppendResult(interp, "aggregation can only be calculated on a single processor", (char *)NULL);
    return TCL_ERROR;
  }

  if (cell_structure.type != CELL_STRUCTURE_DOMDEC) {
    Tcl_AppendResult(interp, "aggregation can only be calculated with the domain decomposition cell system", (char *)NULL);
    return TCL_ERROR;
  }

  if ( (s_mol_id < 0) || (s_mol_id > n_molecules) || (f_mol_id < 0) || (f_mol_id > n_molecules) ) {
    Tcl_AppendResult(interp, "check your start and finish molecule id's", (char *)NULL);
    return TCL_ERROR;
  }

  if ( max_range_non_bonded2 < dist_criteria2) {
    Tcl_AppendResult(interp, "dist_criteria is larger than max_range_non_bonded.", (char *)NULL);
    return TCL_ERROR;    
    
  }

  if (argc == 4) {
      if (!ARG_IS_I(3,min_contact)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "usage: analyze aggregation <dist_criteria> <start mol_id> <finish mol_id> [<min_contact>] [<charge_criteria>]", (char *)NULL);
	  return (TCL_ERROR);
      }
  } else {
      min_contact = 1;
  }

  if (argc == 5) {
      if (!ARG_IS_I(4, charge_criteria)) {
	  Tcl_ResetResult(interp);
	  Tcl_AppendResult(interp, "usage: analyze aggregation <dist_criteria> <start mol_id> <finish mol_id> [<min_contact>] [<charge_criteria>]", (char *)NULL);
	  return (TCL_ERROR);
      }
  } else {
      charge_criteria = 0;
  }


  aggregation(dist_criteria2, min_contact, s_mol_id, f_mol_id, head_list, link_list, agg_id_list, 
	      &agg_num, agg_size, &agg_max, &agg_min, &agg_avg, &agg_std, charge_criteria);

  fagg_avg = (float) (agg_avg)/agg_num;
  sprintf (buffer, " MAX %d MIN %d AVG %f STD %f AGG_NUM %d AGGREGATES", 
	   agg_max, agg_min, fagg_avg, sqrt( (float) (agg_std/(float)(agg_num)-fagg_avg*fagg_avg)), agg_num);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  
  for (i = s_mol_id ; i <= f_mol_id ; i++) {
    if (head_list[i] != -2) {
      target1= head_list[i];
      sprintf(buffer, " { %d ", target1); 
      Tcl_AppendResult(interp, buffer, (char *)NULL);
      while( link_list[target1] != -1) {
	target1= link_list[target1];
	sprintf(buffer, "%d ", target1); 
	Tcl_AppendResult(interp, buffer, (char *)NULL); 
      }
      sprintf(buffer, "} ");
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
  }

  free(agg_id_list);
  free(head_list);
  free(link_list);
  free(agg_size);

  return TCL_OK;
}


static int parse_mindist(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze mindist [<type_list_a> <type_list_b>]' */
  double result;
  char buffer[TCL_DOUBLE_SPACE];
  IntList p1, p2;

  init_intlist(&p1); init_intlist(&p2);

  if (n_total_particles <= 1) {
    Tcl_AppendResult(interp, "(not enough particles)",
		     (char *)NULL);
    return (TCL_OK);
  }
  if (argc == 0)
    result = mindist(NULL, NULL);
  else {
    /* parse arguments */
    if (argc < 2) {
      Tcl_AppendResult(interp, "usage: analyze mindist [<type_list> <type_list>]", (char *)NULL);
      return (TCL_ERROR);
    }

    if (!ARG0_IS_INTLIST(p1)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze mindist [<type_list> <type_list>]", (char *)NULL);
      return (TCL_ERROR);
    }
    if (!ARG1_IS_INTLIST(p2)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze mindist [<type_list> <type_list>]", (char *)NULL);
      return (TCL_ERROR);
    }
    result = mindist(&p1, &p2);
  }

  Tcl_PrintDouble(interp, result, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

static int parse_centermass(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze centermass [<type>]' */
  double com[3];
  char buffer[3*TCL_DOUBLE_SPACE+3];
  int p1;
  
  /* parse arguments */
  if (argc != 1) {
    Tcl_AppendResult(interp, "usage: analyze centermass [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze centermass [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  
  centermass(p1, com);
  
  sprintf(buffer,"%f %f %f",com[0],com[1],com[2]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

static int parse_angularmomentum(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze angularmomentum [<type>]' */
  double com[3];
  char buffer[3*TCL_DOUBLE_SPACE+3];
  int p1;

  /* parse arguments */
  if (argc != 1) {
    Tcl_AppendResult(interp, "usage: analyze angularmomentum [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze angularmomentum [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  angularmomentum(p1, com);

  sprintf(buffer,"%f %f %f",com[0],com[1],com[2]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

  //subfunction: mark all neighbors of a particle and their neighbors (recursiv!)
void mark_neighbours(int type,int pa_nr,double dist,int *list){
  int k;
  for (k=0;k<n_total_particles;k++){
     //only unmarked and particles with right distance
     if ( (partCfg[k].p.type == type) && (list[k] == 0) && (min_distance(partCfg[pa_nr].r.p,partCfg[k].r.p) < dist) ){
        //mark particle with same number as calling particle
        list[k]=list[pa_nr];
        mark_neighbours(type,k,dist,list);
     }
  }
}

static int parse_cluster_size_dist(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze cluster_size_dist [<type>]' */
  char buffer[3*TCL_DOUBLE_SPACE+3];
  int p1;
  double dist;
  int i,j,max_cluster_number,size;
  int cluster_number[n_total_particles];//cluster number of particle number
  int cluster_size[n_total_particles+1]; //number of clusters with some size

  /* parse arguments */
  if (argc != 2) {
    Tcl_AppendResult(interp, "usage: analyze cluster_size_dist <type> <dist>", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "ERROR: type must be int !", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG1_IS_D(dist)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "ERROR: dist must be double !", (char *)NULL);
    return (TCL_ERROR);
  }

  updatePartCfg(WITHOUT_BONDS);
  for (i=0;i<n_total_particles;i++){
     cluster_number[i]=0;
     cluster_size[i]=0;
  }
  cluster_size[n_total_particles]=0;

  max_cluster_number=1;
  for (i=0;i<n_total_particles-1;i++){
     //if particle was not marked till now
     if ( (partCfg[i].p.type == p1) && (cluster_number[i]==0) ){
       //mark current particle with max_cluster_number
       cluster_number[i]=max_cluster_number;
       mark_neighbours(p1,i,dist,cluster_number);
       max_cluster_number++;
     }
  }

  for (i=0;i<n_total_particles;i++){
    if ( (partCfg[i].p.type == p1) && (cluster_number[i]==0) ) {
       Tcl_ResetResult(interp);
       Tcl_AppendResult(interp, "ERROR: at least one particle is not marked !", (char *)NULL);
       return (TCL_ERROR);
    }
  }

  for (i=1;i<max_cluster_number;i++){
     size=0;
     for (j=0;j<n_total_particles;j++){
        //Finding particles belonging to cluster i
        if ( (partCfg[j].p.type == p1) && (cluster_number[j]==i)) {size++;}
     }
     cluster_size[size]++;
  }

  sprintf(buffer,"%i %f",p1,dist);
  Tcl_AppendResult(interp, "{ analyze cluster_size_dist ",buffer,"} {\n",(char *)NULL);
  for (i=0;i<n_total_particles+1;i++){
      if (cluster_size[i]!=0){
         sprintf(buffer,"%i %i",i,cluster_size[i]);
         Tcl_AppendResult(interp, "{ ",buffer," }\n", (char *)NULL);
      }
  }
  Tcl_AppendResult(interp, "}",(char *)NULL);

  return TCL_OK;
}

static int parse_momentofinertiamatrix(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze  momentofinertiamatrix [<type>]' */
  double MofImatrix[9];
  char buffer[9*TCL_DOUBLE_SPACE+9];
  int p1;

  /* parse arguments */
  if (argc != 1) {
    Tcl_AppendResult(interp, "usage: analyze momentofinertiamatrix [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze momentofinertiamatrix [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  momentofinertiamatrix(p1, MofImatrix);
  
  sprintf(buffer,"%f %f %f %f %f %f %f %f %f",
	  MofImatrix[0],MofImatrix[1],MofImatrix[2],MofImatrix[3],MofImatrix[4],
	  MofImatrix[5],MofImatrix[6],MofImatrix[7],MofImatrix[8]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

static int parse_find_principal_axis(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze find_principal_axis [<type0>]' */
  double MofImatrix[9],eva[3],eve[3];
  char buffer[4*TCL_DOUBLE_SPACE+20];
  int p1,i,j;

  /* parse arguments */
  if (argc != 1) {
    Tcl_AppendResult(interp, "usage: analyze find_principal_axis [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze find_principal_axis [<type>]", (char *)NULL);
    return (TCL_ERROR);
  }

  momentofinertiamatrix(p1, MofImatrix);
  i=calc_eigenvalues_3x3(MofImatrix, eva);
  
  sprintf(buffer,"{eigenval eigenvector} ");
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  for (j= 0; j < 3; j++) {
    i=calc_eigenvector_3x3(MofImatrix,eva[j],eve);
    sprintf(buffer," { %f { %f %f %f } }",eva[j],eve[0],eve[1],eve[2]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  return TCL_OK;
}

static int parse_nbhood(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze nbhood [-planar <x> <y> <z>] { <partid> | <posx> <posy> <posz> } <r_catch> ' */
  int p, i;
  double pos[3];
  double r_catch;
  char buffer[TCL_INTEGER_SPACE + 2];  
  IntList il;
  int planedims[3], change;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }
  planedims[0] = planedims[1] = planedims[2] = 1;

  /* Optional planar argument */
  if ( ARG0_IS_S("planar")) {      
    if ( argc < 6 ) {
      Tcl_AppendResult(interp, "planar option requires exactly 3 integer arguments",
		       (char *)NULL); 
      return TCL_ERROR;
    }
      
    /* Now read in the dimensions to be used */
    if (! ARG_IS_I(1,planedims[0] ))
      return TCL_ERROR;
    if (! ARG_IS_I(2,planedims[1] ))
      return TCL_ERROR;
    if (! ARG_IS_I(3,planedims[2] ))
      return TCL_ERROR;
      
    /* Check to make sure that we actually specified a plane */
    if ( (planedims[0] + planedims[1] + planedims[2] ) != 2 ) {
      Tcl_AppendResult(interp, "you specified the planar option. Please specify exactly 2 dimensions eg 1 1 0 for xy plane",
		       (char *)NULL); 
      return TCL_ERROR;
    }
    change = 4;
    argc -= change;
    argv += change;
    
  }

  /* Process obligatory arguments */
  get_reference_point(interp, &argc, &argv, pos, &p);
  if (!ARG0_IS_D(r_catch)) {
    Tcl_AppendResult(interp, "usage: nbhood [planar <x> <y> <z>] { <partid> | <posx> <posy> <posz> } <r_catch> ",
		     (char *)NULL);    
    return (TCL_ERROR);
  }
  argc--;
  argv++;

  updatePartCfg(WITHOUT_BONDS);
  nbhood(pos, r_catch, &il, planedims );

  
  for (i = 0; i < il.n; i++) {
    sprintf(buffer, "%d ", il.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  realloc_intlist(&il, 0);
  return (TCL_OK);
}

static int parse_distto(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze distto { <part_id> | <posx> <posy> <posz> }' */
  double result;
  int p;
  double pos[3];
  char buffer[TCL_DOUBLE_SPACE], usage[150];
  sprintf(usage, "distto { <partid> | <posx> <posy> <posz> }");

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp, "(no particles)",
		     (char *)NULL);
    return (TCL_OK);
  }

  if (argc == 0) {
    Tcl_AppendResult(interp, "usage: ", usage, (char*)NULL);
    return TCL_ERROR;
  }

  get_reference_point(interp, &argc, &argv, pos, &p);
  if (argc != 0) {
    Tcl_AppendResult(interp, "usage: ", usage, (char *)NULL);
    return TCL_ERROR;
  }

  updatePartCfg(WITHOUT_BONDS);

  result = distto(pos, p);

  Tcl_PrintDouble(interp, result, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (TCL_OK);
}


static int parse_cell_gpb(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze cell_gpb <Manning parameter> <outer cell radius> <inner cell radius> [<accuracy> [<# of interations>]]' */
  double result[3] = {0, 0, 0}, xi_m, Rc, ro;
  double gacc = 1e-6;
  int maxtry  = 30000;
  char buffer[3*TCL_DOUBLE_SPACE+20], usage[150];
  sprintf(usage,"analyze cell_gpb <Manning parameter> <outer cell radius> <inner cell radius> [<accuracy> [<# of interations>]]");

  if ((argc < 3) || (argc > 5)) { 
    Tcl_AppendResult(interp, "usage: ",usage,(char *)NULL); return TCL_ERROR; }
  else if (!ARG_IS_D(0,xi_m) || !ARG_IS_D(1,Rc) || !ARG_IS_D(2,ro))
    return TCL_ERROR;
  if (argc == 4) if (!ARG_IS_D(3,gacc))
    return TCL_ERROR;
  if (argc == 5) if (!ARG_IS_I(4,maxtry))
    return TCL_ERROR;
  if ((xi_m < 0) || !((Rc > 0) && (ro > 0))) {
    Tcl_ResetResult(interp); sprintf(buffer,"%f %f %f",xi_m,Rc,ro);
    Tcl_AppendResult(interp, "usage: ",usage,"\n",(char *)NULL);
    Tcl_AppendResult(interp, "ERROR: All three arguments must be positive, the latter two even non-zero (got: ",buffer,")! Aborting...", (char*)NULL);
    return(TCL_ERROR);
  }

  calc_cell_gpb(xi_m,Rc,ro,gacc,maxtry,result);

  if (result[2] == -2.0) {
    Tcl_ResetResult(interp); sprintf(buffer,"%d",maxtry); Tcl_AppendResult(interp, "ERROR: Maximum number of iterations exceeded (",buffer,")! ");
    sprintf(buffer,"%f and %f",result[0],result[1]);      Tcl_AppendResult(interp, "Got ",buffer," so far, aborting now...", (char*)NULL);
    return(TCL_ERROR); 
  } else if (result[2] == -3.0) {
    Tcl_ResetResult(interp); sprintf(buffer,"%f and %f",result[0],result[1]);
    Tcl_AppendResult(interp, "ERROR: gamma is not bracketed by the programs initial guess (",buffer,")! Aborting...", (char*)NULL);
    return(TCL_ERROR); 
  } else if (result[2] == -4.0) {
    Tcl_ResetResult(interp); sprintf(buffer,"%f and %f",result[0],result[1]);
    Tcl_AppendResult(interp, "ERROR: lower boundary on wrong side of the function (",buffer,")! Aborting...", (char*)NULL);
    return(TCL_ERROR); 
  } else if (result[2] == -5.0) {
    Tcl_ResetResult(interp); Tcl_AppendResult(interp, "ERROR: Something went wrong! Aborting...", (char*)NULL);
    return(TCL_ERROR); 
  }
  sprintf(buffer,"%f %f %f",result[0],result[1],result[2]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (TCL_OK);
}


static int parse_Vkappa(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze Vkappa [{ reset | read | set <Vk1> <Vk2> <avk> }]' */
  double result = 0.0;
  char buffer[3*TCL_DOUBLE_SPACE+3];

  if (argc > 0)
    if (ARG0_IS_S("reset"))
      Vkappa.Vk1 = Vkappa.Vk2 = Vkappa.avk = 0.0;
    else if (ARG0_IS_S("read")) {
      sprintf(buffer,"%f %f %f ",Vkappa.Vk1,Vkappa.Vk2,Vkappa.avk);
      Tcl_AppendResult(interp, buffer, (char *)NULL); return (TCL_OK); }
    else if (ARG0_IS_S("set")) {
      if (argc < 4 || !ARG_IS_D(1,Vkappa.Vk1) || !ARG_IS_D(2,Vkappa.Vk2) || !ARG_IS_D(3,Vkappa.avk)) {
        Tcl_AppendResult(interp, "usage: analyze Vkappa [{ reset | read | set <Vk1> <Vk2> <avk> }] ", (char *)NULL);  return TCL_ERROR;  }
      if (Vkappa.avk <= 0.0) {
        Tcl_AppendResult(interp, "ERROR: # of averages <avk> must be positiv! Resetting values...", (char *)NULL);
        result = Vkappa.Vk1 = Vkappa.Vk2 = Vkappa.avk = 0.0; return TCL_ERROR; }
      result = Vkappa.Vk2/Vkappa.avk - SQR(Vkappa.Vk1/Vkappa.avk); }
    else {
      Tcl_AppendResult(interp, "usage: analyze Vkappa [{ reset | read | set <Vk1> <Vk2> <avk> }] ", (char *)NULL);  return TCL_ERROR;  }
  else {
    Vkappa.Vk1 += box_l[0]*box_l[1]*box_l[2];
    Vkappa.Vk2 += SQR(box_l[0]*box_l[1]*box_l[2]);
    Vkappa.avk += 1.0;
    result = Vkappa.Vk2/Vkappa.avk - SQR(Vkappa.Vk1/Vkappa.avk);
  }

  Tcl_PrintDouble(interp, result, buffer);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (TCL_OK);
}


static int parse_distribution(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze distribution { <part_type_list_a> } { <part_type_list_b> } [<r_min> [<r_max> [<r_bins> [<log_flag> [<int_flag>]]]]]' */
  /*********************************************************************************************************************************/
  char buffer[2*TCL_DOUBLE_SPACE+3*TCL_INTEGER_SPACE+256];
  IntList p1,p2;
  double r_min=0, r_max=-1.0;
  int r_bins=0, log_flag=0, int_flag=0;
  int i;
  double *distribution, low;

  init_intlist(&p1); init_intlist(&p2);

  if (argc < 2) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze distribution [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_INTLIST(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze distribution [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG1_IS_INTLIST(p2)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze distribution [<type_list> <type_list>]", (char *)NULL);
    return (TCL_ERROR);
  }

  argc -= 2; argv += 2;

  if( argc>0 ) { if (!ARG0_IS_D(r_min)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_D(r_max)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_I(r_bins)) return (TCL_ERROR);   argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_I(log_flag)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_I(int_flag)) return (TCL_ERROR); argc--; argv++; }

  /* if not given set defaults */
  if(r_max == -1.) r_max = min_box_l/2.0;
  if(r_bins < 0 )  r_bins = n_total_particles / 20;

  /* give back what you do */
  Tcl_AppendResult(interp, "{ analyze distribution { ", (char *)NULL);
  for(i=0; i<p1.max; i++) { 
    sprintf(buffer,"%d ",p1.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  Tcl_AppendResult(interp,"} { ", (char *)NULL);
  for(i=0; i<p2.max; i++) {
    sprintf(buffer,"%d ",p2.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  sprintf(buffer,"} %f %f %d %d %d",r_min,r_max,r_bins,log_flag,int_flag);
  Tcl_AppendResult(interp, buffer," }", (char *)NULL);
  /* some sanity checks */
  if(r_min < 0.0 || (log_flag==1 && r_min ==0.0 )) return TCL_ERROR;
  if(r_max <= r_min) return TCL_ERROR;
  if(r_bins < 1) return TCL_ERROR;
  /* calculate distribution */
  distribution = malloc(r_bins*sizeof(double));
  updatePartCfg(WITHOUT_BONDS);
  calc_part_distribution(p1.e, p1.max, p2.e, p2.max, r_min, r_max, r_bins, log_flag,&low,distribution);
  if(int_flag==1) {
    distribution[0] += low;
    for(i=0; i<r_bins-1; i++) distribution[i+1] += distribution[i]; 
  }
  /* append result */
  {
    double log_fac=0.0, bin_width=0.0, r=0.0;
    if(log_flag == 1) {
      log_fac       = pow((r_max/r_min),(1.0/(double)r_bins));
      r = r_min * sqrt(log_fac);
    } 
    else {
      bin_width     = (r_max-r_min) / (double)r_bins;
      r = r_min + bin_width/2.0;
    }
    Tcl_AppendResult(interp, " {\n", (char *)NULL);
    for(i=0; i<r_bins; i++) {
      sprintf(buffer,"%f %f",r,distribution[i]);
      Tcl_AppendResult(interp, "{ ",buffer," }\n", (char *)NULL);
      if(log_flag == 1) r *= log_fac; else r += bin_width;
    }
    Tcl_AppendResult(interp, "}\n", (char *)NULL);
  }
  free(distribution);
  return (TCL_OK);
}

static int parse_vel_distr(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze vel_distr [<type>]' */
  char buffer[3*TCL_DOUBLE_SPACE+3];
  int p1;
  int bins=100;
  double max=0.0;

  /* parse arguments */
  if (argc == 0) {
    Tcl_AppendResult(interp, "usage: analyze vel_distr <type> [bins max]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG0_IS_I(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze vel_distr <type> [bins max]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (argc > 1){
    if (!ARG1_IS_I(bins)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze vel_distr <type> [bins max]", (char *)NULL);
      return (TCL_ERROR);
    }
  }

  if (argc > 2){
    if (!ARG_IS_D(2,max)) {
      Tcl_ResetResult(interp);
      Tcl_AppendResult(interp, "usage: analyze vel_distr <type> [bins max]", (char *)NULL);
      return (TCL_ERROR);
    }
  }

  sprintf(buffer,"%i %i %f",p1,bins,max);
  Tcl_AppendResult(interp, "{ analyze vel_distr ",buffer,"} ",(char *)NULL);
  updatePartCfg(WITHOUT_BONDS);
  calc_vel_distr(interp,p1,bins,max);

  return TCL_OK;
}

static int parse_rdf(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze rdf' (radial distribution function) */
  /************************************************/
  char buffer[2*TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE+256];
  IntList p1,p2;
  double r_min=0, r_max=-1.0;
/* added */
  double x_min=0, x_max=-1.0;
/*Up to here*/  
  int r_bins=-1, n_conf=1, i;
  double *rdf;

  init_intlist(&p1); init_intlist(&p2);

  if (argc < 2) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze {rdf|<rdf>|<rdf-intermol> [<type_list> <type_list>]", (char *)NULL);
    Tcl_AppendResult(interp, "usage: analyze <rdf-adress> [<type_list> <type_list>] [x_min] [x_max]", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG0_IS_INTLIST(p1)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze {rdf|<rdf>|<rdf-intermol>} [<type_list> <type_list>]", (char *)NULL);
    Tcl_AppendResult(interp, "usage: analyze <rdf-adress> [<type_list> <type_list>] [x_min] [x_max]", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG1_IS_INTLIST(p2)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze {rdf|<rdf>|<rdf-intermol>} [<type_list> <type_list>]", (char *)NULL);
    Tcl_AppendResult(interp, "usage: analyze <rdf-adress> [<type_list> <type_list>] [x_min] [x_max]", (char *)NULL);
    return (TCL_ERROR);
  }
  argc-=2; argv+=2;

  if( average==3 ) {
    if( argc>0 ) { if (!ARG0_IS_D(x_min)) return (TCL_ERROR); argc--; argv++; }
    if( argc>0 ) { if (!ARG0_IS_D(x_max)) return (TCL_ERROR); argc--; argv++; }
  }

  if( argc>0 ) { if (!ARG0_IS_D(r_min)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_D(r_max)) return (TCL_ERROR); argc--; argv++; }
  if( argc>0 ) { if (!ARG0_IS_I(r_bins)) return (TCL_ERROR); argc--; argv++; }

  if(average)
    {
      if (n_configs == 0) {
	Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
	Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze rdf' to only look at current RDF!", (char *)NULL);
	return TCL_ERROR;
      }
      if( argc>0 ) {
	if (!ARG0_IS_I(n_conf)) return (TCL_ERROR); argc--; argv++;
      }
      else
        n_conf  = n_configs;
    }

  /* if not given use default */
  if(r_max  == -1.0)  r_max = min_box_l/2.0;
  if(r_bins == -1  )  r_bins = n_total_particles / 20;
/* added */
  if(x_max  == -1.0)  x_max = min_box_l/2.0;
  if(r_bins == -1  )  r_bins = n_total_particles / 20;
/*Up to here*/

  /* give back what you do */
  if(average==0)
    Tcl_AppendResult(interp, "{ analyze rdf { ", (char *)NULL);
  else if(average==1)
    Tcl_AppendResult(interp, "{ analyze <rdf> { ", (char *)NULL);
  else if(average==2)
    Tcl_AppendResult(interp, "{ analyze <rdf-intermol> { ", (char *)NULL);
/* added */
  else if(average==3)
    Tcl_AppendResult(interp, "{ analyze <rdf-adress> { ", (char *)NULL);
/*Up to here*/
  else
    {
      Tcl_AppendResult(interp, "WRONG PARAMETER PASSED ", (char *)NULL);
      return TCL_ERROR;
    }

  for(i=0; i<p1.max; i++) {
    sprintf(buffer,"%d ",p1.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  Tcl_AppendResult(interp,"} { ", (char *)NULL);
  for(i=0; i<p2.max; i++) {
    sprintf(buffer,"%d ",p2.e[i]);
    Tcl_AppendResult(interp, buffer, (char *)NULL);
  }
  sprintf(buffer,"} %f %f %d",r_min,r_max,r_bins);
/* added */
  if(average==3) {
  sprintf(buffer,"} %f %f %f %f %d",x_min,x_max,r_min,r_max,r_bins);
  }
/*Up to here*/
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  if(average) {
    sprintf(buffer," %d",n_conf);
    Tcl_AppendResult(interp, buffer, " }",(char *)NULL);
  }
  else
    Tcl_AppendResult(interp, " }", (char *)NULL);
  rdf = malloc(r_bins*sizeof(double));

  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }

  if(average==0)
    calc_rdf(p1.e, p1.max, p2.e, p2.max, r_min, r_max, r_bins, rdf);
  else if(average==1)
    calc_rdf_av(p1.e, p1.max, p2.e, p2.max, r_min, r_max, r_bins, rdf, n_conf);
  else if(average==2)
    calc_rdf_intermol_av(p1.e, p1.max, p2.e, p2.max, r_min, r_max, r_bins, rdf, n_conf);
/* added */
  else if(average==3)
    calc_rdf_adress(p1.e, p1.max, p2.e, p2.max, x_min, x_max, r_min, r_max, r_bins, rdf, n_conf);
/*Up to here*/

  else ;

  /* append result */
  {
    double bin_width=0.0, r=0.0;
    bin_width     = (r_max-r_min) / (double)r_bins;
    r = r_min + bin_width/2.0;
    Tcl_AppendResult(interp, " {\n", (char *)NULL);
    for(i=0; i<r_bins; i++) {
      sprintf(buffer,"%f %f",r,rdf[i]);
      Tcl_AppendResult(interp, "{ ",buffer," }\n", (char *)NULL);
      r += bin_width;
    }
    Tcl_AppendResult(interp, "}\n", (char *)NULL);
  }
  free(rdf);
  return (TCL_OK);
}


int parse_structurefactor(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze { stucturefactor } <type> <order>' */
  /***********************************************************************************************************/
  char buffer[2*TCL_DOUBLE_SPACE+4];
  int i, type, order;
  double qfak, *sf;
  if (argc < 2) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze structurefactor <type> <order> [<chain_start> <n_chains> <chain_length>]",
		     (char *)NULL);
    return (TCL_ERROR);
  } else {
    if (!ARG0_IS_I(type))
      return (TCL_ERROR);
    if (!ARG1_IS_I(order))
      return (TCL_ERROR);
    argc-=2; argv+=2;
  }
  updatePartCfg(WITHOUT_BONDS);
  calc_structurefactor(type, order, &sf); 
  
  qfak = 2.0*PI/box_l[0];
  for(i=0; i<order*order; i++) { 
    if (sf[2*i+1]> 0) { 
      sprintf(buffer,"{%f %f} ",qfak*sqrt(i+1),sf[2*i]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
  }
  free(sf);
  return (TCL_OK);
}

static int parse_density_profile_av(Tcl_Interp *interp, int argc, char **argv)
{
   /* 'analyze <density_profile> [<n_bin> <density> <dir> <number of conf> <type>]' */
  int n_conf;
  int n_bin;
  double density;
  int dir; 
  double *rho_ave;
  int type;
  int i;
  char buffer[2*TCL_DOUBLE_SPACE+TCL_INTEGER_SPACE+256];
  
  /* parse arguments */
  if (argc < 5) {
    Tcl_AppendResult(interp, "usage: analyze <density_profile> [<n_bin> <density> <dir> <number of conf> <type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if (!ARG0_IS_I(n_bin)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze <density_profile> [<n_bin> <density> <dir> <number of conf> <type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if (!ARG1_IS_D(density)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze <density_profile> [<n_bin> <density> <dir> <number of conf> <type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  argc-=2; argv+=2;

  if( argc>0 ) { if (!ARG0_IS_I(dir)) return (TCL_ERROR); argc--; argv++; }
  if ( argc>0 ) 
     { if (!ARG0_IS_I(n_conf)) return (TCL_ERROR); argc--; argv++; }
  else 
    n_conf  = n_configs;
  
  if (!ARG0_IS_I(type)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze <density_profile> [<n_bin> <density> <dir> <number of conf> <type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  
  rho_ave = malloc(n_bin*sizeof(double));
  for(i=0;i<n_bin;i++)
    rho_ave[i]=0.0;

  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }

  density_profile_av(n_conf, n_bin, density, dir, rho_ave, type);
  /* append result */
  double r_bin, r;
  r_bin = box_l[dir]/(double)(n_bin);
  r=r_bin/2.0;
  Tcl_AppendResult(interp, " {\n", (char *)NULL);
  for(i=0; i<n_bin; i++) {
    sprintf(buffer,"%f %f",r,rho_ave[i]);
    Tcl_AppendResult(interp, "{ ",buffer," }\n", (char *)NULL);
    r += r_bin;
  }
  Tcl_AppendResult(interp, "}\n", (char *)NULL);
  
  free(rho_ave);
  
  return TCL_OK;
}


static int parse_diffusion_profile(Tcl_Interp *interp, int argc, char **argv )
{
  int i;
  int nbins, n_part, n_conf, time, type, dir;
  double xmin, xmax;
  double *bins;  
  char buffer[TCL_DOUBLE_SPACE];
  
  /* parse arguments */
  if (argc < 8) {
    Tcl_AppendResult(interp, "usage: analyze <diffusion_profile> [ <dir> <xmin> <xmax> <nbins> <n_part> <n_conf> <time> <type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG0_IS_I(dir)) {
    Tcl_AppendResult(interp, "usage: analyze <diffusion_profile> [ <dir> <xmin> <xmax> <nbins> <n_part> <n_conf> <time> <type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG1_IS_D(xmin)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze <diffusion_profile> [ <dir> <xmin> <xmax> <nbins> <n_part> <n_conf> <time> <type>]", (char *)NULL);
    return (TCL_ERROR);
  }
  //if(!ARG1_IS_D(xmax)) {
  //Tcl_ResetResult(interp);
  //Tcl_AppendResult(interp, "usage: analyze <diffusion_profile> [<xmin> <xmax> <nbins> <n_part> <n_conf> <time> <type>]", (char *)NULL);
  //return (TCL_ERROR);
  //}
  argc-=2; argv+=2;
  if(argc>0){if(!ARG0_IS_D(xmax)) return (TCL_ERROR); argc--;argv++;}
  if(argc>0){if(!ARG0_IS_I(nbins)) return (TCL_ERROR); argc--;argv++;}
  if(argc>0){if(!ARG0_IS_I(n_part)) return (TCL_ERROR); argc--;argv++;}
  if(argc>0){if(!ARG0_IS_I(n_conf)) return (TCL_ERROR); argc--;argv++;}
  if(argc>0){if(!ARG0_IS_I(time)) return (TCL_ERROR); argc--;argv++;}
  if(argc>0){if(!ARG0_IS_I(type)) return (TCL_ERROR); argc--;argv++;}
  
  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  
  bins = malloc(nbins*sizeof(double));
  for (i =0; i<nbins;i++) { bins[i]=0; }
  
  calc_diffusion_profile(dir, xmin, xmax, nbins, n_part, n_conf, time, type, bins);
  
  double r_bin, r=0;
  r_bin = box_l[0]/(double)(nbins);
  Tcl_AppendResult(interp, " {\n", (char *)NULL);
  for(i=0; i<nbins; i++) {
    sprintf(buffer,"%f %f",r,bins[i]);
    Tcl_AppendResult(interp, "{ ",buffer," }\n", (char *)NULL);
    r += r_bin;
  }
  Tcl_AppendResult(interp, "}\n", (char *)NULL);
  
  free(bins);
  return TCL_OK;
}




static int parse_vanhove(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze vanhove' (van Hove Auto correlation function) */
  /**********************************************************/

  char buffer[2*TCL_DOUBLE_SPACE+4];
  int c,i,ptype=0, rbins=0, np=0, tmax=0;
  double rmin=0, rmax=0;
  double **vanhove=NULL;
  double *msd=NULL;

  /* checks */
  if (argc < 4) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "Wrong # of args! usage: analyze vanhove <part_type> <r_min> <r_max> <r_bins> [<t_max>]", (char *)NULL);
    return (TCL_ERROR);
  }
  
  if (!ARG0_IS_I(ptype)) return (TCL_ERROR); argc--; argv++;
  if (!ARG0_IS_D(rmin))  return (TCL_ERROR); argc--; argv++;
  if (!ARG0_IS_D(rmax))  return (TCL_ERROR); argc--; argv++;
  if (!ARG0_IS_I(rbins)) return (TCL_ERROR); argc--; argv++;
  if (argc==1) {
    if (!ARG0_IS_I(tmax)) return (TCL_ERROR); argc--; argv++;
  } else if (argc>1) {
    Tcl_ResetResult(interp);
    sprintf(buffer, "%d", argc);
    Tcl_AppendResult(interp, "Wrong # of args! usage: analyze vanhove <part_type> <r_min> <r_max> <r_bins> [<t_max>]",(char *)NULL);
    return (TCL_ERROR);
  }

  if (n_configs == 0) {
	Tcl_AppendResult(interp, "analyze vanhove: no configurations found! (This is a dynamic quantity!)", (char *)NULL);
	return TCL_ERROR;
  }

  if (tmax>=n_configs) { 
     Tcl_ResetResult(interp);
     Tcl_AppendResult(interp, "analyze vanhove: setting tmax >= n_configs is not allowed", (char *)NULL);
     return (TCL_ERROR);
  } else if (tmax==0) { tmax=n_configs-1; }

  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }

  /* allocate space */
  vanhove = (double **) malloc((tmax)*sizeof(double *));
  for(c=0; c<(tmax); c++) { 
    vanhove[c] = (double *) malloc(rbins*sizeof(double));
    for(i=0; i<rbins; i++) { vanhove[c][i] = 0; }
  }
  msd = (double *) malloc((tmax)*sizeof(double));
  for(i=0; i<(tmax); i++) { msd[i] = 0; }
 
  /* calculation */
  np = calc_vanhove(ptype,rmin,rmax,rbins,tmax,msd,vanhove);
 
  /* return results */
  if(np==0) {
    Tcl_AppendResult(interp, "{ no particles }", (char *)NULL);
  } else {
    Tcl_AppendResult(interp, "{ msd { ", (char *)NULL);
    for(c=0; c<(tmax); c++) {
      sprintf(buffer,"%f ",msd[c]);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    Tcl_AppendResult(interp, "} } { vanhove { ", (char *)NULL);
    for(c=0; c<(tmax); c++) {
      Tcl_AppendResult(interp, "{ ", (char *)NULL);
      for(i=0; i<rbins; i++) {
	sprintf(buffer,"%f ",vanhove[c][i]);
	Tcl_AppendResult(interp, buffer, (char *)NULL);
      }
      Tcl_AppendResult(interp, "} ", (char *)NULL);
    }
    Tcl_AppendResult(interp, "} } ", (char *)NULL);
  }


  // free space of times and vanhove
  for(c=0; c<(tmax); c++) { free(vanhove[c]); } 
  free(vanhove);
  free(msd);

  if(np>0) { return (TCL_OK); } else { return (TCL_ERROR); }

}


/****************************************************************************************
 *                                 parser for config storage stuff
 ****************************************************************************************/

static int parse_append(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze append' */
  /********************/
  char buffer[2*TCL_INTEGER_SPACE+256];

  if (argc != 0) { Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze append", (char *)NULL); return TCL_ERROR; }
  if (n_total_particles == 0) {
    Tcl_AppendResult(interp,"No particles to append! Use 'part' to create some, or 'analyze configs' to submit a bunch!",(char *) NULL); 
    return (TCL_ERROR); }
  if ((n_configs > 0) && (n_part_conf != n_total_particles)) {
    sprintf(buffer,"All configurations stored must have the same length (previously: %d, now: %d)!", n_part_conf, n_total_particles);
    Tcl_AppendResult(interp,buffer,(char *) NULL); return (TCL_ERROR); 
  }
  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  analyze_append();
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_push(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze push [<size>]' */
  /*****************************/
  char buffer[2*TCL_INTEGER_SPACE+256];
  int i, j;

  if (n_total_particles == 0) {
    Tcl_AppendResult(interp,"No particles to append! Use 'part' to create some, or 'analyze configs' to submit a bunch!",(char *) NULL); 
    return (TCL_ERROR); }
  if ((n_configs > 0) && (n_part_conf != n_total_particles)) {
    sprintf(buffer,"All configurations stored must have the same length (previously: %d, now: %d)!", n_part_conf, n_total_particles);
    Tcl_AppendResult(interp,buffer,(char *) NULL); return (TCL_ERROR); 
  }
  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  if (argc == 1) { 
    if(!ARG0_IS_I(i)) return (TCL_ERROR);
    if (n_configs < i) analyze_append(); else analyze_push();
    if (n_configs > i) for(j=0; j < n_configs-i; j++) analyze_remove(0);
  }
  else if (argc != 0) { Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze push [<size>]", (char *)NULL); return TCL_ERROR; }
  else if (n_configs > 0) analyze_push();
  else analyze_append();
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_replace(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze replace <index>' */
  /*****************************/
  char buffer[2*TCL_INTEGER_SPACE+256];
  int i;

  if (argc != 1) { Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze replace <index>", (char *)NULL); return TCL_ERROR; }
  if (n_total_particles == 0) {
    Tcl_AppendResult(interp,"No particles to append! Use 'part' to create some, or 'analyze configs' to submit a bunch!",(char *) NULL); 
    return (TCL_ERROR); }
  if ((n_configs > 0) && (n_part_conf != n_total_particles)) {
    sprintf(buffer,"All configurations stored must have the same length (previously: %d, now: %d)!", n_part_conf, n_total_particles);
    Tcl_AppendResult(interp,buffer,(char *) NULL); return (TCL_ERROR); 
  }
  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  if (!ARG0_IS_I(i)) return (TCL_ERROR);
  if((n_configs == 0) && (i==0)) analyze_append();
  else if ((n_configs == 0) && (i!=0)) {
    Tcl_AppendResult(interp, "Nice try, but there are no stored configurations that could be replaced!", (char *)NULL); return TCL_ERROR; }
  else if((i < 0) || (i > n_configs-1)) {
    sprintf(buffer,"Index %d out of range (must be in [0,%d])!",i,n_configs-1);
    Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_ERROR; }
  else analyze_replace(i);
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_remove(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze remove [<index>]' */
  /******************************/
  char buffer[2*TCL_INTEGER_SPACE+256];
  int i;

  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  if (argc == 0) { for (i = n_configs-1; i >= 0; i--) analyze_remove(i); }
  else if (argc == 1) {
    if (!ARG0_IS_I(i)) return (TCL_ERROR);
    if(n_configs == 0) {
      Tcl_AppendResult(interp, "Nice try, but there are no stored configurations that could be removed!", (char *)NULL); return TCL_ERROR; }
    else if((i < 0) || (i > n_configs-1)) {
      sprintf(buffer,"Index %d out of range (must be in [0,%d])!",i,n_configs-1);
      Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_ERROR; }
    analyze_remove(i);
  }
  else {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze remove [<index>]", (char *)NULL); return TCL_ERROR; 
  }
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_stored(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze stored' */
  /********************/
  char buffer[TCL_INTEGER_SPACE];
  if (argc != 0) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze stored", (char *)NULL);
    return TCL_ERROR; 
  }
  sprintf(buffer,"%d",n_configs);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return TCL_OK;
}

static int parse_configs(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze configs [ { <which> | <configuration> } ]' */
  /*******************************************************/
  char buffer[3*TCL_DOUBLE_SPACE+4*TCL_INTEGER_SPACE+256];
  double *tmp_config;
  int i, j;

  if (argc == 0) {
    for(i=0; i < n_configs; i++) {
      Tcl_AppendResult(interp,"{ ", (char *)NULL);
      for(j=0; j < n_part_conf; j++) {
	sprintf(buffer,"%f %f %f ",configs[i][3*j],configs[i][3*j+1],configs[i][3*j+2]);
	Tcl_AppendResult(interp, buffer,(char *)NULL);
      }
      Tcl_AppendResult(interp,"} ",(char *)NULL);
    }
    return (TCL_OK); }
  else if (argc == 1) {
    if (!ARG0_IS_I(i)) return (TCL_ERROR);
    if ((i<0) || (i>n_configs-1)) {
      sprintf(buffer,"The configs[%d] you requested does not exist, argument must be in [0,%d]!",i,n_configs-1);
      Tcl_AppendResult(interp,buffer,(char *)NULL); return TCL_ERROR; }
    for(j=0; j < n_part_conf; j++) {
      sprintf(buffer,"%f %f %f ",configs[i][3*j],configs[i][3*j+1],configs[i][3*j+2]);
      Tcl_AppendResult(interp, buffer,(char *)NULL);
    }
    return (TCL_OK); }
  else if ((argc == 3*n_part_conf) || (n_part_conf == 0)) {
    if ((n_part_conf == 0) && (argc % 3 == 0)) n_part_conf = argc/3;
    else if (argc != 3*n_part_conf) {
      sprintf(buffer,"Wrong # of args(%d)! Usage: analyze configs [x0 y0 z0 ... x%d y%d z%d]",argc,n_part_conf,n_part_conf,n_part_conf);
      Tcl_AppendResult(interp,buffer,(char *)NULL); return TCL_ERROR; }
    tmp_config = malloc(3*n_part_conf*sizeof(double));
    for(j=0; j < argc; j++)
      if (!ARG_IS_D(j, tmp_config[j])) return (TCL_ERROR);
    analyze_configs(tmp_config, n_part_conf); free(tmp_config);
    sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
  }
  /* else */
  sprintf(buffer,"Wrong # of args(%d)! Usage: analyze configs [x0 y0 z0 ... x%d y%d z%d]",argc,n_part_conf,n_part_conf,n_part_conf);
  Tcl_AppendResult(interp,buffer,(char *)NULL);
  return TCL_ERROR;
}

static int parse_activate(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze replace <index>' */
  /*****************************/
  char buffer[2*TCL_INTEGER_SPACE+256];
  int i;

  if (argc != 1) { Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze activate <index>", (char *)NULL); return TCL_ERROR; }
  if (n_total_particles == 0) {
    Tcl_AppendResult(interp,"No particles to append! Use 'part' to create some, or 'analyze configs' to submit a bunch!",(char *) NULL); 
    return (TCL_ERROR); }
  if ((n_configs > 0) && (n_part_conf != n_total_particles)) {
    sprintf(buffer,"All configurations stored must have the same length (previously: %d, now: %d)!", n_part_conf, n_total_particles);
    Tcl_AppendResult(interp,buffer,(char *) NULL); return (TCL_ERROR); 
  }
  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }
  if (!ARG0_IS_I(i)) return (TCL_ERROR);
  if((n_configs == 0) && (i==0)) analyze_append();
  else if ((n_configs == 0) && (i!=0)) {
    Tcl_AppendResult(interp, "Nice try, but there are no stored configurations that could be replaced!", (char *)NULL); return TCL_ERROR; }
  else if((i < 0) || (i > n_configs-1)) {
    sprintf(buffer,"Index %d out of range (must be in [0,%d])!",i,n_configs-1);
    Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_ERROR; }
  else analyze_activate(i);
  sprintf(buffer,"%d",n_configs); Tcl_AppendResult(interp, buffer, (char *)NULL); return TCL_OK;
}

static int parse_mol(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze mol <"com" or "force"> <molecule id>' */
  /**************************************************/

  /* returns either the center of mass of a molecule or the force that a trap applies to a molecule */
  /* mol force returns the summed force over all time_steps since it was last called along with the */
  /* number of time steps */

#ifdef MOLFORCES

  int mol, i;
  char buffer[3*TCL_DOUBLE_SPACE+256];

  if (argc < 1){
    Tcl_AppendResult(interp, "Wrong # of args! At least 1 required", (char *)NULL);
    return TCL_ERROR;
  }
  if (ARG0_IS_S("force")) {
    argc -= 1; argv += 1;
    if (argc != 1){
      Tcl_AppendResult(interp, "Wrong # of args! Only mol num is required", (char *)NULL);
      return TCL_ERROR;
    }
    if (!ARG0_IS_I(mol)) return (TCL_ERROR); argc--; argv++;
    if (mol > n_molecules) return (TCL_ERROR);
    sprintf(buffer,"%e %e %e %d", topology[mol].fav[0],topology[mol].fav[1],topology[mol].fav[2],topology[mol].favcounter );
    for (i=0;i<3;i++) {
      topology[mol].fav[i]=0;
    }
    topology[mol].favcounter = 0;
    Tcl_AppendResult(interp,buffer, (char *)NULL);
    return TCL_OK;
  }
  else if (ARG0_IS_S("com")) {
    argc -= 1; argv += 1;
    if (argc != 1){
      Tcl_AppendResult(interp, "Wrong # of args! Only mol num is required", (char *)NULL);
      return TCL_ERROR;
    }
    if (!ARG0_IS_I(mol)) return (TCL_ERROR); argc--; argv++;
    if (mol > n_molecules) return (TCL_ERROR);
    sprintf(buffer,"%e %e %e",topology[mol].com[0],topology[mol].com[1],topology[mol].com[2]);
    Tcl_AppendResult(interp,buffer, (char *)NULL);
    return TCL_OK;
  }
  else {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "The operation \"mol ", argv[1], "\" you requested is not implemented.", (char *)NULL);
    return TCL_ERROR;
  }
#else
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "The operation \"mol\" you requested requires the MOLFORCES option.  Activate it in config.h", (char *)NULL);
    return TCL_ERROR;
#endif
}

static int parse_and_print_momentum(Tcl_Interp *interp, int argc, char **argv)
{
    char buffer[TCL_DOUBLE_SPACE];
    double momentum[3] = { 0., 0., 0. };

    momentum_calc(momentum);

    if (argc == 0) {
      Tcl_PrintDouble(interp, momentum[0], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      Tcl_PrintDouble(interp, momentum[1], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      Tcl_PrintDouble(interp, momentum[2], buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    else if (ARG0_IS_S("particles")) {
      mpi_gather_stats(4, momentum, NULL, NULL, NULL);
      Tcl_PrintDouble(interp, momentum[0], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      Tcl_PrintDouble(interp, momentum[1], buffer);
      Tcl_AppendResult(interp, buffer, " ", (char *)NULL);
      Tcl_PrintDouble(interp, momentum[2], buffer);
      Tcl_AppendResult(interp, buffer, (char *)NULL);
    }
    else {
	Tcl_AppendResult(interp, "unknown feature of: analyze momentum",
			 (char *)NULL);
	return TCL_ERROR;
    }

    return TCL_OK;
}

/** EXPERIMENTAL [uschille] */
static int acf_parse_init(Tcl_Interp *interp, int argc, char **argv, int slot_num) {
  Tcl_AppendResult(interp, "TODO: This is not implemented yet!", (char *)NULL) ;
  return TCL_ERROR ;
}

/** EXPERIMENTAL [uschille] */
static int acf_parse_append(Tcl_Interp *interp, int argc, char **argv, int slot_num) {
  Tcl_AppendResult(interp, "TODO: This is not implemented yet!", (char *)NULL) ;  return TCL_ERROR ;
}

static int acf_print(Tcl_Interp *interp, int slot_num) {
  Tcl_AppendResult(interp, "TODO: This is not implemented yet!", (char *)NULL) ;  return TCL_ERROR ;
}

/** EXPERIMENTAL [uschille] */
int acf_cmd(ClientData data, Tcl_Interp *interp, int argc, char ** argv) {

  int err = TCL_OK ;
  int slot_num = 0 ;

  if (argc < 2) {
    Tcl_AppendResult(interp, "Wrong # of args!", (char *)NULL) ;
    return TCL_ERROR ;
  }
    
  if (!ARG1_IS_I(slot_num)) {
    Tcl_AppendResult(interp, "Error while parsing arg 1 of acf!", (char *)NULL) ;
    return TCL_ERROR ;
  }

  if (slot_num < 0) {
    Tcl_AppendResult(interp, "Error: slot number must be positive", (char *)NULL) ;
    return TCL_ERROR ;
  }
  
  if (argc == 2) {
    err = acf_print(interp, slot_num) ;
  }
  else if (ARG_IS_S(2,"init")) {
    err = acf_parse_init(interp, argc-3, argv+3, slot_num) ;
  }
  else if (ARG_IS_S(2,"append")) {
    err = acf_parse_append(interp, argc-3, argv+3, slot_num) ;
  }
  else {
    Tcl_AppendResult(interp, "Error: unknown acf instruction \"",argv[2],"\"", (char *)NULL) ;
    err = TCL_ERROR ;
  }

  if (err==TCL_ERROR) {
    Tcl_AppendResult(interp, "Usage:", (char *)NULL) ;
  }

  return err ;

}

void centermass_conf(int k, int type_1, double *com)
{
  int i, j;
  double M = 0.0;
  com[0]=com[1]=com[2]=0.;

  for (j=0; j<n_total_particles; j++) {
    if ((partCfg[j].p.type == type_1) || (type_1 == -1))
    {
      for (i=0; i<3; i++)
      {
         com[i] += configs[k][3*j+i]*PMASS(partCfg[j]);
      }
      M += PMASS(partCfg[j]);
    }
  }
  for (i=0; i<3; i++) 
  {
    com[i] /= M;
  }
  return;
}

double calc_diffusion_coef(Tcl_Interp *interp,int type_m, int n_time_steps,int n_conf)
{
  int i,j,k;
  double  p1[3],p2[3],p_com[3],p_x[n_configs],p_y[n_configs],p_z[n_configs];
  double MSD[n_configs],MSD_time;
  double D;
  char buffer[TCL_DOUBLE_SPACE];
  int MSD_particles=0;
  int start_value=n_configs-n_conf;

  if (!sortPartCfg()) { Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",(char *) NULL); return (TCL_ERROR); }

  for(i=start_value;i<n_configs;i++)
  {
      MSD[i]=0.0;//MSD for all saved confs
      centermass_conf(i, type_m, p_com); //COM for all saved confs
      p_x[i]=p_com[0];
      p_y[i]=p_com[1];
      p_z[i]=p_com[2];
  }

  for(j=0; j<n_total_particles; j++) {
     if((partCfg[j].p.type == type_m)||(type_m == -1)) {
       MSD_particles++;//count particles for MSD
       for(i=start_value;i<n_configs;i++) {
          p1[0]=configs[i][3*j  ]-p_x[i];
          p1[1]=configs[i][3*j+1]-p_y[i];
          p1[2]=configs[i][3*j+2]-p_z[i];
          for (k=i;k<n_configs;k++)
          {
             p2[0]=configs[k][3*j  ]-p_x[k];
             p2[1]=configs[k][3*j+1]-p_y[k];
             p2[2]=configs[k][3*j+2]-p_z[k];
             MSD[k-i]+=distance2(p1, p2);
          }
        }
     }
  }

 // normalization
  if (MSD_particles!=0) 
  {
      //average over all com particles and time origins
      for (i=start_value;i<n_configs;i++)
      {
          MSD[i]/=(double)(MSD_particles*(n_configs-i));
          MSD_time=time_step*n_time_steps*(i-start_value);
          sprintf(buffer,"{ %e %e }",MSD_time,MSD[i]);
          Tcl_AppendResult(interp,buffer,"\n",(char *)NULL);
      }
      MSD_time=time_step*n_time_steps*(n_configs-1-start_value);
      D=(MSD[n_configs-1]-MSD[start_value])/(6.0*MSD_time);
  }
  else
  {
      D=0;
  }
  return D;
}

static int parse_and_print_dipole(Tcl_Interp *interp,int argc, char **argv)
{
   int i,k;
   char buffer[TCL_DOUBLE_SPACE];
   double dipole[3],total_q=0.0;
   updatePartCfg(WITHOUT_BONDS);
   if (!sortPartCfg()) {
      char *errtxt = runtime_error(128);
      ERROR_SPRINTF(errtxt, "{059 parse_and_print_dipole: could not sort particle config, particle ids not consecutive?} ");
      return TCL_ERROR;
   }
   for (i=0;i<3;i++)
   {
       dipole[i]=0;
   }
   for (i=0;i<n_total_particles;i++)
   {
#ifdef ELECTROSTATICS
       total_q+=partCfg[i].p.q;
       for (k=0;k<3;k++){
            dipole[k]+=partCfg[i].r.p[k]*partCfg[i].p.q;
       }
#endif       
   }
   Tcl_AppendResult(interp,"{ dipolemoment_normal ",(char *)NULL);
   for (k=0;k<3;k++)
   {
       sprintf(buffer,"%e ",dipole[k]);
       Tcl_AppendResult(interp, buffer,(char *)NULL);
   }
   sprintf(buffer,"%e",total_q);
   Tcl_AppendResult(interp,buffer,"}",(char *)NULL);
   return TCL_OK;
}

static int parse_MSD(Tcl_Interp *interp, int argc, char **argv)
{
  /* 'analyze MSD [ <type_m> <n_time_steps>]' */
  int n_time_steps;
  int type_m;
  int n_conf;
  char buffer[3*TCL_DOUBLE_SPACE];
  double D;
  
  /* parse arguments */
  if (argc < 2) {
    Tcl_AppendResult(interp, "usage: analyze MSD {<type_m> <n_time_steps>} [<number of conf>]", (char *)NULL);
    return (TCL_ERROR);
  }


  if (!ARG0_IS_I(type_m)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze MSD {<type_m> <n_time_steps>} [<number of conf>]", (char *)NULL);
    return (TCL_ERROR);
  }

  if (!ARG1_IS_I(n_time_steps)) {
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "usage: analyze MSD {<type_m> <n_time_steps>} [<number of conf>]", (char *)NULL);
    return (TCL_ERROR);
  }
  argc-=2; argv+=2;

  if (n_configs == 0) 
  {
    Tcl_AppendResult(interp, "No configurations found! ", (char *)NULL);
    Tcl_AppendResult(interp, "Use 'analyze append' to save some !", (char *)NULL);
    return TCL_ERROR;
  }
  if( argc>0 ) {
    if (!ARG0_IS_I(n_conf)) return (TCL_ERROR);
    argc--;
    argv++;
  }
  else
  {
    n_conf  = n_configs;
  }
  
  sprintf(buffer,"%i %i %i",type_m,n_time_steps,n_configs);
  Tcl_AppendResult(interp, "{ analyze MSD ",buffer," } {\n",(char *)NULL);
  D=calc_diffusion_coef(interp,type_m, n_time_steps,n_conf);
  sprintf(buffer,"%e",D);
  Tcl_AppendResult(interp, "}\n{approx. D=",buffer,"}", (char *)NULL);
  return TCL_OK;
}

static int parse_and_print_energy_kinetic(Tcl_Interp *interp,int argc, char **argv)
{
   int i,type;
   char buffer[TCL_DOUBLE_SPACE];
   double E_kin=0;

  /* parse arguments */
  if (argc < 1) {
    Tcl_AppendResult(interp, "usage: analyze energy_kinetic <type>", (char *)NULL);
    return (TCL_ERROR);
  }
  if (!ARG0_IS_I(type)) {
     Tcl_AppendResult(interp, "usage: analyze energy_kinetic <type> where type is int", (char *)NULL);
     return (TCL_ERROR);
  }
  updatePartCfg(WITHOUT_BONDS);
  for (i=0;i<n_total_particles;i++)
  {
      if (partCfg[i].p.type == type )
      {
         E_kin+=PMASS(partCfg[i])*sqrlen(partCfg[i].m.v);
      }
  }
  E_kin*=0.5/time_step/time_step;
  Tcl_PrintDouble(interp, E_kin, buffer);;
  Tcl_AppendResult(interp, buffer,(char *)NULL);
  return TCL_OK;
}

/****************************************************************************************
 *                                 main parser for analyze
 ****************************************************************************************/

int analyze(ClientData data, Tcl_Interp *interp, int argc, char **argv)
{
  int err = TCL_OK;
  if (argc < 2) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze <what> ...", (char *)NULL);
    return (TCL_ERROR);
  }

  /* for general options */
#define REGISTER_ANALYZE_OPTION(name, parser)				\
  else if (ARG1_IS_S(name)) err = parser(interp, argc - 2, argv + 2)

  /* for commands of the config storage */
#define REGISTER_ANALYZE_STORAGE(name, parser) \
  else if (ARG1_IS_S(name)) err = parser(interp, argc - 2, argv + 2)

  /* for actual observables */
#define REGISTER_ANALYSIS(name, parser)					\
  else if (ARG1_IS_S(name)) err = parser(interp, argc - 2, argv + 2)
#define REGISTER_ANALYSIS_W_ARG(name, parser, arg)			\
  else if (ARG1_IS_S(name)) err = parser(interp, arg, argc - 2, argv + 2)

  /* for the elses below */
  if (0);

  REGISTER_ANALYZE_OPTION("set", parse_analyze_set_topology);
#ifdef LB
  REGISTER_ANALYZE_OPTION("fluid", parse_analyze_fluid);
#endif
  REGISTER_ANALYSIS("get_folded_positions", parse_get_folded_positions);
#ifdef MODES
  REGISTER_ANALYZE_OPTION("set_bilayer", parse_bilayer_set);
  REGISTER_ANALYSIS("modes2d", parse_modes2d);
  REGISTER_ANALYSIS("bilayer_density_profile", parse_bilayer_density_profile);
  REGISTER_ANALYSIS("radial_density_map", parse_radial_density_map);
  REGISTER_ANALYSIS("get_lipid_orients", parse_get_lipid_orients);
  REGISTER_ANALYSIS("lipid_orient_order", parse_lipid_orient_order);
#endif
  REGISTER_ANALYSIS("mol", parse_mol);
  REGISTER_ANALYSIS("cluster_size_dist", parse_cluster_size_dist);
  REGISTER_ANALYSIS("mindist", parse_mindist);
  REGISTER_ANALYSIS("aggregation", parse_aggregation);
  REGISTER_ANALYSIS("centermass", parse_centermass);
  REGISTER_ANALYSIS("angularmomentum",parse_angularmomentum);
  REGISTER_ANALYSIS("MSD",parse_MSD);
  REGISTER_ANALYSIS("dipmom_normal",parse_and_print_dipole);
  REGISTER_ANALYSIS("momentofinertiamatrix", parse_momentofinertiamatrix);
  REGISTER_ANALYSIS("find_principal_axis", parse_find_principal_axis);
  REGISTER_ANALYSIS("nbhood", parse_nbhood);
  REGISTER_ANALYSIS("distto", parse_distto);
  REGISTER_ANALYSIS("cell_gpb", parse_cell_gpb);
  REGISTER_ANALYSIS("Vkappa", parse_Vkappa);
  REGISTER_ANALYSIS("energy", parse_and_print_energy);
  REGISTER_ANALYSIS("energy_kinetic", parse_and_print_energy_kinetic);
  REGISTER_ANALYSIS_W_ARG("pressure", parse_and_print_pressure, 0);
#ifdef VIRTUAL_SITES
  REGISTER_ANALYSIS("energy_kinetic_mol",parse_and_print_energy_kinetic_mol);
  REGISTER_ANALYSIS("pressure_mol",  parse_and_print_pressure_mol);
  REGISTER_ANALYSIS("check_mol",  parse_and_check_mol_pos);
  REGISTER_ANALYSIS("dipmom_mol",  parse_and_print_dipole_mol);
#endif
  REGISTER_ANALYSIS_W_ARG("stress_tensor", parse_and_print_stress_tensor, 0);
  REGISTER_ANALYSIS("local_stress_tensor", parse_local_stress_tensor);
  REGISTER_ANALYSIS_W_ARG("p_inst", parse_and_print_pressure, 1);
  REGISTER_ANALYSIS("momentum", parse_and_print_momentum);
  REGISTER_ANALYSIS("bins", parse_bins);
  REGISTER_ANALYSIS("p_IK1", parse_and_print_p_IK1);
  REGISTER_ANALYSIS_W_ARG("re", parse_re, 0);
  REGISTER_ANALYSIS_W_ARG("<re>", parse_re, 1);
  REGISTER_ANALYSIS_W_ARG("rg", parse_rg, 0);
  REGISTER_ANALYSIS_W_ARG("<rg>", parse_rg, 1);
  REGISTER_ANALYSIS_W_ARG("rh", parse_rh, 0);
  REGISTER_ANALYSIS_W_ARG("<rh>", parse_rh, 1);
  REGISTER_ANALYSIS_W_ARG("internal_dist", parse_intdist, 0);
  REGISTER_ANALYSIS_W_ARG("<internal_dist>", parse_intdist, 1);
  REGISTER_ANALYSIS_W_ARG("bond_l", parse_bond_l, 0);
  REGISTER_ANALYSIS_W_ARG("<bond_l>", parse_bond_l, 1);
  REGISTER_ANALYSIS_W_ARG("bond_dist", parse_bond_dist, 0);
  REGISTER_ANALYSIS_W_ARG("<bond_dist>", parse_bond_dist, 1);
  REGISTER_ANALYSIS_W_ARG("g123", parse_g123, 1);    
  REGISTER_ANALYSIS_W_ARG("<g1>", parse_g_av, 1);    
  REGISTER_ANALYSIS_W_ARG("<g2>", parse_g_av, 2);    
  REGISTER_ANALYSIS_W_ARG("<g3>", parse_g_av, 3);
  REGISTER_ANALYSIS_W_ARG("formfactor", parse_formfactor, 0);
  REGISTER_ANALYSIS_W_ARG("<formfactor>", parse_formfactor, 1);    
  REGISTER_ANALYSIS("necklace", parse_necklace_analyzation);   
  REGISTER_ANALYSIS("holes", parse_hole_cluster_analyzation);   
  REGISTER_ANALYSIS("distribution", parse_distribution);
  REGISTER_ANALYSIS("vel_distr", parse_vel_distr);
  REGISTER_ANALYSIS_W_ARG("rdf", parse_rdf, 0);
  REGISTER_ANALYSIS_W_ARG("<rdf>", parse_rdf, 1);
  REGISTER_ANALYSIS_W_ARG("<rdf-intermol>", parse_rdf, 2);
  REGISTER_ANALYSIS_W_ARG("<rdf-adress>", parse_rdf, 3);
  REGISTER_ANALYSIS("rdfchain", parse_rdfchain);
#ifdef ELECTROSTATICS
  REGISTER_ANALYSIS("cwvac", parse_cwvac);
#endif
  REGISTER_ANALYSIS("structurefactor", parse_structurefactor);
  REGISTER_ANALYSIS("<density_profile>", parse_density_profile_av);
  REGISTER_ANALYSIS("<diffusion_profile>", parse_diffusion_profile);
  REGISTER_ANALYSIS("vanhove", parse_vanhove);
  REGISTER_ANALYZE_STORAGE("append", parse_append);
  REGISTER_ANALYZE_STORAGE("push", parse_push);
  REGISTER_ANALYZE_STORAGE("replace", parse_replace);
  REGISTER_ANALYZE_STORAGE("activate", parse_activate);
  REGISTER_ANALYZE_STORAGE("remove", parse_remove);
  REGISTER_ANALYZE_STORAGE("stored", parse_stored);
  REGISTER_ANALYZE_STORAGE("configs", parse_configs);
  else {
    /* the default */
    /***************/
    Tcl_ResetResult(interp);
    Tcl_AppendResult(interp, "The operation \"", argv[1],
		     "\" you requested is not implemented.", (char *)NULL);
    err = (TCL_ERROR);
  }
  return mpi_gather_runtime_errors(interp, err);
}

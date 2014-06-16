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
/** \file statistics.cpp
    This is the place for analysis (so far...).
    Implementation of statistics.hpp
*/
#include <cstdlib>
#include <cstring>
#include "utils.hpp"
#include "statistics.hpp"
#include "statistics_chain.hpp"
#include "statistics_molecule.hpp"
#include "statistics_cluster.hpp"
#include "statistics_fluid.hpp"
//#include "statistics_correlation.hpp"
#include "energy.hpp"
#include "modes.hpp"
#include "pressure.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "particle_data.hpp"
#include "interaction_data.hpp"
#include "domain_decomposition.hpp"
#include "verlet.hpp"
#include "lb.hpp"
#include "virtual_sites.hpp"
#include "initialize.hpp"

/** Previous particle configurations (needed for offline analysis and
    correlation analysis in \ref tclcommand_analyze) */
double **configs = NULL; int n_configs = 0; int n_part_conf = 0;

/****************************************************************************************
 *                                 helper functions
 ****************************************************************************************/

double min_distance2(double pos1[3], double pos2[3])
{
  double diff[3];
  get_mi_vector(diff, pos1, pos2);
  return sqrlen(diff);
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
  for (j=0; j<n_part-1; j++) {
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

    for (i=j+1; i<n_part; i++)
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

  MPI_Reduce(momentum, result, 3, MPI_DOUBLE, MPI_SUM, 0, comm_cart);
}

/** Calculate total momentum of the system (particles & LB fluid)
 * @param momentum Rsult for this processor (Output)
 */
void momentum_calc(double *momentum) 
{
    double momentum_fluid[3] = { 0., 0., 0. };
    double momentum_particles[3] = { 0., 0., 0. };

    mpi_gather_stats(4, momentum_particles, NULL, NULL, NULL);
#ifdef LB
    mpi_gather_stats(6, momentum_fluid, NULL, NULL, NULL);
#endif

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
  for (j=0; j<n_part; j++) {
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
  for (j=0; j<n_part; j++) {
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
  for (j=0; j<n_part; j++) 
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
  for (j=0; j<n_part; j++) {
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

void calc_gyration_tensor(int type, double **_gt)
{
  int i, j, count;
  double com[3];
  double eva[3],eve0[3],eve1[3],eve2[3];
  double *gt=NULL, tmp;
  double Smatrix[9],p1[3];

  for (i=0; i<9; i++) Smatrix[i] = 0;
  /* 3*ev, rg, b, c, kappa, eve0[3], eve1[3], eve2[3]*/
  *_gt = gt = (double*)realloc(gt,16*sizeof(double)); 

  updatePartCfg(WITHOUT_BONDS);

  /* Calculate the position of COM */
  centermass(type,com);

  /* Calculate the gyration tensor Smatrix */
  count=0;
  for (i=0;i<n_part;i++) {
    if ((partCfg[i].p.type == type) || (type == -1)) {
      for ( j=0; j<3 ; j++ ) { 
        p1[j] = partCfg[i].r.p[j] - com[j];
      }
      count ++;
      Smatrix[0] += p1[0]*p1[0];
      Smatrix[1] += p1[0]*p1[1];
      Smatrix[2] += p1[0]*p1[2];
      Smatrix[4] += p1[1]*p1[1];
      Smatrix[5] += p1[1]*p1[2];
      Smatrix[8] += p1[2]*p1[2];
    }
  }
  /* use symmetry */
  Smatrix[3]=Smatrix[1];
  Smatrix[6]=Smatrix[2];
  Smatrix[7]=Smatrix[5];
  for (i=0;i<9;i++){
    Smatrix[i] /= count;
  }

  /* Calculate the eigenvalues of Smatrix */
  i=calc_eigenvalues_3x3(Smatrix, eva);
  tmp=0.0;
  for (i=0;i<3;i++) {
    /* Eigenvalues */
    gt[i] = eva[i];
    tmp += eva[i];
  }
  
  i=calc_eigenvector_3x3(Smatrix,eva[0],eve0); 
  i=calc_eigenvector_3x3(Smatrix,eva[1],eve1); 
  i=calc_eigenvector_3x3(Smatrix,eva[2],eve2); 
  gt[3] = tmp; /* Squared Radius of Gyration */
  gt[4] = eva[0]-0.5*(eva[1]+eva[2]);  /* Asphericity */
  gt[5] = eva[1]-eva[2];  /* Acylindricity */
  gt[6] = (gt[4]*gt[4]+0.75*gt[5]*gt[5])/(gt[3]*gt[3]); /* Relative shape anisotropy */
  /* Eigenvectors */
  for (j=0;j<3;j++) {
    gt[7+j]=eve0[j]; 
    gt[10+j]=eve1[j];
    gt[13+j]=eve2[j];
  }
}


void nbhood(double pt[3], double r, IntList *il, int planedims[3] )
{
  double d[3];
  int i,j;
  double r2;

  r2 = r*r;

  init_intlist(il);
 
  updatePartCfg(WITHOUT_BONDS);

  for (i = 0; i<n_part; i++) {
    if ( (planedims[0] + planedims[1] + planedims[2]) == 3 ) {
      get_mi_vector(d, pt, partCfg[i].r.p);
    } else {
      /* Calculate the in plane distance */
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

  updatePartCfg(WITHOUT_BONDS);

  /* larger than possible */
  mindist=SQR(box_l[0] + box_l[1] + box_l[2]);
  for (i=0; i<n_part; i++) {
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
      fprintf(stderr,"WARNING: Lower boundary is actually larger than l.hpp.s, flipping!\n");
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
  for(i=0; i<n_part; i++) {
    for(t1=0; t1<n_p1; t1++) {
      if(partCfg[i].p.type == p1_types[t1]) {
	min_dist2 = start_dist2;
	/* particle loop: p2_types*/
	for(j=0; j<n_part; j++) {
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


void calc_rdf(int *p1_types, int n_p1, int *p2_types, int n_p2, 
	      double r_min, double r_max, int r_bins, double *rdf)
{
  long int cnt=0;
  int i,j,t1,t2,ind;
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
  for(i=0; i<n_part; i++) {
    for(t1=0; t1<n_p1; t1++) {
      if(partCfg[i].p.type == p1_types[t1]) {
	/* distinguish mixed and identical rdf's */
	if(mixed_flag == 1) start = 0;
	else                start = (i+1);
	/* particle loop: p2_types*/
	for(j=start; j<n_part; j++) {
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
  long int cnt=0;
  int i,j,k,l,t1,t2,ind,cnt_conf=1;
  int mixed_flag=0,start;
  double inv_bin_width=0.0,bin_width=0.0, dist;
  double volume, bin_volume, r_in, r_out;
  double *rdf_tmp, p1[3],p2[3];

  rdf_tmp = (double*)malloc(r_bins*sizeof(double));

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
    for(i=0; i<n_part; i++) {
      for(t1=0; t1<n_p1; t1++) {
	if(partCfg[i].p.type == p1_types[t1]) {
	  // distinguish mixed and identical rdf's
	  if(mixed_flag == 1) start = 0;
	  else                start = (i+1);
	  //particle loop: p2_types
	  for(j=start; j<n_part; j++) {
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

  rdf_tmp = (double*)malloc(r_bins*sizeof(double));

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
    for(i=0; i<n_part; i++) {
      for(t1=0; t1<n_p1; t1++) {
	if(partCfg[i].p.type == p1_types[t1]) {
	  // distinguish mixed and identical rdf's
	  if(mixed_flag == 1) start = 0;
	  else                start = (i+1);
	  //particle loop: p2_types
	  for(j=start; j<n_part; j++) {
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

void calc_structurefactor(int type, int order, double **_ff) {
  int i, j, k, n, qi, p, order2;
  double qr, twoPI_L, C_sum, S_sum, *ff=NULL;
  
  order2 = order*order;
  *_ff = ff = (double*)realloc(ff,2*order2*sizeof(double));
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
	    for(p=0; p<n_part; p++) {
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
    for(p=0; p<n_part; p++) {
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
      for(i=0; i<n_part; i++) {
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
  label = (int*)malloc(n_part*sizeof(int));
  
  /* calculation over last n_conf configurations */
  t=n_configs-n_conf;
  
  while (t<n_configs-time) {
    /* check initial condition */
    count = 0;
    
    for (i=0;i<n_part;i++) {
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
  for (i = 0 ; i < n_part ; i++) {
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,0);
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,1);
    fold_coordinate(partCfg[i].r.p,partCfg[i].l.i,2);
  }

  beadcount = 0;
  xav = 0.0;
  yav = 0.0;
  for ( pi = 0 ; pi < n_part ; pi++ ) {
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
    thetaradii = (double*)malloc(thetabins*nbeadtypes*sizeof(double));
    thetacounts = (int*)malloc(thetabins*nbeadtypes*sizeof(int));
    for ( bi = 0 ; bi < nbeadtypes ; bi++ ) {
      for ( t = 0 ; t < thetabins ; t++ ) {
	thetaradii[bi*thetabins+t] = 0.0;
	thetacounts[bi*thetabins+t] = 0.0;
      }
    }
    /* Maybe there is a nicer way to do this but now I will just repeat the loop over all particles */
      for ( pi = 0 ; pi < n_part ; pi++ ) {
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
  return ES_OK;
}

double calc_vanhove(int ptype, double rmin, double rmax, int rbins, int tmax, double *msd, double **vanhove) 
{
  int i, c1, c3, c3_max, ind, np=0;
  double p1[3],p2[3],dist;
  double bin_width, inv_bin_width;
  IntList p;

  /* create particle list */
  init_intlist(&p);
  for(i=0; i<n_part; i++) { if( partCfg[i].p.type == ptype ) { np ++; } }
  if(np==0) { return 0; }
  alloc_intlist(&p,np);
  for(i=0; i<n_part; i++) { if( partCfg[i].p.type == ptype ) { p.e[p.n]=i; p.n++; } }

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
  n_part_conf = n_part;
  configs = (double**)realloc(configs,(n_configs+1)*sizeof(double *));
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
  n_part_conf = n_part;
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
  n_part_conf = n_part;
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
  configs = (double**)realloc(configs,n_configs*sizeof(double *));
  if (n_configs == 0) n_part_conf = 0;
}

void analyze_configs(double *tmp_config, int count) {
  int i;
  n_part_conf = count;
  configs = (double**)realloc(configs,(n_configs+1)*sizeof(double *));
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
  n_part_conf = n_part;

  for(i=0; i<n_part_conf; i++) {
    pos[0] = configs[ind][3*i];
    pos[1] = configs[ind][3*i+1];
    pos[2] = configs[ind][3*i+2];
    if (place_particle(i, pos)==ES_ERROR) {
      char *errtxt = runtime_error(128 + ES_INTEGER_SPACE);
      ERROR_SPRINTF(errtxt, "{057 failed upon replacing particle %d in Espresso} ", i); 
    }
  }
}


/****************************************************************************************
 *                                 Observables handling
 ****************************************************************************************/

void obsstat_realloc_and_clear(Observable_stat *stat, int n_pre, int n_bonded, int n_non_bonded,
			       int n_coulomb, int n_dipolar, int n_vsr, int c_size)
{
  
  int i;
  // Number of doubles to store pressure in
  int total = c_size*(n_pre + n_bonded_ia + n_non_bonded + n_coulomb + n_dipolar + n_vsr);

  // Allocate mem for the double list
  realloc_doublelist(&(stat->data), stat->data.n = total);
  // Number of doubles per interaction (pressure=1, stress tensor=9,...)
  stat->chunk_size = c_size;
  
  // Number of chunks for different interaction types
  stat->n_coulomb    = n_coulomb;
  stat->n_dipolar    = n_dipolar;
  stat->n_non_bonded = n_non_bonded;
  stat->n_vs_relative = n_vsr; // virtual sites relative (rigid bodies)
  // Pointers to the start of different contributions
  stat->bonded     = stat->data.e + c_size*n_pre;
  stat->non_bonded = stat->bonded + c_size*n_bonded_ia;
  stat->coulomb    = stat->non_bonded + c_size*n_non_bonded;
  stat->dipolar    = stat->coulomb    + c_size*n_coulomb;
  stat->vs_relative    = stat->dipolar    + c_size*n_dipolar;

  // Set all obseravables to zero
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


  //subfunction: mark all neighbors of a particle and their neighbors (recursiv!)
void mark_neighbours(int type,int pa_nr,double dist,int *list){
  int k;
  for (k=0;k<n_part;k++){
     //only unmarked and particles with right distance
     if ( (partCfg[k].p.type == type) && (list[k] == 0) && (min_distance(partCfg[pa_nr].r.p,partCfg[k].r.p) < dist) ){
        //mark particle with same number as calling particle
        list[k]=list[pa_nr];
        mark_neighbours(type,k,dist,list);
     }
  }
}


void centermass_conf(int k, int type_1, double *com)
{
  int i, j;
  double M = 0.0;
  com[0]=com[1]=com[2]=0.;

  for (j=0; j<n_part; j++) {
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

// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
/** \file statistics_chain.c
    Implementation of \ref statistics_chain.h "statistics_chain.h".
*/
#include "statistics.h"
#include "parser.h"
#include "debug.h"
#include "topology.h"
#include "communication.h"
#include "cells.h"

/** Particles' initial positions (needed for g1(t), g2(t), g3(t) in \ref #analyze) */
/*@{*/
float *partCoord_g=NULL, *partCM_g=NULL;
int n_part_g = 0, n_chains_g = 0;
/*@}*/

/** data for a system consisting of chains. TBRS. */
/*@{*/
int chain_start = 0;
int chain_n_chains = 0;
int chain_length = 0;
/*@}*/

void calc_re(double **_re)
{
  int i;
  double dx,dy,dz;
  double dist=0.0,dist2=0.0,dist4=0.0;
  double *re=NULL, tmp;
  *_re = re = realloc(re,4*sizeof(double));

  for (i=0; i<chain_n_chains; i++) {
    dx = partCfg[chain_start+i*chain_length + chain_length-1].r.p[0]
      - partCfg[chain_start+i*chain_length].r.p[0];
    dy = partCfg[chain_start+i*chain_length + chain_length-1].r.p[1]
      - partCfg[chain_start+i*chain_length].r.p[1];
    dz = partCfg[chain_start+i*chain_length + chain_length-1].r.p[2]
      - partCfg[chain_start+i*chain_length].r.p[2];
    tmp = (SQR(dx) + SQR(dy) + SQR(dz));
    dist  += sqrt(tmp);
    dist2 += tmp;
    dist4 += tmp*tmp;
  }
  tmp = (double)chain_n_chains;
  re[0] = dist/tmp;
  re[2] = dist2/tmp;
  re[1] = sqrt(re[2] - re[0]*re[0]);
  re[3] = sqrt(dist4/tmp - re[2]*re[2]);
}

void calc_re_av(double **_re)
{
  int i,j;
  double dx,dy,dz;
  double dist=0.0,dist2=0.0,dist4=0.0;
  double *re=NULL, tmp;
  *_re = re = realloc(re,4*sizeof(double));

  for (j=0; j<n_configs; j++) {
    for (i=0; i<chain_n_chains; i++) {
      dx = configs[j][3*(chain_start+i*chain_length + chain_length-1)]     
	- configs[j][3*(chain_start+i*chain_length)];
      dy = configs[j][3*(chain_start+i*chain_length + chain_length-1) + 1] 
	- configs[j][3*(chain_start+i*chain_length) + 1];
      dz = configs[j][3*(chain_start+i*chain_length + chain_length-1) + 2] 
	- configs[j][3*(chain_start+i*chain_length) + 2];
      tmp = (SQR(dx) + SQR(dy) + SQR(dz));
      dist  += sqrt(tmp);
      dist2 += tmp;
      dist4 += tmp*tmp;
    }
  }
  tmp = (double)chain_n_chains*n_configs;
  re[0] = dist/tmp;
  re[2] = dist2/tmp;
  re[1] = sqrt(re[2] - re[0]*re[0]);
  re[3] = sqrt(dist4/tmp - re[2]*re[2]);
}

void calc_rg(double **_rg)
{
  int i, j;
  double dx,dy,dz, r_CM_x,r_CM_y,r_CM_z;
  double r_G=0.0,r_G2=0.0,r_G4=0.0;
  double *rg=NULL, IdoubMPC, tmp;
  *_rg = rg = realloc(rg,4*sizeof(double));

  IdoubMPC = 1./(double)chain_length;
  for (i=0; i<chain_n_chains; i++) {
    r_CM_x = r_CM_y = r_CM_z = 0.0;
    for (j=0; j<chain_length; j++) {
      r_CM_x += partCfg[chain_start+i*chain_length + j].r.p[0];
      r_CM_y += partCfg[chain_start+i*chain_length + j].r.p[1];
      r_CM_z += partCfg[chain_start+i*chain_length + j].r.p[2];
    }
    r_CM_x *= IdoubMPC;
    r_CM_y *= IdoubMPC;
    r_CM_z *= IdoubMPC;
    tmp = 0.0;
    for (j=0; j<chain_length; ++j) {
      dx = partCfg[chain_start+i*chain_length + j].r.p[0] - r_CM_x;
      dy = partCfg[chain_start+i*chain_length + j].r.p[1] - r_CM_y;
      dz = partCfg[chain_start+i*chain_length + j].r.p[2] - r_CM_z;
      tmp += (SQR(dx) + SQR(dy) + SQR(dz));
    }
    tmp *= IdoubMPC;
    r_G  += sqrt(tmp);
    r_G2 += tmp;
    r_G4 += tmp*tmp;
  }
  tmp = (double)chain_n_chains;
  rg[0] = r_G/tmp;
  rg[2] = r_G2/tmp;
  rg[1] = sqrt(rg[2] - rg[0]*rg[0]);
  rg[3] = sqrt(r_G4/tmp - rg[2]*rg[2]);
}

void calc_rg_av(double **_rg)
{
  int i, j, k;
  double dx,dy,dz, r_CM_x,r_CM_y,r_CM_z;
  double r_G=0.0,r_G2=0.0, r_G4=0.0;
  double *rg=NULL, IdoubMPC, tmp;
  *_rg = rg = realloc(rg,4*sizeof(double));

  IdoubMPC = 1./(double)chain_length;
  for (k=0; k<n_configs; k++) {
    for (i=0; i<chain_n_chains; i++) {
      r_CM_x = r_CM_y = r_CM_z = 0.0;
      for (j=0; j<chain_length; j++) {
	r_CM_x += configs[k][3*(chain_start+i*chain_length + j)];
	r_CM_y += configs[k][3*(chain_start+i*chain_length + j) + 1];
	r_CM_z += configs[k][3*(chain_start+i*chain_length + j) + 2];
      }
      r_CM_x *= IdoubMPC; r_CM_y *= IdoubMPC; r_CM_z *= IdoubMPC;
      tmp = 0.0;
      for (j=0; j<chain_length; ++j) {
	dx = configs[k][3*(chain_start+i*chain_length + j)]     - r_CM_x;
	dy = configs[k][3*(chain_start+i*chain_length + j) + 1] - r_CM_y;
	dz = configs[k][3*(chain_start+i*chain_length + j) + 2] - r_CM_z;
	tmp += (SQR(dx) + SQR(dy) + SQR(dz));
      }
      tmp *= IdoubMPC;
      r_G  += sqrt(tmp);
      r_G2 += tmp;
      r_G4 += tmp*tmp;
    }
  }
  tmp = (double)(chain_n_chains*n_configs);
  rg[0] = r_G/tmp;
  rg[2] = r_G2/tmp;
  rg[1] = sqrt(rg[2] - rg[0]*rg[0]);
  rg[3] = sqrt(r_G4/tmp - rg[2]*rg[2]);
}

void calc_rh(double **_rh)
{
  int i, j, p;
  double dx,dy,dz, r_H=0.0,r_H2=0.0, *rh=NULL, ri=0.0, prefac, tmp;
  *_rh = rh = realloc(rh,2*sizeof(double));

  prefac = 0.5*(chain_length*(chain_length-1));
  for(p=0;p<chain_n_chains;p++) {
    ri=0.0;
    for(i=chain_start+chain_length*p;i<chain_start+chain_length*(p+1);i++) {
      for(j=i+1;j<chain_start+chain_length*(p+1);j++) {
	dx = partCfg[i].r.p[0]-partCfg[j].r.p[0]; dx*=dx;
	dy = partCfg[i].r.p[1]-partCfg[j].r.p[1]; dy*=dy;
	dz = partCfg[i].r.p[2]-partCfg[j].r.p[2]; dz*=dz;
	ri += 1.0/sqrt(dx+dy+dz);
      }
    }
    tmp = prefac/ri;
    r_H  += tmp;
    r_H2 += tmp*tmp;
  }
  tmp = (double)chain_n_chains;
  rh[0] = r_H/tmp;
  rh[1] = sqrt(r_H2/tmp - rh[0]*rh[0]);
}

void calc_rh_av(double **_rh)
{
  int i, j, p, k;
  double dx,dy,dz, r_H=0.0,r_H2=0.0, *rh=NULL, ri=0.0, prefac, tmp;
  *_rh = rh = realloc(rh,2*sizeof(double));

  prefac = 0.5*(chain_length*(chain_length-1));
  for(k=0; k<n_configs; k++) {
    for(p=0;p<chain_n_chains;p++) {
      ri=0.0;
      for(i=chain_start+chain_length*p;i<chain_start+chain_length*(p+1);i++) {
	for(j=i+1;j<chain_start+chain_length*(p+1);j++) {
	  dx = configs[k][3*i]  -configs[k][3*j];
	  dy = configs[k][3*i+1]-configs[k][3*j+1];
	  dz = configs[k][3*i+2]-configs[k][3*j+2];
	  ri += 1.0/sqrt(dx*dx + dy*dy + dz*dz);
	}
      }
      tmp = prefac/ri;
      r_H  += tmp;
      r_H2 += tmp*tmp;
    }
  }
  tmp = (double)chain_n_chains*n_configs;
  rh[0] = r_H/tmp;
  rh[1] = sqrt(r_H2/tmp - rh[0]*rh[0]);
}

void calc_internal_dist(double **_idf) {
  int i,j,k;
  double dx,dy,dz;
  double *idf=NULL;
  *_idf = idf = realloc(idf,chain_length*sizeof(double));

  idf[0] = 0.0;
  for (k=1; k < chain_length; k++) {
    idf[k] = 0.0;
    for (i=0; i<chain_n_chains; i++) {
      for (j=0; j < chain_length-k; j++) {
	dx = partCfg[chain_start+i*chain_length + j+k].r.p[0]
	  - partCfg[chain_start+i*chain_length + j].r.p[0];
	dy = partCfg[chain_start+i*chain_length + j+k].r.p[1]
	  - partCfg[chain_start+i*chain_length + j].r.p[1];
	dz = partCfg[chain_start+i*chain_length + j+k].r.p[2]
	  - partCfg[chain_start+i*chain_length +j].r.p[2];
	idf[k] += (SQR(dx) + SQR(dy) + SQR(dz));
      }
    }
    idf[k] = sqrt(idf[k] / (1.0*(chain_length-k)*chain_n_chains));
  }
}

void calc_internal_dist_av(double **_idf) {
  int i,j,k,n, i1,i2;
  double dx,dy,dz;
  double *idf=NULL;
  *_idf = idf = realloc(idf,chain_length*sizeof(double));

  idf[0] = 0.0;
  for (k=1; k < chain_length; k++) {
    idf[k] = 0.0;
    for (n=0; n<n_configs; n++) {
      for (i=0; i<chain_n_chains; i++) {
	for (j=0; j < chain_length-k; j++) {
	  i2 = chain_start+i*chain_length + j; i1 = i2 + k;
	  dx = configs[n][3*i1]   - configs[n][3*i2];
	  dy = configs[n][3*i1+1] - configs[n][3*i2+1];
	  dz = configs[n][3*i1+2] - configs[n][3*i2+2];
	  idf[k] += (SQR(dx) + SQR(dy) + SQR(dz));
	}
      }
    }
    idf[k] = sqrt(idf[k] / (1.0*(chain_length-k)*chain_n_chains*n_configs));
  }
}

void calc_bond_l(double **_bond_l) {
  int i,j;
  double dx,dy,dz, tmp;
  double *bond_l=NULL;
  *_bond_l = bond_l = realloc(bond_l,4*sizeof(double));

  bond_l[0] = bond_l[1] = bond_l[2] = 0.0; bond_l[3] = 20.0;
  for (i=0; i<chain_n_chains; i++) {
    for (j=0; j < chain_length-1; j++) {
      dx = partCfg[chain_start+i*chain_length + j+1].r.p[0]
	- partCfg[chain_start+i*chain_length + j].r.p[0];
      dy = partCfg[chain_start+i*chain_length + j+1].r.p[1]
	- partCfg[chain_start+i*chain_length + j].r.p[1];
      dz = partCfg[chain_start+i*chain_length + j+1].r.p[2]
	- partCfg[chain_start+i*chain_length +j].r.p[2];
      tmp = SQR(dx) + SQR(dy) + SQR(dz);
      bond_l[0] += sqrt(tmp);
      bond_l[1] += tmp;
      if (tmp > bond_l[2]) {bond_l[2] = tmp;}
      if (tmp < bond_l[3]) {bond_l[3] = tmp;}
    }
  }
  tmp = (double)(chain_length-1)*chain_n_chains;
  bond_l[0] = bond_l[0] / tmp;
  bond_l[1] = sqrt(bond_l[1]/tmp - SQR(bond_l[0]));
  bond_l[2] = sqrt(bond_l[2]);
  bond_l[3] = sqrt(bond_l[3]);
}

void calc_bond_l_av(double **_bond_l) {
  int i,j,n, i1,i2;
  double dx,dy,dz, tmp;
  double *bond_l=NULL;
  *_bond_l = bond_l = realloc(bond_l,4*sizeof(double));

  bond_l[0] = bond_l[1] = bond_l[2] = 0.0; bond_l[3] = 20.0;
  for (n=0; n<n_configs; n++) {
    for (i=0; i<chain_n_chains; i++) {
      for (j=0; j < chain_length-1; j++) {
	i2 = chain_start+i*chain_length + j; i1 = i2 + 1;
	dx = configs[n][3*i1]   - configs[n][3*i2];
	dy = configs[n][3*i1+1] - configs[n][3*i2+1];
	dz = configs[n][3*i1+2] - configs[n][3*i2+2];
	tmp = SQR(dx) + SQR(dy) + SQR(dz);
	bond_l[0] += sqrt(tmp);
	bond_l[1] += tmp;
	if (tmp > bond_l[2]) {bond_l[2] = tmp;}
      	if (tmp < bond_l[3]) {bond_l[3] = tmp;}
      }
    }
  }
  tmp = (double)(chain_length-1)*chain_n_chains*n_configs;
  bond_l[0] = bond_l[0] / tmp;
  bond_l[1] = sqrt(bond_l[1]/tmp - SQR(bond_l[0]));
  bond_l[2] = sqrt(bond_l[2]);
  bond_l[3] = sqrt(bond_l[3]);
}

void calc_bond_dist(double **_bdf, int ind_n) {
  int i,j=ind_n,k;
  double dx,dy,dz;
  double *bdf=NULL;
  *_bdf = bdf = realloc(bdf,chain_length*sizeof(double));

  bdf[0] = 0.0;
  for (k=1; k < chain_length; k++) {
    bdf[k] = 0.0;
    for (i=0; i<chain_n_chains; i++) {
      if (j < chain_length-k) {
	dx = partCfg[chain_start+i*chain_length + j+k].r.p[0]
	  - partCfg[chain_start+i*chain_length + j].r.p[0];
	dy = partCfg[chain_start+i*chain_length + j+k].r.p[1]
	  - partCfg[chain_start+i*chain_length + j].r.p[1];
	dz = partCfg[chain_start+i*chain_length + j+k].r.p[2]
	  - partCfg[chain_start+i*chain_length +j].r.p[2];
	bdf[k] += (SQR(dx) + SQR(dy) + SQR(dz));
      }
    }
    bdf[k] = sqrt(bdf[k] / (1.0*chain_n_chains));
  }
}

void calc_bond_dist_av(double **_bdf, int ind_n) {
  int i,j=ind_n,k,n, i1,i2;
  double dx,dy,dz;
  double *bdf=NULL;
  *_bdf = bdf = realloc(bdf,chain_length*sizeof(double));

  bdf[0] = 0.0;
  for (k=1; k < chain_length; k++) {
    bdf[k] = 0.0;
    for (n=0; n<n_configs; n++) {
      for (i=0; i<chain_n_chains; i++) {
	if (j < chain_length-k) {
	  i2 = chain_start+i*chain_length + j; i1 = i2 + k;
	  dx = configs[n][3*i1]   - configs[n][3*i2];
	  dy = configs[n][3*i1+1] - configs[n][3*i2+1];
	  dz = configs[n][3*i1+2] - configs[n][3*i2+2];
	  bdf[k] += (SQR(dx) + SQR(dy) + SQR(dz));
	}
      }
    }
    bdf[k] = sqrt(bdf[k] / (1.0*chain_n_chains*n_configs));
  }
}

void init_g123()
{
  int i, j, p;
  double cm_tmp[3];

  /* Save particles' current positions 
     (which'll be used as initial position later on) */
  partCoord_g = (float *) realloc(partCoord_g, 3*n_total_particles*sizeof(float));
  partCM_g = (float *) realloc(partCM_g, 3*chain_n_chains*sizeof(float));
  n_part_g = n_total_particles;
  n_chains_g = chain_n_chains;
  for(j=0; j<chain_n_chains; j++) {
    cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
    for(i=0; i<chain_length; i++) {
      p = chain_start+j*chain_length + i;
      partCoord_g[3*p]   = partCfg[p].r.p[0]; cm_tmp[0]+=partCfg[p].r.p[0];
      partCoord_g[3*p+1] = partCfg[p].r.p[1]; cm_tmp[1]+=partCfg[p].r.p[1];
      partCoord_g[3*p+2] = partCfg[p].r.p[2]; cm_tmp[2]+=partCfg[p].r.p[2];
    }
    partCM_g[3*j]   = cm_tmp[0]/(1.*chain_length);
    partCM_g[3*j+1] = cm_tmp[1]/(1.*chain_length);
    partCM_g[3*j+2] = cm_tmp[2]/(1.*chain_length);
  }
}

void calc_g123(double *_g1, double *_g2, double *_g3)
{
  /* - Mean square displacement of a monomer
     - Mean square displacement in the center of gravity of the chain itself
     - Motion of the center of mass */
  int i, j, p;
  double g1=0.0, g2=0.0, g3=0.0, cm_tmp[3];

  for(j=0; j<chain_n_chains; j++) {
    cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
    for(i=0; i<chain_length; i++) {
      p = chain_start+j*chain_length + i;
      cm_tmp[0]+=partCfg[p].r.p[0];
      cm_tmp[1]+=partCfg[p].r.p[1];
      cm_tmp[2]+=partCfg[p].r.p[2];
    }
    cm_tmp[0] /= (1.*chain_length);
    cm_tmp[1] /= (1.*chain_length);
    cm_tmp[2] /= (1.*chain_length);
    for(i=0; i<chain_length; i++) {
      p = chain_start+j*chain_length + i;
      g1 += SQR(partCfg[p].r.p[0]-partCoord_g[3*p])
	+ SQR(partCfg[p].r.p[1]-partCoord_g[3*p+1])
	+ SQR(partCfg[p].r.p[2]-partCoord_g[3*p+2]);
      g2 += SQR( (partCfg[p].r.p[0]-partCoord_g[3*p])
		 - (cm_tmp[0]-partCM_g[3*j]  ) ) 
	+ SQR( (partCfg[p].r.p[1]-partCoord_g[3*p+1])
	       - (cm_tmp[1]-partCM_g[3*j+1]) ) 
	+ SQR( (partCfg[p].r.p[2]-partCoord_g[3*p+2])
	       - (cm_tmp[2]-partCM_g[3*j+2]) );
    }
    g3 += SQR(cm_tmp[0]-partCM_g[3*j])
      + SQR(cm_tmp[1]-partCM_g[3*j+1])
      + SQR(cm_tmp[2]-partCM_g[3*j+2]);
  }
  *_g1 = g1 / (1.*chain_n_chains*chain_length);
  *_g2 = g2 / (1.*chain_n_chains*chain_length);
  *_g3 = g3 / (1.*chain_n_chains);
}

void calc_g1_av(double **_g1) {
  int i, j, p, t,k;
  double *g1=NULL;
  *_g1 = g1 = realloc(g1,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g1[k] = 0.0;
    for(t=0; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  g1[k] += SQR(configs[t+k][3*p]-configs[t][3*p])
	    + SQR(configs[t+k][3*p+1]-configs[t][3*p+1])
	    + SQR(configs[t+k][3*p+2]-configs[t][3*p+2]);
	}
      }
    }
    g1[k] /= ((double)chain_n_chains*chain_length*(n_configs-k));
  }
}

void calc_g2_av(double **_g2) {
  int i, j, p, t,k;
  double *g2=NULL, cm_tmp[3];
  *_g2 = g2 = realloc(g2,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g2[k] = 0.0;
    for(t=0; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  cm_tmp[0] += configs[t+k][3*p]   - configs[t][3*p];
	  cm_tmp[1] += configs[t+k][3*p+1] - configs[t][3*p+1];
	  cm_tmp[2] += configs[t+k][3*p+2] - configs[t][3*p+2];
	}
	cm_tmp[0] /= (1.*chain_length);	cm_tmp[1] /= (1.*chain_length);	cm_tmp[2] /= (1.*chain_length);
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  g2[k] += SQR( (configs[t+k][3*p]-configs[t][3*p]) - cm_tmp[0] )
	    + SQR( (configs[t+k][3*p+1]-configs[t][3*p+1])  - cm_tmp[1] ) 
	    + SQR( (configs[t+k][3*p+2]-configs[t][3*p+2])  - cm_tmp[2] );
	}
      }
    }
    g2[k] /= ((double)chain_n_chains*chain_length*(n_configs-k));
  }
}

void calc_g3_av(double **_g3) {
  int i, j, p, t,k;
  double *g3=NULL, cm_tmp[3];
  *_g3 = g3 = realloc(g3,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g3[k] = 0.0;
    for(t=0; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  cm_tmp[0] += configs[t+k][3*p]   - configs[t][3*p];
	  cm_tmp[1] += configs[t+k][3*p+1] - configs[t][3*p+1];
	  cm_tmp[2] += configs[t+k][3*p+2] - configs[t][3*p+2];
	}
	g3[k] += SQR(cm_tmp[0] / (1.*chain_length)) 
	  + SQR(cm_tmp[1] / (1.*chain_length)) 
	  + SQR(cm_tmp[2] / (1.*chain_length));
      }
    }
    g3[k] /= ((double)chain_n_chains*chain_length*(n_configs-k));
  }
}


void analyze_formfactor(double qmin, double qmax, int qbins, double **_ff) {
  int i,j,k,qi, cnt,cnt_max;
  double q,qfak, qr, dx,dy,dz, *r_ij=NULL, *ff=NULL;
  *_ff = ff = realloc(ff,(qbins+1)*sizeof(double));
  r_ij = realloc(r_ij,chain_length*(chain_length-1)/2*sizeof(double));

  qfak = pow((qmax/qmin),(1.0/qbins));
  for(qi=0; qi<=qbins; qi++) ff[qi] = 0.0;
  for(k=0; k < chain_n_chains; k++) {
    cnt = 0;
    /* Prepare distance matrice r_ij for current chain k */
    for(i=chain_start+k*chain_length; i<chain_start+(k+1)*chain_length; i++) {
      for(j=i+1; j<chain_start+(k+1)*chain_length;j++) {
	dx = partCfg[i].r.p[0]-partCfg[j].r.p[0]; dx*=dx;
	dy = partCfg[i].r.p[1]-partCfg[j].r.p[1]; dy*=dy;
	dz = partCfg[i].r.p[2]-partCfg[j].r.p[2]; dz*=dz;
	r_ij[cnt++] = sqrt(dx + dy + dz);
      }
    }
    q = qmin; cnt_max = cnt;
    /* Derive spherically averaged S(q) = 1/chain_length * Sum(i,j=1..chain_length)[sin(q*r_ij)/q*r_ij] for current chain k */
    for(qi=0; qi<=qbins; qi++) {
      ff[qi] += chain_length;
      for(cnt=0; cnt<cnt_max; cnt++) {
	qr = q*r_ij[cnt];
	ff[qi] += 2*sin(qr)/qr;
      }
      q *= qfak;
    }
  }
  for(qi=0; qi<=qbins; qi++) ff[qi] /= ((double)chain_length*chain_n_chains);
  free(r_ij);
}


void analyze_formfactor_av(double qmin, double qmax, int qbins, double **_ff) {
  int i,j,k,n,qi, cnt,cnt_max;
  double q,qfak, qr, dx,dy,dz, *r_ij=NULL, *ff=NULL;
  *_ff = ff = realloc(ff,(qbins+1)*sizeof(double));
  r_ij = realloc(r_ij,chain_length*(chain_length-1)/2*sizeof(double));

  qfak = pow((qmax/qmin),(1.0/qbins));
  for(qi=0; qi<=qbins; qi++) ff[qi] = 0.0;
  for(n=0; n < n_configs; n++) {
    for(k=0; k < chain_n_chains; k++) {
      cnt = 0;
      /* Prepare distance matrice r_ij for current chain k of configuration n */
      for(i=chain_start+k*chain_length; i<chain_start+(k+1)*chain_length; i++) {
	for(j=i+1; j<chain_start+(k+1)*chain_length;j++) {
	  dx = configs[n][3*i]  -configs[n][3*j];   dx*=dx;
	  dy = configs[n][3*i+1]-configs[n][3*j+1]; dy*=dy;
	  dz = configs[n][3*i+2]-configs[n][3*j+2]; dz*=dz;
	  r_ij[cnt++] = sqrt(dx + dy + dz);
	}
      }
      q = qmin; cnt_max = cnt;
      /* Derive spherically averaged S(q) = 1/chain_length * Sum(i,j=1..chain_length)[sin(q*r_ij)/q*r_ij] for current chain k */
      for(qi=0; qi<=qbins; qi++) {
	ff[qi] += chain_length;
	for(cnt=0; cnt<cnt_max; cnt++) {
	  qr = q*r_ij[cnt];
	  ff[qi] += 2*sin(qr)/qr;
	}
	q *= qfak;
      }
    }
  }
  for(qi=0; qi<=qbins; qi++) ff[qi] /= ((double)chain_length*chain_n_chains*n_configs);
  free(r_ij);
}

/****************************************************************************************
 *                                 chain structure commands parsing
 ****************************************************************************************/

int parse_chain_structure_info(Tcl_Interp *interp, int argc, char **argv)
{
  int m, i, pc;
  /* 'analyze set chains <chain_start> <n_chains> <chain_length>' */
  
  if (argc < 3) {
    Tcl_AppendResult(interp, "chain structure info consists of <start> <n> <length>", (char *)NULL);    
    return TCL_ERROR;
  }

  if (! (ARG0_IS_I(chain_start) && ARG1_IS_I(chain_n_chains) && ARG_IS_I(2, chain_length)))
    return TCL_ERROR;

  realloc_topology(chain_n_chains);
  pc = 0;
  for (m = 0; m < n_molecules; m++) {
    topology[m].type = 0;
    realloc_intlist(&topology[m].part, topology[m].part.n = chain_length);
    for (i = 0; i < chain_length; i++)
      topology[m].part.e[i] = pc++;
  }
 
  return TCL_OK;
}


void update_mol_ids_setchains() {
  Particle *p;
  int i, np, c;
  Cell *cell;
  
  for (c = 0; c < local_cells.n; c++) {
    cell = local_cells.cell[c];
    p  = cell->part;
    np = cell->n;
    for(i = 0; i < np; i++) {
      p[i].p.mol_id = floor((p[i].p.identity - chain_start)/(double)chain_length);
    }
  }
}


int check_and_parse_chain_structure_info(Tcl_Interp *interp, int argc, char **argv)
{
  if (argc > 0)
    if (parse_chain_structure_info(interp, argc, argv) != TCL_OK)
      return TCL_ERROR;
  
  if (!sortPartCfg()) {
    Tcl_AppendResult(interp, "for analyze, store particles consecutively starting with 0.",
		     (char *) NULL);
    return (TCL_ERROR);      
  }

  return TCL_OK;
}

int parse_re(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { re | <re> } [<chain_start> <n_chains> <chain_length>]' */
  char buffer[4*TCL_DOUBLE_SPACE+4];
  double *re;

  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR)
    return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) {
    Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
    return TCL_ERROR;
  }

  if (!average)
    calc_re(&re); 
  else {
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze re' to only look at current state!", (char *)NULL);
      return TCL_ERROR;
    }
    calc_re_av(&re);
  }

  sprintf(buffer,"%f %f %f %f",re[0],re[1],re[2],re[3]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  free(re);
  return (TCL_OK);
}


int parse_rg(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { rg | <rg> } [<chain_start> <n_chains> <chain_length>]' */
  char buffer[4*TCL_DOUBLE_SPACE+4];
  double *rg;
  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR)
    return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) {
    Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
    return TCL_ERROR;
  }
  if (!average)
    calc_rg(&rg); 
  else {
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze rg' to only look at current state!", (char *)NULL);
      return TCL_ERROR;
    }
    calc_rg_av(&rg);
  }

  sprintf(buffer,"%f %f %f %f",rg[0],rg[1],rg[2],rg[3]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  free(rg);
  return (TCL_OK);
}

int parse_rh(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { rh | <rh> } [<chain_start> <n_chains> <chain_length>]' */
  char buffer[2*TCL_DOUBLE_SPACE+2];
  double *rh;
  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR)
    return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) {
    Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL);
    return TCL_ERROR;
  }
  if (!average)
    calc_rh(&rh); 
  else {
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze rh' to only look at current state!", (char *)NULL);
      return TCL_ERROR;
    }
    else
      calc_rh_av(&rh);
  }

  sprintf(buffer,"%f %f",rh[0],rh[1]);
  Tcl_AppendResult(interp, buffer, (char *)NULL);

  free(rh);
  return (TCL_OK);
}

int parse_intdist(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { internal_dist | <internal_dist> } [<chain_start> <n_chains> <chain_length>]' */
  char buffer[TCL_DOUBLE_SPACE+2];
  int i;
  double *idf;

  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  if (!average)
    calc_internal_dist(&idf); 
  else {
    if (n_configs == 0) {
      Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
      Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze internal_dist' to only look at current state!", (char *)NULL);
      return TCL_ERROR;
    }
    else
      calc_internal_dist_av(&idf);
  }

  for (i=0; i<chain_length; i++) { 
    sprintf(buffer,"%f ",idf[i]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  }

  free(idf);
  return (TCL_OK);
}

int parse_bond_l(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { bond_l | <bond_l> } [<chain_start> <n_chains> <chain_length>]' */
  /*****************************************************************************/
  char buffer[4*TCL_DOUBLE_SPACE+2];
  double *bond_l;

  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  if (!average) calc_bond_l(&bond_l); 
  else if (n_configs == 0) {
    Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
    Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze bond_l' to only look at current state!", (char *)NULL);
    return TCL_ERROR; }
  else calc_bond_l_av(&bond_l);
  sprintf(buffer,"%f %f %f %f",bond_l[0],bond_l[1],bond_l[2],bond_l[3]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  free(bond_l);
  return (TCL_OK);
}

int parse_bond_dist(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { bond_dist | <bond_dist> } [index <index>] [<chain_start> <n_chains> <chain_length>]' */
  /***************************************************************************************************/
  char buffer[TCL_INTEGER_SPACE + TCL_DOUBLE_SPACE+2];
  double *bdf; int ind_n=0, i;

  if (argc >= 1 && !strncmp(argv[0], "index", strlen(argv[0]))) { ind_n = atoi(argv[1]); argc-=2; argv+=2; }
  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  if (ind_n < 0 || ind_n > chain_length-1) { 
    sprintf(buffer,"%d!",chain_length-1);
    Tcl_AppendResult(interp, "ERROR: <index> must be between 0 and ", buffer, (char *)NULL);  return TCL_ERROR; }
  if (ind_n >= chain_length/2) ind_n = (chain_length-1) - ind_n;
  if (!average) calc_bond_dist(&bdf,ind_n); 
  else if (n_configs == 0) {
    Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
    Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze internal_dist' to only look at current state!", (char *)NULL);
    return TCL_ERROR; }
  else calc_bond_dist_av(&bdf,ind_n);
  for (i=0; i<chain_length-ind_n; i++) { 
    sprintf(buffer,"%f ",bdf[i]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  }
  free(bdf);
  return (TCL_OK);
}

	   
int parse_g123(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze g123 [-init] [<chain_start> <n_chains> <chain_length>]' */
  /********************************************************************/
  char buffer[3*TCL_DOUBLE_SPACE+6];
  int init = 0;
  double g1, g2, g3;

  if (argc > 0 && ARG0_IS_S("-init")) {
    init = 1; argc--; argv++; 
  }
  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR)
    return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  
  if (init) { init_g123(); return TCL_OK; }
  if (partCoord_g == NULL || partCM_g == NULL) {
    Tcl_AppendResult(interp, "please call with -init first", (char *)NULL); return TCL_ERROR; }
  if (chain_n_chains != n_chains_g || n_total_particles != n_part_g) {
    fprintf(stderr, "%d %d %d %d\n", chain_n_chains, n_chains_g, n_total_particles, n_part_g);
    Tcl_AppendResult(interp, "initial config has different topology", (char *)NULL);
    return TCL_ERROR;      
  }
  calc_g123(&g1, &g2, &g3);
  sprintf(buffer,"{ %f %f %f }",g1, g2, g3);
  Tcl_AppendResult(interp, buffer, (char *)NULL);
  return (TCL_OK);
}

int parse_g_av(Tcl_Interp *interp, int what, int argc, char **argv)
{
  /* 'analyze { <g1> | <g2> | <g3> } [<chain_start> <n_chains> <chain_length>]' */
  /******************************************************************************/
  int i;
  char buffer[TCL_DOUBLE_SPACE+2];
  double *gx;

  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if ((argc != 0) && (argc != 3)) { Tcl_AppendResult(interp, "only chain structure info required", (char *)NULL); return TCL_ERROR; }
  if (n_configs == 0) { Tcl_AppendResult(interp, "no configurations found! Use 'analyze append' to save some!", (char *)NULL); return TCL_ERROR; }
  switch (what) {
  case 1:
    calc_g1_av(&gx); break;
  case 2:
    calc_g2_av(&gx); break;
  case 3:
    calc_g3_av(&gx); break;
  default: ;
  }
  for (i=0; i<n_configs; i++) { 
    sprintf(buffer,"%f ",gx[i]); Tcl_AppendResult(interp, buffer, (char *)NULL); 
  }
  free(gx);
  return (TCL_OK);
}


int parse_formfactor(Tcl_Interp *interp, int average, int argc, char **argv)
{
  /* 'analyze { formfactor | <formfactor> } <qmin> <qmax> <qbins> [<chain_start> <n_chains> <chain_length>]' */
  /***********************************************************************************************************/
  char buffer[2*TCL_DOUBLE_SPACE+4];
  int i;
  double qmin,qmax, q,qfak, *ff; int qbins;
  if (argc < 3) {
    Tcl_AppendResult(interp, "Wrong # of args! Usage: analyze formfactor <qmin> <qmax> <qbins> [<chain_start> <n_chains> <chain_length>]",
		     (char *)NULL);
    return (TCL_ERROR);
  } else {
    if (!ARG0_IS_D(qmin))
      return (TCL_ERROR);
    if (!ARG1_IS_D(qmax))
      return (TCL_ERROR);
    if (!ARG_IS_I(2, qbins))
      return (TCL_ERROR);
    argc-=3; argv+=3;
  }
  if (check_and_parse_chain_structure_info(interp, argc, argv) == TCL_ERROR) return TCL_ERROR;
  if (qbins <=0) {
    Tcl_AppendResult(interp, "Nothing to be done - choose <qbins> greater zero to get S(q)!",(char *)NULL); return TCL_ERROR;
  }

  if (qmin <= 0.) {
    Tcl_AppendResult(interp, "formfactor S(q) requires qmin > 0", (char *)NULL);
    return TCL_ERROR;
  }
  if (qmax <= qmin) {
    Tcl_AppendResult(interp, "formfactor S(q) requires qmin < qmax", (char *)NULL);
    return TCL_ERROR;
  }

  if (!average) analyze_formfactor(qmin, qmax, qbins, &ff);
  else if (n_configs == 0) {
    Tcl_AppendResult(interp, "no configurations found! ", (char *)NULL);
    Tcl_AppendResult(interp, "Use 'analyze append' to save some, or 'analyze formfactor ...' to only look at current state!",
		     (char *)NULL);
    return TCL_ERROR; }
  else analyze_formfactor_av(qmin, qmax, qbins, &ff);
  
  q = qmin; qfak = pow((qmax/qmin),(1.0/qbins));
  for(i=0; i<=qbins; i++) { sprintf(buffer,"{%f %f} ",q,ff[i]); q*=qfak; Tcl_AppendResult(interp, buffer, (char *)NULL); }
  free(ff);
  return (TCL_OK);
}


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
/** \file statistics_chain.cpp
    Implementation of \ref statistics_chain.hpp "statistics_chain.hpp".
*/
#include "statistics.hpp"
#include "utils.hpp"
#include "topology.hpp"
#include "communication.hpp"
#include "cells.hpp"
#include "grid.hpp"

/** Particles' initial positions (needed for g1(t), g2(t), g3(t) in \ref tclcommand_analyze) */
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



void calc_re(double **_re)
{
  int i;
  double dx,dy,dz;
  double dist=0.0,dist2=0.0,dist4=0.0;
  double *re=NULL, tmp;
  *_re = re = (double*)realloc(re,4*sizeof(double));

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
  *_re = re = (double*)realloc(re,4*sizeof(double));

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
  int i, j, p;
  double dx,dy,dz, r_CM_x,r_CM_y,r_CM_z;
  double r_G=0.0,r_G2=0.0,r_G4=0.0;
  double *rg=NULL,IdoubMPC, tmp;
  double M;
  *_rg = rg = (double*)realloc(rg,4*sizeof(double));


  for (i=0; i<chain_n_chains; i++) {
    M=0.0;
    r_CM_x = r_CM_y = r_CM_z = 0.0;
    IdoubMPC = 1./(double)chain_length;
    for (j=0; j<chain_length; j++) {
        p = chain_start+i*chain_length + j;
      r_CM_x += partCfg[p].r.p[0]*PMASS(partCfg[p]);
      r_CM_y += partCfg[p].r.p[1]*PMASS(partCfg[p]);
      r_CM_z += partCfg[p].r.p[2]*PMASS(partCfg[p]);
      M += PMASS(partCfg[p]);
    }
    r_CM_x /= M;
    r_CM_y /= M;
    r_CM_z /= M;
    tmp = 0.0;
    for (j=0; j<chain_length; ++j) {
      p = chain_start+i*chain_length + j;
      dx = partCfg[p].r.p[0] - r_CM_x;
      dy = partCfg[p].r.p[1] - r_CM_y;
      dz = partCfg[p].r.p[2] - r_CM_z;
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
  int i, j, k, p;
  double dx,dy,dz, r_CM_x,r_CM_y,r_CM_z;
  double r_G=0.0,r_G2=0.0, r_G4=0.0;
  double *rg=NULL, IdoubMPC, tmp;
  double M;
  *_rg = rg = (double*)realloc(rg,4*sizeof(double));

  IdoubMPC = 1./(double)chain_length;
  for (k=0; k<n_configs; k++) {
    for (i=0; i<chain_n_chains; i++) {
      M=0.0;
      r_CM_x = r_CM_y = r_CM_z = 0.0;
      for (j=0; j<chain_length; j++) {
        p = chain_start+i*chain_length + j;
	r_CM_x += configs[k][3*p]*PMASS(partCfg[p]);
	r_CM_y += configs[k][3*p + 1]*PMASS(partCfg[p]);
	r_CM_z += configs[k][3*p + 2]*PMASS(partCfg[p]);
        M += PMASS(partCfg[p]);
      }
      r_CM_x /= M; r_CM_y /= M; r_CM_z /= M;
      tmp = 0.0;
      for (j=0; j<chain_length; ++j) {
        p = chain_start+i*chain_length + j;
	dx = configs[k][3*p]     - r_CM_x;
	dy = configs[k][3*p + 1] - r_CM_y;
	dz = configs[k][3*p + 2] - r_CM_z;
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
  *_rh = rh = (double*)realloc(rh,2*sizeof(double));

  prefac = 0.5*chain_length*chain_length; /* 1/N^2 is not a normalization factor */
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
  *_rh = rh = (double*)realloc(rh,2*sizeof(double));

  prefac = 0.5*chain_length*chain_length;
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
  *_idf = idf = (double*)realloc(idf,chain_length*sizeof(double));

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
  *_idf = idf = (double*)realloc(idf,chain_length*sizeof(double));

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
  *_bond_l = bond_l = (double*)realloc(bond_l,4*sizeof(double));

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
  *_bond_l = bond_l = (double*)realloc(bond_l,4*sizeof(double));

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
  *_bdf = bdf = (double*)realloc(bdf,chain_length*sizeof(double));

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
  *_bdf = bdf = (double*)realloc(bdf,chain_length*sizeof(double));

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
  double cm_tmp[3], M;

  /* Save particles' current positions 
     (which'll be used as initial position later on) */
  partCoord_g = (float *) realloc(partCoord_g, 3*n_part*sizeof(float));
  partCM_g = (float *) realloc(partCM_g, 3*chain_n_chains*sizeof(float));
  n_part_g = n_part;
  n_chains_g = chain_n_chains;
  for(j=0; j<chain_n_chains; j++) {
    cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
    M=0.0;
    for(i=0; i<chain_length; i++) {
      p = chain_start+j*chain_length + i;
      partCoord_g[3*p]   = partCfg[p].r.p[0]; cm_tmp[0]+=partCfg[p].r.p[0]*PMASS(partCfg[p]);
      partCoord_g[3*p+1] = partCfg[p].r.p[1]; cm_tmp[1]+=partCfg[p].r.p[1]*PMASS(partCfg[p]);
      partCoord_g[3*p+2] = partCfg[p].r.p[2]; cm_tmp[2]+=partCfg[p].r.p[2]*PMASS(partCfg[p]);
      M += PMASS(partCfg[p]);
    }
    partCM_g[3*j]   = cm_tmp[0]/M;
    partCM_g[3*j+1] = cm_tmp[1]/M;
    partCM_g[3*j+2] = cm_tmp[2]/M;
  }
}

void calc_g123(double *_g1, double *_g2, double *_g3)
{
  /* - Mean square displacement of a monomer
     - Mean square displacement in the center of gravity of the chain itself
     - Motion of the center of mass */
  int i, j, p;
  double g1=0.0, g2=0.0, g3=0.0, cm_tmp[3];
  double M;

  for(j=0; j<chain_n_chains; j++) {
    cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
    M=0.0;
    for(i=0; i<chain_length; i++) {
      p = chain_start+j*chain_length + i;
      cm_tmp[0]+=partCfg[p].r.p[0]*PMASS(partCfg[p]);
      cm_tmp[1]+=partCfg[p].r.p[1]*PMASS(partCfg[p]);
      cm_tmp[2]+=partCfg[p].r.p[2]*PMASS(partCfg[p]);
      M += PMASS(partCfg[p]);
    }
    cm_tmp[0] /= M;
    cm_tmp[1] /= M;
    cm_tmp[2] /= M;
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

void calc_g1_av(double **_g1, int window, double weights[3]) {
  int i, j, p, t,k, cnt;
  double *g1=NULL;
  *_g1 = g1 = (double*)realloc(g1,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g1[k] = 0.0;
    cnt = window?n_configs-k:1;
    for(t=window?0:n_configs-k-1; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  g1[k] += weights[0]*SQR(configs[t+k][3*p]-configs[t][3*p])
	    + weights[1]*SQR(configs[t+k][3*p+1]-configs[t][3*p+1])
	    + weights[2]*SQR(configs[t+k][3*p+2]-configs[t][3*p+2]);
	}
      }
    }
    g1[k] /= ((double)chain_n_chains*chain_length*cnt);
  }
}

void calc_g2_av(double **_g2, int window, double weights[3]) {
  int i, j, p, t,k,cnt;
  double *g2=NULL, cm_tmp[3];
  double M;
  *_g2 = g2 = (double*)realloc(g2,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g2[k] = 0.0;
    cnt = window?n_configs-k:1;
    for(t=window?0:n_configs-k-1; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
        M=0.0;
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  cm_tmp[0] += (configs[t+k][3*p]   - configs[t][3*p])*PMASS(partCfg[p]);
	  cm_tmp[1] += (configs[t+k][3*p+1] - configs[t][3*p+1])*PMASS(partCfg[p]);
	  cm_tmp[2] += (configs[t+k][3*p+2] - configs[t][3*p+2])*PMASS(partCfg[p]);
          M += PMASS(partCfg[p]);
	}
	cm_tmp[0] /= M;	cm_tmp[1] /= M;	cm_tmp[2] /= M;
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  g2[k] += weights[0]*SQR( (configs[t+k][3*p]-configs[t][3*p]) - cm_tmp[0] )
	    + weights[1]*SQR( (configs[t+k][3*p+1]-configs[t][3*p+1])  - cm_tmp[1] ) 
	    + weights[2]*SQR( (configs[t+k][3*p+2]-configs[t][3*p+2])  - cm_tmp[2] );
	}
      }
    }
    g2[k] /= ((double)chain_n_chains*chain_length*cnt);
  }
}

void calc_g3_av(double **_g3, int window, double weights[3]) {
  int i, j, p, t,k,cnt;
  double *g3=NULL, cm_tmp[3];
  double M;
  *_g3 = g3 = (double*)realloc(g3,n_configs*sizeof(double));

  for(k=0; k < n_configs; k++) {
    g3[k] = 0.0;
    cnt = window?n_configs-k:1;
    for(t=window?0:n_configs-k-1; t < n_configs-k; t++) {
      for(j=0; j<chain_n_chains; j++) {
	cm_tmp[0] = cm_tmp[1] = cm_tmp[2] = 0.0;
        M=0.0;
	for(i=0; i<chain_length; i++) {
	  p = chain_start+j*chain_length + i;
	  cm_tmp[0] += (configs[t+k][3*p]   - configs[t][3*p])*PMASS(partCfg[p]);
	  cm_tmp[1] += (configs[t+k][3*p+1] - configs[t][3*p+1])*PMASS(partCfg[p]);
	  cm_tmp[2] += (configs[t+k][3*p+2] - configs[t][3*p+2])*PMASS(partCfg[p]);
          M += PMASS(partCfg[p]);
	}
	g3[k] += (weights[0]*SQR(cm_tmp[0]) + weights[1]*SQR(cm_tmp[1]) + weights[2]*SQR(cm_tmp[2]))/SQR(M);
      }
    }
    g3[k] /= ((double)chain_n_chains*cnt);
  }
}


void analyze_formfactor(double qmin, double qmax, int qbins, double **_ff) {
  int i,j,k,qi, cnt,cnt_max;
  double q,qfak, qr, dx,dy,dz, *r_ij=NULL, *ff=NULL;
  *_ff = ff = (double*)realloc(ff,(qbins+1)*sizeof(double));
  r_ij = (double*)realloc(r_ij,chain_length*(chain_length-1)/2*sizeof(double));

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
  *_ff = ff = (double*)realloc(ff,(qbins+1)*sizeof(double));
  r_ij = (double*)realloc(r_ij,chain_length*(chain_length-1)/2*sizeof(double));

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


void analyze_rdfchain(double r_min, double r_max, int r_bins, double **_f1, double **_f2, double **_f3) {
  int i, j, ind, c_i, c_j, mon;
  double bin_width, inv_bin_width, factor, r_in, r_out, bin_volume, dist, chain_mass,
    *cm=NULL, *min_d=NULL, *f1=NULL, *f2=NULL, *f3=NULL;
  
  *_f1 = f1 = (double*)realloc(f1,r_bins*sizeof(double));
  *_f2 = f2 = (double*)realloc(f2,r_bins*sizeof(double));
  *_f3 = f3 = (double*)realloc(f3,r_bins*sizeof(double));
  cm = (double*)realloc(cm,(chain_n_chains*3)*sizeof(double));
  min_d = (double*)realloc(min_d,(chain_n_chains*chain_n_chains)*sizeof(double));
  for(i=0; i<r_bins; i++) {
    f1[i] = f2[i] = f3[i] = 0.0;
  }
  bin_width     = (r_max-r_min) / (double)r_bins;
  inv_bin_width = 1.0 / bin_width;
  for(c_i=0; c_i<chain_n_chains; c_i++) {
    cm[3*c_i] = cm[3*c_i+1] = cm[3*c_i+2] = 0.0;
    for(c_j=(c_i+1); c_j<chain_n_chains; c_j++) {
      min_d[c_i*chain_n_chains+c_j] = box_l[0] + box_l[1] + box_l[2];
    }
  }    
  for(c_i=0; c_i<chain_n_chains; c_i++) {
    for(i=0; i<chain_length; i++) {
      mon = chain_start + c_i*chain_length + i;
      cm[3*c_i]+= partCfg[mon].r.p[0]*PMASS(partCfg[mon]);
      cm[3*c_i+1]+= partCfg[mon].r.p[1]*PMASS(partCfg[mon]);
      cm[3*c_i+2]+= partCfg[mon].r.p[2]*PMASS(partCfg[mon]);	  
      for(c_j=(c_i+1); c_j<chain_n_chains; c_j++) {
        for(j=0; j<chain_length; j++) {
          dist = min_distance(partCfg[mon].r.p, partCfg[chain_start + c_j*chain_length+j].r.p);		
          if ((dist > r_min) && (dist < r_max)) {
            ind = (int) ( (dist - r_min)*inv_bin_width ); 
            f1[ind]+= 1.0; 
          }        
          if (dist<min_d[c_i*chain_n_chains+c_j]) min_d[c_i*chain_n_chains+c_j] = dist;
        }
      }
    }       
  }
  chain_mass = 0;
  for(i=0; i<chain_length; i++) chain_mass+= PMASS(partCfg[i]);
  for(c_i=0; c_i<chain_n_chains; c_i++) {  
    cm[3*c_i]/= chain_mass;
    cm[3*c_i+1]/= chain_mass;
    cm[3*c_i+2]/= chain_mass;
    for(c_j=(c_i+1); c_j<chain_n_chains; c_j++) {
      dist = min_d[c_i*chain_n_chains+c_j];
      if ((dist > r_min) && (dist < r_max)) {
        ind = (int) ( (dist - r_min)*inv_bin_width );
        f3[ind]+= 1.0; 
      }  
    }
  }    
  for(c_i=0; c_i<chain_n_chains; c_i++) {  
    for(c_j=(c_i+1); c_j<chain_n_chains; c_j++) {
      dist = min_distance(&cm[3*c_i],&cm[3*c_j]);
      if ((dist > r_min) && (dist < r_max)) {
        ind = (int) ( (dist - r_min)*inv_bin_width );
        f2[ind]+= 1.0; 
      }  
    }
  }
  /* normalization */
  factor = box_l[0]*box_l[1]*box_l[2] *1.5/PI/(chain_n_chains*(chain_n_chains-1));
  for(i=0; i<r_bins; i++) {
    r_in       = i*bin_width + r_min; 
    r_out      = r_in + bin_width;
    bin_volume = r_out*r_out*r_out - r_in*r_in*r_in;
    f1[i] *= factor / (bin_volume * chain_length*chain_length);
    f2[i] *= factor / bin_volume;
    f3[i] *= factor / bin_volume;
  }  
}

#ifdef ELECTROSTATICS
void analyze_cwvac(int maxtau, int interval, double **_avac, double **_evac) {
    int i,j,k,ic,n,mon,tau;
    double *vac=NULL, *avac=NULL, *evac=NULL;
    double *varr=NULL;
    double vcm[3];
    
    // result
    vac = (double*)realloc(vac,(maxtau+1)*sizeof(double));
    *_avac = avac = (double*)realloc(avac,(maxtau+1)*sizeof(double));
    *_evac = evac = (double*)realloc(evac,(maxtau+1)*sizeof(double));
    
    // set up velarray
    varr = (double*)realloc(varr,(3*(n_configs-1)*n_part)*sizeof(double));
    
    // fill varr
    for(k=0;k<n_configs-1;k++) {
      for(i=0;i<n_part;i++) {
    	for(j=0;j<3;j++) {
          varr[k*3*n_part+i*3+j]=configs[k+1][3*i+j]-configs[k][3*i+j];	
    	}
      }
    } 
    
    // loop over chain
    for(tau=0;tau<=maxtau;tau++) { avac[tau] = 0.0; evac[tau] = 0.0; }
    for(n=0;n<chain_n_chains;n++) {
      for(tau=0;tau<=maxtau;tau++) { vac[tau] = 0.0; }
      // loop over configuration (interval)
      for(ic=(n_configs-2);ic>=0;ic -= interval) {
      	if(ic < maxtau) break;
    	// calc center of mass velocity
    	vcm[0] = vcm[1] = vcm[2] = 0.0;
    	for(k=0;k<chain_length;k++) {
    	  mon = chain_start+n*chain_length+k;
    	  for(j=0;j<3;j++) {
    	    vcm[j] +=  varr[ic*3*n_part+mon*3+j];
    	  }
    	}
    	for(j=0;j<3;j++) {
    	  vcm[j] /=  chain_length;
    	}
        // loop over tau (maxtau)
        for(tau=0;tau<=maxtau;tau++) {
          // loop over particles 
          for(k=0;k<n_part;k++) {
            vac[tau] += partCfg[k].p.q*(
                            varr[(ic-tau)*3*n_part+k*3  ]*vcm[0]
                         +  varr[(ic-tau)*3*n_part+k*3+1]*vcm[1]
                         +  varr[(ic-tau)*3*n_part+k*3+2]*vcm[2]);
            }
        }
      }
    
      // normalize vac and sum up to average
      for(tau=0;tau<=maxtau;tau++) { 
      	vac[tau] /= ((n_configs-maxtau)/interval); 
      	avac[tau] += vac[tau];
      	evac[tau] += vac[tau]*vac[tau];
      }
    }
   
    // normalize
    // NOTE: normalization with time step and KbT has to be done on script level 
    for(tau=0;tau<=maxtau;tau++) { 
    	evac[tau] = sqrt((evac[tau]-avac[tau]*avac[tau]/chain_n_chains)/(chain_n_chains*(chain_n_chains-1)));
    	avac[tau] /= chain_n_chains;
    }
    
    free(vac); free(varr);
 }
#endif


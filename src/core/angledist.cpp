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
/** \file angledist.cpp
 *
 *  Implementation of \ref angledist.hpp
 */
#include "angledist.hpp"

#ifdef BOND_ANGLEDIST
#include "communication.hpp"
#include "constraint.hpp"
#include "grid.hpp"

int angledist_set_params(int bond_type, double bend, double phimin, double distmin, double phimax, double distmax)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.angledist.bend = bend;
  bonded_ia_params[bond_type].p.angledist.phimin = phimin;
  bonded_ia_params[bond_type].p.angledist.distmin = distmin;
  bonded_ia_params[bond_type].p.angledist.phimax = phimax;
  bonded_ia_params[bond_type].p.angledist.distmax = distmax;
#ifdef BOND_ANGLEDIST_COSINE
#error angledist not implemented for BOND_ANGLEDIST_COSINE
#endif
#ifdef BOND_ANGLEDIST_COSSQUARE
#error angledist not implemented for BOND_ANGLEDIST_COSSQUARE
#endif
  bonded_ia_params[bond_type].type = BONDED_IA_ANGLEDIST;
  bonded_ia_params[bond_type].num = 2;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

/** Function to calculate wall distance and phi0(dist).
    Called by \ref calc_angledist_force */
static double calc_angledist_param(Particle *p_mid, Particle *p_left,
				   Particle *p_right, Bonded_ia_parameters *iaparams)
{
  double cosine=0.0, vec1[3], vec2[3], d1i=0.0, d2i=0.0, dist1=0.0, dist2=0.0, phi0=0.0;
  //  double pwdist=0.0, pwdist0=0.0, pwdist1=0.0;
  double normal, folded_pos[3], phimn=0.0, distmn=0.0, phimx=0.0, distmx=0.0, drange=0.0;
  double pwdist[n_constraints],pwdistmin=0.0;
  Constraint_wall wall;
  int j, k;
  int img[3];

  wall.d = 0;

  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist1 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist1);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* vectors are normalised so cosine is just cos(angle_between_vec1_and_vec2) */
  cosine = scalar(vec1, vec2);
  if ( cosine >  TINY_COS_VALUE)  cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  //  fac    = iaparams->p.angledist.bend; /* spring constant from .tcl file */
  phimn  = iaparams->p.angledist.phimin;
  distmn = iaparams->p.angledist.distmin;
  phimx  = iaparams->p.angledist.phimax;
  distmx = iaparams->p.angledist.distmax;

  /* folds coordinates of p_mid into original box */
  memcpy(folded_pos, p_mid->r.p, 3*sizeof(double));
  memcpy(img, p_mid->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);

  /* Calculates distance between p_mid and constraint */
  for(k=0;k<n_constraints;k++) {
    pwdist[k]=0.0;
  }
  for(k=0;k<n_constraints;k++) {
    switch(constraints[k].type) {
      case CONSTRAINT_WAL: 

      /* dist is distance of wall from origin */
      wall=constraints[k].c.wal;

      /* check that constraint vector is normalised */
      normal=0.0;
      for(j=0;j<3;j++) normal += wall.n[j] * wall.n[j];
      if (sqrt(normal) != 1.0) {
        for(j=0;j<3;j++) wall.n[j]=wall.n[j]/normal;
      }

      /* pwdist is distance of wall from p_mid */
      pwdist[k]=-1.0 * constraints[k].c.wal.d;
      for(j=0;j<3;j++) {
        pwdist[k] += folded_pos[j] * constraints[k].c.wal.n[j];
      }
      if (k==0) {
        pwdistmin=pwdist[k];
      }
      if (pwdist[k] <= pwdistmin) {
        pwdistmin = pwdist[k];
      }
      break;
    }
  }

  /*get phi0(z)*/
  if (pwdistmin <= distmn) {
    phi0 = phimn;
    //    fprintf(stdout,"\nIn angledist_set_params:  z_p_mid=%f  pwdistmin=%f  distmn=%f  ",folded_pos[2],pwdistmin,distmn);
    //    fprintf(stdout,"  phi0=%f\n",phi0*180.0/PI);
  }
  else if (pwdistmin >= distmx && pwdistmin <= box_l[2]-wall.d-distmx) {
    phi0 = phimx;
  }
  else {
    drange = (pwdistmin-distmn)*PI/(distmx-distmn);
    phi0 = ((cos(drange-PI)+1.0)*(phimx-phimn))*0.5+phimn;
  //  fprintf(stdout,"\nIn angledist_set_params:  z_p_mid=%f  pwdistmin=%f  box_lz/2=%f  ",folded_pos[2],pwdistmin,box_l[2]/2.0);
  //  fprintf(stdout,"  phi0=%f\n",phi0*180.0/PI);
  }

  return phi0;
}

int calc_angledist_force(Particle *p_mid, Particle *p_left,
				  Particle *p_right,
                                  Bonded_ia_parameters *iaparams,
				  double force1[3], double force2[3])
{
  double cosine=0.0, vec1[3], vec2[3], d1i=0.0, d2i=0.0, dist2, fac=0.0, f1=0.0, f2=0.0;
  int j;

  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar product of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  fac    = iaparams->p.angledist.bend;

#ifdef BOND_ANGLEDIST_HARMONIC
  {
    double phi,sinphi;
    /* NOTE: The angledist is ONLY implemented for the HARMONIC case */
    double phi0 = calc_angledist_param(p_mid, p_left, p_right, iaparams);

      if ( cosine >  TINY_COS_VALUE) cosine =  TINY_COS_VALUE;
    if ( cosine < -TINY_COS_VALUE) cosine = -TINY_COS_VALUE;
    phi =  acos(-cosine);
    sinphi = sin(phi);
    if ( sinphi < TINY_SIN_VALUE ) sinphi = TINY_SIN_VALUE;
    fac *= (phi - phi0)/sinphi;
    //    fprintf(stdout,"\n force:  z_pmid=%f, phi0=%f  phi=%f fac=%f",p_mid->r.p[2],phi0*180.0/PI,phi*180.0/PI,fac);
  }
#endif

#ifdef BOND_ANGLEDIST_COSINE
  #error angledist ONLY implemented for harmonic case!
#endif
#ifdef BOND_ANGLEDIST_COSSQUARE
  #error angledist ONLY implemented for harmonic case!
#endif

  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    force1[j] = (f1-f2);
    force2[j] = -f1;
  }
  return 0;
}

int angledist_energy(Particle *p_mid, Particle *p_left, Particle *p_right, 
		     Bonded_ia_parameters *iaparams, double *_energy)
{
  int j;
  double cosine=0.0, d1i, d2i, dist1, dist2;
  double vec1[3], vec2[3];

  //  if (phi0 < PI) {
  //    fprintf(stdout,"\nIn angledist_energy:  z_p_mid=%f, phi0=%f\n",p_mid->r.p[2],phi0*180.0/PI);
  //  }

  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist1 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist1);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_right to p_mid */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  if ( cosine >  TINY_COS_VALUE)  cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;

#ifdef BOND_ANGLEDIST_HARMONIC
  {
    double phi;
    double phi0=calc_angledist_param(p_mid, p_left, p_right, iaparams);
    phi =  acos(-cosine);
    *_energy = 0.5*iaparams->p.angledist.bend*SQR(phi - phi0);
    //    fprintf(stdout,"\n energy:  z_pmid=%f  bend=%f  phi0=%f  phi=%f energy=%f",p_mid->r.p[2],iaparams->p.angledist.bend,phi0*180.0/PI,phi*180.0/PI,0.5*iaparams->p.angledist.bend*SQR(phi - phi0));

  }
#endif
#ifdef BOND_ANGLEDIST_COSINE
  #error angledist ONLY implemented for harmonic case!
#endif
#ifdef BOND_ANGLEDIST_COSSQUARE
  #error angledist ONLY implemented for harmonic case!
#endif
  return 0;
}

#endif


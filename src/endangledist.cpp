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
/** \file endangledist.cpp
 *
 *  Implementation of \ref endangledist.hpp
 */
#include "utils.hpp"
#include "interaction_data.hpp"
#include "communication.hpp"
#include "constraint.hpp"
#include "grid.hpp"

#ifdef BOND_ENDANGLEDIST

int endangledist_set_params(int bond_type, double bend, double phi0 ,double distmin, double distmax)
{
  if(bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.endangledist.bend = bend;
  bonded_ia_params[bond_type].p.endangledist.phi0 = phi0;
  bonded_ia_params[bond_type].p.endangledist.distmin = distmin;
  bonded_ia_params[bond_type].p.endangledist.distmax = distmax;

  bonded_ia_params[bond_type].type = BONDED_IA_ENDANGLEDIST;
  /* Normally LENGTH=1 ANGLE=2 DIHEDRAL=3 
   * Here angle only requires one particle (other reference is wall constraint)
   */
  bonded_ia_params[bond_type].num  = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return ES_OK;
}

#ifdef CONSTRAINTS
/// Calculates the minimum distance between a particle to any wall constraint
static double calc_pwdist(Particle *p1, Bonded_ia_parameters *iaparams, int *clconstr)
{
  int j,k,img[3];
  double distwallmin=0.0, normal=0.0;
  double folded_pos_p1[3];
  double pwdist[n_constraints];
  Constraint_wall wall;

  /*fprintf(stdout,"  Entering calc_pwdist:\n");*/

  /* folds coordinates of p_left into original box */
  memcpy(folded_pos_p1, p1->r.p, 3*sizeof(double));
  memcpy(img, p1->l.i, 3*sizeof(int));
  fold_position(folded_pos_p1, img);
  /*fprintf(stdout,"        p1= %9.6f %9.6f %9.6f\n",p1->r.p[0],p1->r.p[1],p1->r.p[2]);*/

  /* Gets and tests wall data */
  for(k=0;k<n_constraints;k++) {
    switch(constraints[k].type) {
      case CONSTRAINT_WAL: 
      wall=constraints[k].c.wal;
      /* check that constraint wall normal is normalised */
      for(j=0;j<3;j++) normal += wall.n[j] * wall.n[j];
      if (sqrt(normal) != 1.0) {
        for(j=0;j<3;j++) wall.n[j]=wall.n[j]/normal;
      }
      break;
    }
  }

  /* Calculate distance of end particle from closest wall */
  for(k=0;k<n_constraints;k++) {
    switch(constraints[k].type) {
      case CONSTRAINT_WAL:
      wall=constraints[k].c.wal;
      /* distwallmin is distance of closest wall from p1 */
      pwdist[k]=-1.0 * wall.d;
      for(j=0;j<3;j++) {
        pwdist[k] += folded_pos_p1[j] * wall.n[j];
      }
      if (k==0) {
        distwallmin=pwdist[k];
      } else {
        if (pwdist[k] < distwallmin) {
          distwallmin = pwdist[k];
          *clconstr =  k;
        }
      }
      /*fprintf(stdout,"  k=%d  clconstr=%d\n",k,*clconstr);*/
      break;
    }
  }
  /*
    if (distwallmin <= iaparams->p.endangledist.distmax) {
    fprintf(stdout,"  clconstr=%d  distwallmin=%f  distmx=%f\n",*clconstr,distwallmin,distmx);
  }
  */
  return distwallmin;
}

/// Calculate angle that p1--p2 makes with wall constraint
static double calc_pwangle(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, int *constr)
{
  int j;
  double dist,di,cosine,phi;
  double vec[3];

  /* vector from p1 to p2 */
  get_mi_vector(vec, p2->r.p, p1->r.p);
  dist = sqrlen(vec);
  di = 1.0/sqrt(dist);
  for(j=0;j<3;j++) vec[j] *= di;
  /*
  fprintf(stdout,"Normalised: p1= %9.6f %9.6f %9.6f   p1= %9.6f %9.6f %9.6f   vec= %9.6f %9.6f %9.6f\n",p1->r.p[0],p1->r.p[1],p1->r.p[2],p2->r.p[0],p2->r.p[1],p2->r.p[2],vec[0],vec[1],vec[2]);
  */
  /* vectors are normalised so cosine is just cos(angle_between_vec1_and_vec2)
   * Wall is closest wall to particle
   */

  cosine = scalar(vec, constraints[*constr].c.wal.n);
  if ( cosine >  TINY_COS_VALUE)  cosine =  TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  phi=acos(cosine);

  /*
  fprintf(stdout,"Angle with wall 0=%f  ",(acos(scalar(vec, constraints[0].c.wal.n)))*180.0/PI);
  fprintf(stdout,"Angle with wall 1=%f  ",(acos(scalar(vec, constraints[1].c.wal.n)))*180.0/PI);
  fprintf(stdout,"dxy=%f  dz=%f  angle=%f\n",sqrt(vec[0]*vec[0]+vec[1]*vec[1]),vec[2],atan(sqrt(vec[0]*vec[0]+vec[1]*vec[1])/vec[2])*180.0/PI);
  fprintf(stdout,"Angle with closest wall %d=%f  ",*constr,(acos(scalar(vec, constraints[*constr].c.wal.n)))*180.0/PI);
  */
  return phi;
}
#endif

///
int calc_endangledist_pair_force(Particle *p1, Particle *p2, 
                                 Bonded_ia_parameters *iaparams, 
                                 double dx[3], double force1[3], double force2[3])
{

  int i=0;
  int clconstr=0;
  //  double distwallmin=0.0, distmx, distmn;
  double bend=0.0,phieq=0.0,phi=0.0,distwallmin=0.0, distmx, distmn, dist, di;
  double smooth, sinphi, cosphi, fac_a, fac_b, gradharm1, gradharm2;
  double vec[3],dsmooth[3],f1a[3],f1b[3],f2a[3];

  /*fprintf(stdout,"\nEntering calc_endangledist_pair_force:\n");*/

  distwallmin = calc_pwdist(p1, iaparams, &clconstr);
  distmx = iaparams->p.endangledist.distmax;
  distmn = iaparams->p.endangledist.distmin;

  if (distwallmin < distmx) {
    /* function which goes smoothly from 0 to 1 as z goes from distmax to distmin */
    if (distwallmin < distmn) {
      smooth = 1.0;
      for(i=0;i<3;i++) {
        dsmooth[i] = 0.0;
      }
    } else {
      smooth = 0.5*(cos((distwallmin-distmn)/(distmx-distmn)*PI)+1.0);
      for(i=0;i<3;i++) {
        dsmooth[i] = -0.5*PI/(distmx-distmn)*sin((distwallmin-distmn)/(distmx-distmn)*PI)*constraints[clconstr].c.wal.n[i];
      }
    }
    /* Get vector from particle 1 to particle 2 */
    get_mi_vector(vec, p2->r.p, p1->r.p);
    dist = sqrlen(vec);
    di = 1.0/sqrt(dist);
    /*
    for(j=0;j<3;j++) vec[j] *= di;
    */
    /* Calculate angle that p1-p2 makes with wall */
    phi = calc_pwangle(p1, p2, iaparams, &clconstr);

    sinphi = sin(phi);
    cosphi = cos(phi);
    bend   = iaparams->p.endangledist.bend;
    phieq  = iaparams->p.endangledist.phi0;

    /*
    fprintf(stdout,"  Bead %4d: Cl.wall=%2d  distwallmin=%9.6f  \n",p1->p.identity,clconstr,distwallmin);
    fprintf(stdout,"    vector=(%f %f %f)\n",vec[0],vec[1],vec[2]);
    fprintf(stdout,"pos1=(%f %f %f)  pos2=(%f %f %f)  distwallmin=%9.6f  angle=%9.6f\n",p1->r.p[0],p1->r.p[1],p1->r.p[2],p2->r.p[0],p2->r.p[1],p2->r.p[2],distwallmin,phi*180.0/PI);
    */

#ifdef BOND_ENDANGLEDIST_HARMONIC
    /* Force = -dU/dr_i= k*smooth*(phi-phi0)/sin(phi)(cosphi*vec + n)/|vec| */
    fac_a = bend*(phi-phieq)/sinphi;
    fac_b = 0.5*bend*SQR(phi-phieq);
    for(i=0;i<3;i++) {
      gradharm1 = -1.0*fac_a*(cosphi*vec[i]-constraints[clconstr].c.wal.n[i])*di;
      gradharm2 = -1.0*gradharm1;
      f1a[i] = smooth*gradharm1;
      f1b[i] = dsmooth[i]*fac_b;
      f2a[i] = smooth*gradharm2;
      force1[i] = -1.0*(f1a[i]+f1b[i]);
      force2[i] = -1.0*f2a[i];
    }
    /*
    fprintf(stdout,"    f1=(% 9.6f % 9.6f % 9.6f)  ",f1[0],f1[1],f1[2]);
    fprintf(stdout,"    f2=(% 9.6f % 9.6f % 9.6f)\n",f2[0],f2[1],f2[2]);
    fprintf(stdout," force=(% 9.6f % 9.6f % 9.6f)  ",force1[0],force1[1],force1[2]);
    fprintf(stdout," force2=(% 9.6f % 9.6f % 9.6f)\n",force2[0],force2[1],force2[2]);
    */
  } else if (distwallmin >= distmx) {
    for(i=0;i<3;i++) {
      force1[i] = 0.0;
      force2[i] = 0.0;
    }
  }
#endif

  ONEPART_TRACE(if(p1->p.identity==check_id)                            \
                  fprintf(stderr,"%d: OPT: ENDANGLEDIST f = (%.3e,%.3e,%.3e) " \
                          "with part id=%d at dist %f fac %.3e %.3e\n",      \
                          this_node,p1->f.f[0],p1->f.f[1],p1->f.f[2],   \
                          p2->p.identity,dist,fac_a,fac_b)                   \
                );
  ONEPART_TRACE(if(p2->p.identity==check_id)                            \
                  fprintf(stderr,"%d: OPT: ENDANGLEDIST f = (%.3e,%.3e,%.3e)" \
                          " with part id=%d at dist %f fac %.3e\n",     \
                          this_node,p2->f.f[0],p2->f.f[1],p2->f.f[2],   \
                          p1->p.identity,dist,fac_a,fac_b));

  return 0;
}

int endangledist_pair_energy(Particle *p1, Particle *p2, Bonded_ia_parameters *iaparams, double dx[3], double *_energy)
{
  int clconstr=0;
  double bend=0.0,phieq=0.0,phi=0.0;
  double distwallmin=0.0, distmx, distmn, smooth;

  /*fprintf(stdout,"Entering endangledist_pair_energy\n");*/

  distwallmin = calc_pwdist(p1, iaparams, &clconstr);
  /* fprintf(stdout,"clconstr=%d\n",clconstr);*/
  /* fprintf(stdout,"Minimum particle-wall distance=%f\n",distwallmin);*/
  distmx = iaparams->p.endangledist.distmax;
  distmn = iaparams->p.endangledist.distmin;

#ifdef BOND_ENDANGLEDIST_HARMONIC
  if (distwallmin < distmx) {
    /* function which goes smoothly from 0 to 1 as z goes from distmax to distmin */
    if (distwallmin < distmn) {
      smooth = 1.0;
    } else {
      smooth = 0.5*(cos((distwallmin-distmn)/(distmx-distmn)*PI)+1);
    }
    /* Calculate angle that p1-p2 makes with wall */
    phi = calc_pwangle(p1, p2, iaparams, &clconstr);
    /*fprintf(stdout,"clconstr=%d smooth=%f\n",clconstr,smooth);*/
    bend   = iaparams->p.endangledist.bend;
    phieq  = iaparams->p.endangledist.phi0;
    *_energy = 0.5*bend*smooth*SQR(phi - phieq);
  } else if (distwallmin >= distmx) {
    *_energy = 0.0;
  }
#endif

  return 0;
}

#endif


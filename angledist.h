// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2006; all rights reserved unless otherwise stated.
#ifndef ANGLEDIST_H
#define ANGLEDIST_H
/** \file angledist.h
 *  Routines to calculate the angle and distance dependent (from a constraint) energy or/and and force
 *  for a particle triple.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:stuehn@mpip-mainz.mpg.de">Torsten Stuehn</a>
*/

#ifdef BOND_ANGLEDIST

#include "utils.h"

/************************************************************/

/** set parameters for the angledist potential. The type of the angledist potential
    is chosen via config.h and cannot be changed at runtime.
**/
MDINLINE int angledist_set_params(int bond_type, double bend, double phimin, double distmin, double phimax, double distmax)
{
  if(bond_type < 0)
    return TCL_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].p.angledist.bend = bend;
  bonded_ia_params[bond_type].p.angledist.phimin = phimin;
  bonded_ia_params[bond_type].p.angledist.distmin = distmin;
  bonded_ia_params[bond_type].p.angledist.phimax = phimax;
  bonded_ia_params[bond_type].p.angledist.distmax = distmax;
#ifdef BOND_ANGLEDIST_COSINE
#error angledist not implemented for BOND_ANGLEDIST_COSINE
  /*  bonded_ia_params[bond_type].p.angledist.cos_phi0 = cos(phi0);
      bonded_ia_params[bond_type].p.angledist.sin_phi0 = sin(phi0);*/
#endif
#ifdef BOND_ANGLEDIST_COSSQUARE
#error angledist not implemented for BOND_ANGLEDIST_COSSQUARE
  /*  bonded_ia_params[bond_type].p.angledist.cos_phi0 = cos(phi0);*/
#endif
  bonded_ia_params[bond_type].type = BONDED_IA_ANGLEDIST;
  bonded_ia_params[bond_type].num = 2;

  /* 
  printf ("bond_type=%d\n",bond_type);
  printf ("bonded_ia_params[%d].p.angledist.bend=%f\n",bond_type,bonded_ia_params[bond_type].p.angledist.bend);
  printf ("bonded_ia_params[%d].p.angledist.phimin=%f\n",bond_type,bonded_ia_params[bond_type].p.angledist.phimin);
  printf ("bonded_ia_params[%d].p.angledist.distmin=%f\n",bond_type,bonded_ia_params[bond_type].p.angledist.distmin);
  printf ("bonded_ia_params[%d].p.angledist.phimax=%f\n",bond_type,bonded_ia_params[bond_type].p.angledist.phimax);
  printf ("bonded_ia_params[%d].p.angledist.distmax=%f\n\n",bond_type,bonded_ia_params[bond_type].p.angledist.distmax);
  */

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1); 

  return TCL_OK;
}

/// parse parameters for the angle potential
MDINLINE int inter_parse_angledist(Tcl_Interp *interp, int bond_type, int argc, char **argv)
{
  double bend, phimin, distmin, phimax, distmax;

  if (argc != 6) {
    Tcl_AppendResult(interp, "angledist needs 5 parameters: "
		     "<bend> <phimin> <distmin> <phimax> <distmax>", (char *) NULL);
    printf ("argc=%d\n",argc);
    return (TCL_ERROR);
  }

  if (! ARG_IS_D(1, bend)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<bend> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(2, phimin)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
                     "<phimin> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(3, distmin)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<distmin> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(4, phimax)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
                     "<phimax> ", (char *) NULL);
    return TCL_ERROR;
  }

  if (! ARG_IS_D(5, distmax)) {
    Tcl_AppendResult(interp, "angledist needs a DOUBLE parameter: "
		     "<distmax> ", (char *) NULL);
    return TCL_ERROR;
  }


  CHECK_VALUE(angledist_set_params(bond_type, bend, phimin, distmin, phimax, distmax), "bond type must be nonnegative");
}

/** Computes the four body angle and distance interaction force and adds this
    force to the particle forces (see \ref #inter). 
    @param p_mid     Pointer to second/middle particle.
    @param p_left    Pointer to first/left particle.
    @param p_right   Pointer to third/right particle.
    @param iaparams  bond type number of the angledist interaction (see \ref #inter).
    @param force1 returns force of particle 1
    @param force2 returns force of particle 2
    @return 0
*/

/***
MDINLINE double calculate_wall_angledist(Particle *p_mid, double ppos[3], Particle *c_p, Constraint_wall *c, double *dist, double *vec)
{
  int i;
  double pwdist;

  printf ("\nIn function calculate_wall_angledist in angledist.h\n");

  *dist = -c->d;
  for(i=0;i<3;i++) *dist += ppos[i]*c->n[i];
  
  for(i=0;i<3;i++) vec[i] = c->n[i] * *dist;
  
  printf ("particle position (%f,%f,%f)\n",ppos[0],ppos[1],ppos[2]);
  printf ("Wall normal (%f,%f,%f)\n",c.n[0],c.n[1],c.n[2]);
  printf ("Distance=%f\n",c.d);

  return pwdist;
}
***/

 /*                                  Constraint_wall *c,*/

MDINLINE double calc_angdist_param(Particle *p_mid, Particle *p_left, Particle *p_right, 
                       Bonded_ia_parameters *iaparams)
{
  double cosine=0.0, vec1[3], vec2[3], d1i=0.0, d2i=0.0, dist1=0.0, dist2=0.0, fac=0.0, phi0=0.0;
  double pwdist=0.0, pwdist0=0.0, pwdist1=0.0, normal, force[3], folded_pos[3], phimn=0.0, distmn=0.0, phimx=0.0, distmx=0.0, drange=0.0;
  Constraint_wall wall;
  int j, k;
  int img[3];

  /*
  FILE *fp;
  if ((fp = fopen("angledist.dat","a+")) == NULL){
    fprintf(stderr,"Can't open angledist.dat\n");
    exit(EXIT_FAILURE);
  }
  */
  cosine=0.0;
  /* vector from p_left to p_mid */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist1 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist1);
  /*printf ("v1[0]=%f  v1[1]=%f  v1[2]=%f dist1=%f\n",vec1[0],vec1[1],vec1[2], dist1);*/
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_mid to p_right */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  /*printf ("v2[0]=%f  v2[1]=%f  v2[2]=%f dist2=%f\n",vec2[0],vec2[1],vec2[2], dist2);*/
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* vectors are normalised so cosine is just cos(angle_between_vec1_and_vec2) */
  cosine = scalar(vec1, vec2);
  if ( cosine >  TINY_COS_VALUE)  cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  /* printf ("cosine=%f  angles=%f %f (deg)  \n",cosine,acos(cosine)*180/PI,acos(-cosine)*180/PI); */
  fac    = iaparams->p.angledist.bend; /* spring constant from .tcl file */
  phimn  = iaparams->p.angledist.phimin;
  distmn = iaparams->p.angledist.distmin;
  phimx  = iaparams->p.angledist.phimax;
  distmx = iaparams->p.angledist.distmax;
  /* printf ("fac=%f  phimin=%f  distmin=%f  phimax=%f  distmax=%f\n",fac,phimn,distmn,phimx,distmx); */

  /* folds coordinates of p_mid into original box */
  memcpy(folded_pos, p_mid->r.p, 3*sizeof(double));
  memcpy(img, p_mid->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);
  /* printf ("Particle %d folded to (%f,%f,%f)\n",p_mid->p.identity, folded_pos[0],folded_pos[1],folded_pos[2]); 
   */
  /*
  printf ("Box dimensions:  "); 
  for (j=0;j<3;j++) printf ("%f  ",box_l[j]);
  printf ("\nFolded_pos=%f  %f  %f\n",folded_pos[0],folded_pos[1],folded_pos[2]); 
  */

  /* Calculates distance between p_mid and constraint */
  /*fprintf (fp,"calc_angdist_param:  %f %f:  ",distmn,distmx);*/
  for(k=0;k<n_constraints;k++) {
    /*printf ("\nWall %d of %d:  ",k+1,n_constraints);*/
    /*printf ("(%3.2f %3.2f %3.2f), %f, box_lz=%f\n",wall.n[0],wall.n[1],wall.n[2],wall.d,box_l[2]);*/
    for (j=0; j<3; j++) {
      force[j] = 0;
    }
    /* printf ("\nConstraint number %d:",k+1);*/
    switch(constraints[k].type) {
    case CONSTRAINT_WAL: 

      /* dist is distance of wall from origin */
      wall=constraints[k].c.wal;
      /*printf ("wall0=%f,wall1=%f,dist=%f\n",constraints[0].c.wal.d,constraints[1].c.wal.d,fabs(constraints[1].c.wal.d)-fabs(constraints[0].c.wal.d));*/

      /* check that constraint vector is normalised */
      normal=0.0;
      for(j=0;j<3;j++) normal += wall.n[j] * wall.n[j];
      if (sqrt(normal) != 1.0) {
        for(j=0;j<3;j++) wall.n[j]=wall.n[j]/normal;
      }
      /* If distance from origin is outside box range fold constraints into box
         However, wall.d is signed so that normal * distance is positive and it
         should not be necessary to fold wall into box. NB: This has not been tested.
      for (j=0;j<3;j++) nwall.n[j] = wall.n[j] * wall.d;
      for (j=0;j<3;j++) {
        while (nwall.n[j] < 0.0 ) nwall.n[j] += box_l[j];
        while (nwall.n[j] > box_l[j] ) nwall.n[j] -= box_l[j];
      }
      nwall.d=0.0;
      for (j=0;j<3;j++) nwall.d += nwall.n[j] * nwall.n[j];
      nwall.d = sqrt(nwall.d);
      printf ("Wall normal   = %f %f %f\n",wall.n[0],wall.n[1],wall.n[2]);
      printf ("Wall distance = %f, normalised %f\n",wall.d,nwall.d);
      */

      /* pwdist is distance of wall from p_mid */
      pwdist0=-1.0 * constraints[0].c.wal.d;
      for(j=0;j<3;j++) {
        pwdist0 += folded_pos[j] * constraints[0].c.wal.n[j];
      }
      pwdist1=-1.0 * constraints[1].c.wal.d;
      for(j=0;j<3;j++) {
        pwdist1 += folded_pos[j] * constraints[1].c.wal.n[j];
      }
      if (pwdist0 <= pwdist1) {
        pwdist = pwdist0;
      }
      else {
        pwdist = pwdist1;
      }
      /*printf ("pos=(%f,%f,%f) dist=%f ",folded_pos[0],folded_pos[1],folded_pos[2],pwdist);*/
      /*printf ("pwdist0=%f pwdist1=%f pwdist=%f\n",pwdist0,pwdist1,pwdist);*/

      /*get phi0(z)*/
      if (pwdist <= distmn) {
        phi0 = phimn;
      }
      else if (pwdist >= distmx && pwdist <= box_l[2]-wall.d-distmx) {
          phi0 = phimx;
      }
      else {
        drange = (pwdist-distmn)*PI/(distmx-distmn);
        phi0 = ((cos(drange-PI)+1.0)*(phimx-phimn))*0.5+phimn;
      }
      break;
    }
  }
  /*
  fprintf(fp,"%f %f  ",pwdist,phi0/PI*180.0);
  fclose(fp);
  */
  return phi0;
}


MDINLINE int calc_angledist_force(Particle *p_mid, Particle *p_left, Particle *p_right,
                                  Bonded_ia_parameters *iaparams, double force1[3], double force2[3])
{
  double cosine, vec1[3], vec2[3], d1i=0.0, d2i=0.0, fac=0.0, f1=0.0, f2=0.0, phi0=0.0;
  int j;

  /* NOTE The angledist is ONLY implemented for the HARMONIC case */
  /* printf ("\nIn function calc_angledist_force in angledist.h\n"); */
  phi0=calc_angdist_param(p_mid, p_left, p_right, iaparams);
  /*
  FILE *fp;
  if ((fp = fopen("angledist.dat","a+")) == NULL){
    fprintf(stderr,"Can't open angledist.dat\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp,"\nIn calc_angledist_force:  phi0=%f\n",phi0*180.0/PI);
  fclose(fp);
  */
#ifdef BOND_ANGLEDIST_HARMONIC
  {
    double phi=0.0,sinphi=0.0;
    if ( cosine >  TINY_COS_VALUE) cosine =  TINY_COS_VALUE;
    if ( cosine < -TINY_COS_VALUE) cosine = -TINY_COS_VALUE;
    phi =  acos(-cosine);
    sinphi = sin(phi);
    if ( sinphi < TINY_SIN_VALUE ) sinphi = TINY_SIN_VALUE;
    fac *= (phi - phi0)/sinphi;
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

/** Computes the three body angle interaction energy (see \ref #inter, \ref #analyze). 
    @param p_mid        Pointer to second/middle particle.
    @param p_left       Pointer to first particle.
    @param p_right      Pointer to third particle.
    @param iaparams  bond type number of the angle interaction (see \ref #inter).
    @param _energy   return energy pointer.
    @return 0.
*/



MDINLINE int angledist_energy(Particle *p_mid, Particle *p_left, Particle *p_right, 
       			      Bonded_ia_parameters *iaparams, double *_energy)
{
  double cosine, vec1[3], vec2[3],  d1i, d2i, dist1, dist2;
  int j;
  double phimn=0.0, phimx=0.0, distmn=0.0, distmx=0.0, folded_pos[3], pwdist=0.0, phi0=PI, drange, force[3], normal;
  Constraint_wall wall;
  int k, img[3];

  phi0=calc_angdist_param(p_mid, p_left, p_right, iaparams);
  /* 
  FILE *fp;
  if ((fp = fopen("angledist.dat","a+")) == NULL){
    fprintf(stderr,"Can't open angledist.dat\n");
    exit(EXIT_FAILURE);
  }
  fprintf(fp,"\nIn angledist_energy:  phi0=%f\n",phi0*180.0/PI);
  fclose(fp);
  */

  cosine=0.0;
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

  /* get phi0 as a function of distance from wall constraint */
  phimn  = iaparams->p.angledist.phimin;
  distmn = iaparams->p.angledist.distmin;
  phimx  = iaparams->p.angledist.phimax;
  distmx = iaparams->p.angledist.distmax;
  memcpy(folded_pos, p_mid->r.p, 3*sizeof(double));
  memcpy(img, p_mid->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);
  for(k=0;k<n_constraints;k++) {
    for (j=0; j<3; j++) {
      force[j] = 0;
    }
    /* printf ("\nConstraint number %d:",k+1);*/
    switch(constraints[k].type) {
    case CONSTRAINT_WAL: 
      wall=constraints[k].c.wal;
      normal=0.0;
      for(j=0;j<3;j++) normal += wall.n[j] * wall.n[j];
      if (sqrt(normal) != 1.0) {
        for(j=0;j<3;j++) wall.n[j]=wall.n[j]/normal;
      }
      pwdist=-1.0 * wall.d;
      for(j=0;j<3;j++) {
        pwdist += folded_pos[j] * wall.n[j];
      }
      /* printf ("P-W distance  = %f\n",pwdist); */
      if (pwdist <= distmn) {
        phi0 = phimn;
        /* printf ("pwdist <= distmn: dist=%f  phi0=%f\n",pwdist, phi0); */
      }
      else if (pwdist >= distmx) {
        phi0 = phimx;
        /* printf ("pwdist >= distmx: dist=%f  phi0=%f\n",pwdist, phi0); */
      }
      else {
        drange = (pwdist-distmn)*PI/(distmx-distmn);
        phi0 = ((cos(drange-PI)+1.0)*(phimx-phimn))*0.5+phimn;
        /* printf ("distmn < pwdist < distmx: dist=%f  phi0=%f\n",pwdist, phi0); */
      }
      break;
    }
  }


#ifdef BOND_ANGLEDIST_HARMONIC
  {
    double phi;
    phi =  acos(-cosine);
    *_energy = 0.5*iaparams->p.angledist.bend*SQR(phi - phi0);
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

#endif /* BOND_ANGLEDIST */
#endif /* ANGLEDIST_H */

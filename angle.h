// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2003; all rights reserved unless otherwise stated.
#ifndef ANGLE_H
#define ANGLE_H
/** \file angle.h
 *  Routines to calculate the angle energy or/and and force 
 *  for a particle pair.
 *  \ref forces.c
 *
 *  <b>Responsible:</b>
 *  <a href="mailto:limbach@mpip-mainz.mpg.de">Hanjo</a>
*/


#define TINY_SIN_VALUE 1e-9
#define TINY_COS_VALUE 0.9999999999

/************************************************************/

/** Computes the three body angle interaction force and adds this
    force to the particle forces (see \ref #inter). 
    @param p_mid     Pointer to second/middle particle.
    @param p_left    Pointer to first/left particle.
    @param p_right   Pointer to third/right particle.
    @param type_num  bond type number of the angle interaction (see \ref #inter).
*/
MDINLINE void add_angle_force(Particle *p_mid, Particle *p_left, Particle *p_right, int type_num)
{
  double cosine, vec1[3], vec2[3], d1i, d2i, dist2,  fac, f1=0.0, f2=0.0;
  int j;

  cosine=0.0;
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
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  fac    = bonded_ia_params[type_num].p.angle.bend;
#ifdef BOND_ANGLE_HARMONIC
  {
    double phi,sinphi;
    phi =  acos(-cosine);
    sinphi = sin(phi);
    if ( sinphi < TINY_SIN_VALUE ) sinphi = TINY_SIN_VALUE;
    fac *= (phi-bonded_ia_params[type_num].p.angle.phi0)/sinphi;
  }
#endif
#ifdef BOND_ANGLE_COSINE
  if ( cosine >  TINY_COS_VALUE ) cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  fac *= bonded_ia_params[type_num].p.angle.sin_phi0 * (cosine/sqrt(1-SQR(cosine))) + bonded_ia_params[type_num].p.angle.cos_phi0;
#endif
#ifdef BOND_ANGLE_COSSQUARE
  fac *= bonded_ia_params[type_num].p.angle.cos_phi0 + cosine;
#endif
  /* apply bend forces */
  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j]) * d2i;
    p_left->f.f[j]  -= f1;
    p_mid->f.f[j]   += (f1-f2);
    p_right->f.f[j] += f2;
  }
}

/** Computes the three body angle interaction energy (see \ref #inter, \ref #analyze). 
    @param p_mid        Pointer to first particle.
    @param p_left        Pointer to second/middle particle.
    @param p_right        Pointer to third particle.
    @param type_num  bond type number of the angle interaction (see \ref #inter).
    @return energy   energy.
*/
MDINLINE double angle_energy(Particle *p_mid, Particle *p_left, Particle *p_right, int type_num)
{
  double cosine, vec1[3], vec2[3],  d1i, d2i, dist2;
  int j;

  cosine=0.0;
  /* vector from p_mid to p_left */
  get_mi_vector(vec1, p_mid->r.p, p_left->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p_right to p_mid */
  get_mi_vector(vec2, p_right->r.p, p_mid->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  /* bond angle energy */
#ifdef BOND_ANGLE_HARMONIC
  {
    double phi;
    phi =  acos(-cosine);
    return 0.5*bonded_ia_params[type_num].p.angle.bend*SQR(phi-bonded_ia_params[type_num].p.angle.phi0);
  }
#endif
#ifdef BOND_ANGLE_COSINE
  {
    return bonded_ia_params[type_num].p.angle.bend*(cosine*bonded_ia_params[type_num].p.angle.cos_phi0 - sqrt(1-SQR(cosine))*bonded_ia_params[type_num].p.angle.sin_phi0+1);
  }
#endif
#ifdef BOND_ANGLE_COSSQUARE
  return 0.5*bonded_ia_params[type_num].p.angle.bend*SQR(cosine+bonded_ia_params[type_num].p.angle.cos_phi0);
#endif
}


#endif

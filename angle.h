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

/************************************************************/

/** Computes the three body angle interaction force and adds this
    force to the particle forces (see \ref #inter). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param p3        Pointer to third particle.
    @param type_num  bond type number of the angle interaction (see \ref #inter).
*/
MDINLINE void add_angle_force(Particle *p1, Particle *p2, Particle *p3, int type_num)
{
  double cosine, vec1[3], vec2[3], d1i, d2i, dist2, f1, f2;
  int j;

  cosine=0.0;
  /* vector from p1 to p2 */
  get_mi_vector(vec1, p2->r.p, p1->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p1 to p3 */
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  /* apply bend forces */
  for(j=0;j<3;j++) {
    f1 = bonded_ia_params[type_num].p.angle.bend * (vec2[j] - cosine * vec1[j]) * d1i;
    f2 = bonded_ia_params[type_num].p.angle.bend * (vec1[j] - cosine * vec2[j]) * d2i;
    p2->f.f[j] -= f1;
    p1->f.f[j] += (f1+f2);
    p3->f.f[j] -= f2;
  }
}

/** Computes the three body angle interaction energy (see \ref #inter, \ref #analyze). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param p3        Pointer to third particle.
    @param type_num  bond type number of the angle interaction (see \ref #inter).
    @return energy   energy.
*/
MDINLINE double angle_energy(Particle *p1, Particle *p2, Particle *p3, int type_num)
{
  double cosine, vec1[3], vec2[3], d1i, d2i, dist2;
  int j;

  cosine=0.0;
  /* vector from p1 to p2 */
  get_mi_vector(vec1, p2->r.p, p1->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p3 to p1 */
  get_mi_vector(vec2, p1->r.p, p3->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  /* bond bond angle energy */
  return bonded_ia_params[type_num].p.angle.bend * ( 1 - cosine );
}


#endif

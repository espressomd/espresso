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

MDINLINE void add_angle_pair_force(Particle *p1, Particle *p2, Particle *p3, int type_num)
{
  double cosine, vec1[3], vec2[3], d1i, d2i, dist2, f1, f2;
  int j;

  cosine=0.0;
  /* vector from p2 to p1 */
  for(j=0;j<3;j++) vec1[j] = p2->r.p[j] - p1->r.p[j];
  dist2 = SQR(vec1[0]) + SQR(vec1[1]) + SQR(vec1[2]);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p3 to p1 */
  for(j=0;j<3;j++) vec2[j] = p3->r.p[j] - p1->r.p[j];
  dist2 = SQR(vec2[0]) + SQR(vec2[1]) + SQR(vec2[2]);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  for(j=0;j<3;j++) cosine += vec1[j] * vec2[j];
  /* apply bend forces */
  for(j=0;j<3;j++) {
    f1 = bonded_ia_params[type_num].p.angle.bend * (vec2[j] - cosine * vec1[j]) * d1i;
    f2 = bonded_ia_params[type_num].p.angle.bend * (vec1[j] - cosine * vec2[j]) * d2i;
    p2->f[j] -= f1;
    p1->f[j] += (f1+f2);
    p3->f[j] -= f2;
  }
}

#endif

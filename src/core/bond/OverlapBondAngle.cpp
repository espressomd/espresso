#include "OverlapBondAngle.hpp"
#include "grid.hpp" //get_mi_vector
#include "interaction_data.hpp"

//force *
/** Computes the three body overlapped angle interaction force.
    Adds this force to the particle forces in forces.hpp (see \ref
   tclcommand_inter).
    @param p1     Pointer to second/middle particle.
    @param p2    Pointer to first/left particle.
    @param p3   Pointer to third/right particle.
    @param iaparams  bond type number of the angle interaction (see \ref
   tclcommand_inter).
    @param force returns force of particle 1
    @param force2 returns force of particle 2
    @return 0
    Needs feature OVERLAPPED compiled in (see \ref config.hpp).
*/
int Bond::OverlapBondAngle::calc_bonded_three_particle_force(Particle *p1, Particle *p2, 
							    Particle *p3, double force[3], 
							    double force2[3]) const
{

  double cosine, vec1[3], vec2[3], d1i, d2i, dist2, fac, f1 = 0.0, f2 = 0.0;
  int j;

  int ig;
  double ev;
  double prob = 0.;
  double Apart = 0.;

  cosine = 0.0;
  /* vector from p2 to p1 */
  get_mi_vector(vec1, p1->r.p, p2->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for (j = 0; j < 3; j++)
    vec1[j] *= d1i;
  /* vector from p1 to p3 */
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for (j = 0; j < 3; j++)
    vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  /* Notice: cosine = - costheta */
  cosine = scalar(vec1, vec2);

  if (cosine > TINY_COS_VALUE)
    cosine = TINY_COS_VALUE;
  if (cosine < -TINY_COS_VALUE)
    cosine = -TINY_COS_VALUE;

  /* compute fac = - dU/d(costheta) */
  /* prob = sum_(i=1,N) [ a_i * exp(-(costheta -b_i)^2/c_i^2))] */
  /* Apart = sum_(i=1,N) [ a_i * exp(-(costheta -b_i)^2/c_i^2)) * (costheta-b_i)
   * / c_i^2] */
  /* fac = -2.0 * ( Apart / prob ) */

  for (ig = 0; ig < m_noverlaps; ig++) {
    ev = (-cosine - m_para_b[ig]) /
         m_para_c[ig];
    prob = prob + m_para_a[ig] * exp(-1.0 * ev * ev);
    Apart =
        Apart +
        m_para_a[ig] * exp(-1.0 * ev * ev) *
            (-cosine - m_para_b[ig]) /
            (m_para_c[ig] * m_para_c[ig]);
  }
  fac = -2. * (Apart / prob);

  /* compute force */
  for (j = 0; j < 3; j++) {
    f1 = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    f2 = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    force[j] = (f1 - f2);
    force2[j] = -f1;
  }
  return 0;

}

//energy *
/** Computes the three body overlapped angle interaction energy (see \ref
   tclcommand_inter, \ref tclcommand_analyze).
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param p3        Pointer to third particle.
    @param iaparams  bond type number of the angle interaction (see \ref
   tclcommand_inter).
    @param _energy   return energy pointer.
    @return 0.
    Needs feature OVERLAPPED compiled in (see \ref config.hpp).
*/
int Bond::OverlapBondAngle::calc_bonded_three_particle_energy(Particle *p1, Particle *p2, 
							     Particle *p3, double *_energy) const
{

  double cosine, vec1[3], vec2[3], d1i, d2i, dist2;
  int j;

  int ig;
  double ev;
  double prob = 0.;

  cosine = 0.0;
  /* vector from p1 to p2 */
  get_mi_vector(vec1, p1->r.p, p2->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for (j = 0; j < 3; j++)
    vec1[j] *= d1i;
  /* vector from p3 to p1 */
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for (j = 0; j < 3; j++)
    vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  /* Notice: cosine = - costheta */
  cosine = scalar(vec1, vec2);

  if (cosine > TINY_COS_VALUE)
    cosine = TINY_COS_VALUE;
  if (cosine < -TINY_COS_VALUE)
    cosine = -TINY_COS_VALUE;

  /* compute overlapped cosangl energy */
  /*prob = sum_(i=1,N) [ a_i * exp(-(costheta -b_i)^2/c_i^2))] */
  /*pot  = -1.0 * Log { prob } */
  for (ig = 0; ig < m_noverlaps; ig++) {
    ev = (-cosine - m_para_b[ig]) /
         m_para_c[ig];
    prob = prob + m_para_a[ig] * exp(-1.0 * ev * ev);
  }
  *_energy = -log(prob);

  return 0;

}

#include "AngleCosine.hpp"
#include "grid.hpp" // get_mi_vector
#include "interaction_data.hpp"

/************************************************************/

/** Computes the three body angle interaction force and adds this
    force to the particle forces (see \ref tclcommand_inter). 
    @param p1     Pointer to second/middle particle.
    @param p2    Pointer to first/left particle.
    @param p3   Pointer to third/right particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param force1 returns force of particle 1
    @param force2 returns force of particle 2
    @return 0
*/
int Bond::AngleCosine::calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3], double force2[3]) const {

  double cosine, vec1[3], vec2[3], d1i, d2i, dist2,  fac, f1=0.0, f2=0.0;
  int j;

  cosine=0.0;
  /* vector from p2 to p1 */
  get_mi_vector(vec1, p1->r.p, p2->r.p);
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
  fac    = m_bend;

  if ( cosine >  TINY_COS_VALUE ) cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  fac *= m_sin_phi0 * (cosine/sqrt(1-Utils::sqr(cosine))) + m_cos_phi0;
  
  for(j=0;j<3;j++) {
    f1               = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    f2               = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    force[j] = (f1-f2);
    force2[j] = -f1;
  }
  return 0;

}


/** Computes the three body angle interaction energy (see \ref tclcommand_inter, \ref tclcommand_analyze). 
    @param p1        Pointer to first particle.
    @param p2        Pointer to second/middle particle.
    @param p3        Pointer to third particle.
    @param iaparams  bond type number of the angle interaction (see \ref tclcommand_inter).
    @param _energy   return energy pointer.
    @return 0.
*/
int Bond::AngleCosine::calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
							  double *_energy) const {
  double cosine, vec1[3], vec2[3],  d1i, d2i, dist2;
  int j;

  cosine=0.0;
  /* vector from p1 to p2 */
  get_mi_vector(vec1, p1->r.p, p2->r.p);
  dist2 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec1[j] *= d1i;
  /* vector from p3 to p1 */
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for(j=0;j<3;j++) vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  if ( cosine >  TINY_COS_VALUE)  cosine = TINY_COS_VALUE;
  if ( cosine < -TINY_COS_VALUE)  cosine = -TINY_COS_VALUE;
  /* bond angle energy */

  *_energy = m_bend*(cosine*m_cos_phi0 - sqrt(1-Utils::sqr(cosine))*m_sin_phi0+1);

  return 0;
}

/* The force on each particle due to a three-body bonded potential
   is computed. */
int Bond::AngleCosine::calc_3body_forces(Particle *p_mid, Particle *p_left,
					  Particle *p_right, double force1[3], 
					  double force2[3], double force3[3]) const 
{


  int j;
  double pot_dep;
  double cos_phi;
  double sin_phi;
  double vec31[3];
  double vec21[3];
  double vec12[3]; // espresso convention
  double vec21_sqr;
  double vec31_sqr;
  double vec21_magn;
  double vec31_magn;
  double fj[3];
  double fk[3];
  double fac;

  get_mi_vector(vec12, p_mid->r.p, p_left->r.p);
  for(j = 0; j < 3; j++)
    vec21[j] = -vec12[j];

  get_mi_vector(vec31, p_right->r.p, p_mid->r.p);
  vec21_sqr = sqrlen(vec21);
  vec21_magn = sqrt(vec21_sqr);
  vec31_sqr = sqrlen(vec31);
  vec31_magn = sqrt(vec31_sqr);
  cos_phi = scalar(vec21, vec31) / (vec21_magn * vec31_magn);
  sin_phi = sqrt(1.0 - Utils::sqr(cos_phi));

  /* uncomment this block if interested in the angle 
  if(cos_phi < -1.0) cos_phi = -TINY_COS_VALUE;
  if(cos_phi >  1.0) cos_phi =  TINY_COS_VALUE; 
  phi = acos(cos_phi);
  */

  {
    double K, sin_phi0, cos_phi0;
    K = m_bend;
    sin_phi0 = m_sin_phi0;
    cos_phi0 = m_cos_phi0;

    // potential dependent term [dU/dphi = K * sin(phi - phi0)]
    // trig identity: sin(a - b) = sin(a)cos(b) - cos(a)sin(b) 
    pot_dep = K * (sin_phi * cos_phi0 - cos_phi * sin_phi0);
  }

  fac = pot_dep / sin_phi;

  for(j = 0; j < 3; j++) {
    fj[j] = vec31[j] / (vec21_magn * vec31_magn) - cos_phi * vec21[j] / vec21_sqr;
    fk[j] = vec21[j] / (vec21_magn * vec31_magn) - cos_phi * vec31[j] / vec31_sqr;
  }

  // note that F1 = -(F2 + F3)
  for(j = 0; j < 3; j++) {
    force1[j] = force1[j] - fac * (fj[j] + fk[j]);
    force2[j] = force2[j] + fac * fj[j];
    force3[j] = force3[j] + fac * fk[j];
  }

  return 0;

}

boost::any Bond::AngleCosine::get_bond_parameters_from_bond() const
{

  Angle_cosine_bond_parameters params = {m_bend, m_phi0, m_cos_phi0, m_sin_phi0};
  return boost::any(params);
  
}

#include "TabulatedBondAngle.hpp"
#include "grid.hpp" //get_mi_vector

//force *
/** Calculate a tabulated bond angle force with number type_num (see
    \ref Bonded_ia_parameters) between particles p2, p1 and
    p3 and add it to the particle forces. The force on p2 and
    p3 acts perpendicular to the connecting vector between the
    particle and p1 and in the plane defined by the three
    particles. The force on the middle particle balances the other two
    forces. The forces are scaled with the invers length of the
    connecting vectors. It is assumed that the potential is tabulated
    for all angles between 0 and Pi. 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
int Bond::TabulatedBondAngle::calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3], double force2[3]) const
{

  double cosine, phi, invsinphi, vec1[3], vec2[3], d1i, d2i, dist2, fac,
      f1 = 0.0, f2 = 0.0;
  int j;

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
  cosine = scalar(vec1, vec2);
#ifdef TABANGLEMINUS
  phi = acos(-cosine);
#else
  phi = acos(cosine);
#endif
  invsinphi = sin(phi);
  if (invsinphi < TINY_SIN_VALUE)
    invsinphi = TINY_SIN_VALUE;
  invsinphi = 1.0 / invsinphi;
  /* look up force factor */
  fac = m_tab_pot.force(phi);
  /* apply bend forces */
  for (j = 0; j < 3; j++) {
    f1 = fac * (cosine * vec1[j] - vec2[j]) * invsinphi * d1i;
    f2 = fac * (cosine * vec2[j] - vec1[j]) * invsinphi * d2i;
    force[j] = (f1 - f2);
    force2[j] = -f1;
  }

  return 0;

}

//energy *
/** Calculate and return tabulated bond angle energy with number
    type_num (see \ref Bonded_ia_parameters) between particles p2,
    p1 and p3. It is assumed that the potential is tabulated
    for all angles between 0 and Pi. 
    Needs feature TABULATED compiled in (see \ref config.hpp). */
int Bond::TabulatedBondAngle::calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
				     double *_energy) const
{

  double phi, vec1[3], vec2[3], vl1, vl2;

  /* vector from p1 to p2 */
  get_mi_vector(vec1, p1->r.p, p2->r.p);
  vl1 = sqrt(sqrlen(vec1));
  /* vector from p3 to p1 */
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  vl2 = sqrt(sqrlen(vec2));
/* calculate phi */
#ifdef TABANGLEMINUS
  phi = acos(-scalar(vec1, vec2) / (vl1 * vl2));
#else
  phi = acos(scalar(vec1, vec2) / (vl1 * vl2));
#endif

  *_energy = m_tab_pot.energy(phi);

  return 0;

}

// for pressure.hpp
/* The force on each particle due to a three-body bonded tabulated
   potential is computed. */
int Bond::TabulatedBondAngle::calc_3body_forces(Particle *p_mid, Particle *p_left,
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
  double phi, dU; // bond angle and d/dphi of U(phi)

  get_mi_vector(vec12, p_mid->r.p, p_left->r.p);
  for (j = 0; j < 3; j++)
    vec21[j] = -vec12[j];

  get_mi_vector(vec31, p_right->r.p, p_mid->r.p);
  vec21_sqr = sqrlen(vec21);
  vec21_magn = sqrt(vec21_sqr);
  vec31_sqr = sqrlen(vec31);
  vec31_magn = sqrt(vec31_sqr);
  cos_phi = scalar(vec21, vec31) / (vec21_magn * vec31_magn);
  sin_phi = sqrt(1.0 - SQR(cos_phi));

  if (cos_phi < -1.0)
    cos_phi = -TINY_COS_VALUE;
  if (cos_phi > 1.0)
    cos_phi = TINY_COS_VALUE;
#ifdef TABANGLEMINUS
  phi = acos(-cos_phi);
#else
  phi = acos(cos_phi);
#endif

  dU = m_tab_pot.force(phi);

  // potential dependent term (dU/dphi * 1 / sin(phi))
  pot_dep = dU / sin_phi;

  for (j = 0; j < 3; j++) {
    fj[j] =
        vec31[j] / (vec21_magn * vec31_magn) - cos_phi * vec21[j] / vec21_sqr;
    fk[j] =
        vec21[j] / (vec21_magn * vec31_magn) - cos_phi * vec31[j] / vec31_sqr;
  }

  // note that F1 = -(F2 + F3) in analytical case
  for (j = 0; j < 3; j++) {
    force1[j] = force1[j] - pot_dep * (fj[j] + fk[j]);
    force2[j] = force2[j] + pot_dep * fj[j];
    force3[j] = force3[j] + pot_dep * fk[j];
  }
  
  return 0;

}

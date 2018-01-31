#include "AngleDist.hpp"

#ifdef BOND_ANGLEDIST

#include "communication.hpp"
#include "grid.hpp" //for get_mi_vector
#include "constraints.hpp" // for constraints stuff
#include "constraints/ShapeBasedConstraint.hpp"
#include "shapes/Wall.hpp" // for shapes


int Bond::AngleDist::calc_bonded_three_particle_force(Particle *p1, Particle *p2, Particle *p3, double force[3], 
				    double force2[3]) const {

  double cosine = 0.0, vec1[3], vec2[3], fac = 0.0, f1 = 0.0, f2 = 0.0;

  /* vector from p2 to p1 */
  get_mi_vector(vec1, p1->r.p, p2->r.p);
  const double dist12 = sqrlen(vec1);
  const double d1i = 1.0 / sqrt(dist12);
  for (int j = 0; j < 3; j++)
    vec1[j] *= d1i;

  /* vector from p1 to p3 */
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  const double dist22 = sqrlen(vec2);
  const double d2i = 1.0 / sqrt(dist22);
  for (int j = 0; j < 3; j++)
    vec2[j] *= d2i;
  /* scalar product of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  fac = m_bend;

#ifdef BOND_ANGLEDIST_HARMONIC
  {
    const double phi0 = calc_angledist_param(p1, p2, p3);

    if (cosine > TINY_COS_VALUE)
      cosine = TINY_COS_VALUE;
    if (cosine < -TINY_COS_VALUE)
      cosine = -TINY_COS_VALUE;
    const double phi = acos(-cosine);

    double sinphi = sin(phi);
    if (sinphi < TINY_SIN_VALUE)
      sinphi = TINY_SIN_VALUE;
    fac *= (phi - phi0) / sinphi;
  }
#endif

  for (int j = 0; j < 3; j++) {
    f1 = fac * (cosine * vec1[j] - vec2[j]) * d1i;
    f2 = fac * (cosine * vec2[j] - vec1[j]) * d2i;

    force[j] = (f1 - f2);
    force2[j] = -f1;
  }
  return 0;

}


int Bond::AngleDist::calc_bonded_three_particle_energy(Particle *p1, Particle *p2, Particle *p3, 
				     double *_energy) const {

  int j;
  double cosine = 0.0, d1i, d2i, dist1, dist2;
  double vec1[3], vec2[3];

  /* vector from p1 to p2 */
  get_mi_vector(vec1, p1->r.p, p2->r.p);
  dist1 = sqrlen(vec1);
  d1i = 1.0 / sqrt(dist1);
  for (j = 0; j < 3; j++)
    vec1[j] *= d1i;
  /* vector from p3 to p1 */
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for (j = 0; j < 3; j++)
    vec2[j] *= d2i;
  /* scalar produvt of vec1 and vec2 */
  cosine = scalar(vec1, vec2);
  if (cosine > TINY_COS_VALUE)
    cosine = TINY_COS_VALUE;
  if (cosine < -TINY_COS_VALUE)
    cosine = -TINY_COS_VALUE;

#ifdef BOND_ANGLEDIST_HARMONIC
  {
    double phi;
    double phi0 = calc_angledist_param(p1, p2, p3);
    phi = acos(-cosine);
    *_energy = 0.5 * m_bend * SQR(phi - phi0);
  }
#endif

  return 0;
}

/** Function to calculate wall distance and phi0(dist).
    Called by \ref calc_angledist_force */
double Bond::AngleDist::calc_angledist_param(Particle *p1, Particle *p2,
			    Particle *p3) const {

  double vec1[3], vec2[3], d2i = 0.0, dist2 = 0.0;

  double normal, folded_pos[3];

  int img[3];

  /* vector from p2 to p1 */
  get_mi_vector(vec1, p1->r.p, p2->r.p);
  const double dist1 = sqrlen(vec1);
  const double d1i = 1.0 / sqrt(dist1);
  for (int j = 0; j < 3; j++)
    vec1[j] *= d1i;
  /* vector from p1 to p3 */
  get_mi_vector(vec2, p3->r.p, p1->r.p);
  dist2 = sqrlen(vec2);
  d2i = 1.0 / sqrt(dist2);
  for (int j = 0; j < 3; j++)
    vec2[j] *= d2i;

  const double phimn = m_phimin;
  const double distmn = m_distmin;
  const double phimx = m_phimax;
  const double distmx = m_distmax;

  /* folds coordinates of p1 into original box */
  memmove(folded_pos, p1->r.p, 3 * sizeof(double));
  memmove(img, p1->l.i, 3 * sizeof(int));
  fold_position(folded_pos, img);

  double pwdistmin = std::numeric_limits<double>::infinity();
  double pwdistmin_d = 0.0;

  for (auto const &c : Constraints::constraints) {
    auto cs = std::dynamic_pointer_cast<const Constraints::ShapeBasedConstraint>(c);
    if (cs) {
      if (dynamic_cast<const Shapes::Wall*>(&(cs->shape()))) {
        const Shapes::Wall* wall = dynamic_cast<const Shapes::Wall*>(&(cs->shape()));
        double dist = -wall->d();
        for (int j = 0; j < 3; j++) {
          dist += folded_pos[j] * (wall->n())[j];
        };

        if (dist < pwdistmin) {
          pwdistmin = dist;
          pwdistmin_d = wall->d();
        };
      };
    };
  };

  /*get phi0(z)*/
  if (pwdistmin <= distmn) {
    return phimn;
  } else if (pwdistmin >= distmx &&
             pwdistmin <= box_l[2] - pwdistmin_d - distmx) {
    return phimx;
  } else {
    const double drange = (pwdistmin - distmn) * PI / (distmx - distmn);
    return ((cos(drange - PI) + 1.0) * (phimx - phimn)) * 0.5 + phimn;
  }

}

#endif //BOND_ANGLE_DIST

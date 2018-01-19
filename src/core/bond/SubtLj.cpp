#include "SubtLj.hpp"
#include "lj.hpp"

int Bond::SubtLj::calc_bonded_pair_force(Particle *p1, Particle *p2,
                                         double dx[3], double force[3]) const {
#ifdef LENNARD_JONES
  double dist = sqrt(sqrlen(dx));
  double lj_force[3] = {0., 0., 0.};

  add_lj_pair_force(p1, p2, get_ia_param(p1->p.type, p2->p.type), dx, dist,
                    lj_force);

  for (int i = 0; i < 3; i++) {
    force[i] -= -lj_force[i];
  }

#endif
  return 0;
}

int Bond::SubtLj::calc_bonded_pair_energy(Particle *p1, Particle *p2,
                                          double dx[3], double *_energy) const {
#ifdef LENNARD_JONES
  double dist = sqrt(sqrlen(dx));

  *_energy =
      -lj_pair_energy(p1, p2, get_ia_param(p1->p.type, p2->p.type), dx, dist);
#endif
  return 0;
}

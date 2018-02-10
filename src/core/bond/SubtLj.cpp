#include "SubtLj.hpp"
#include "lj.hpp"

int Bond::SubtLj::calc_bonded_pair_force(Particle *p1, Particle *p2,
                                         double dx[3], double force[3]) const {
#ifdef LENNARD_JONES
  for (int i = 0; i < 3; i++) {
    dx[i] *= -1;
  };
  add_lj_pair_force(p1, p2, get_ia_param(p1->p.type, p2->p.type), dx,
		    Utils::veclen(dx), force);
#endif
  return 0;
}

int Bond::SubtLj::calc_bonded_pair_energy(Particle *p1, Particle *p2,
                                          double dx[3], double *_energy) const {
#ifdef LENNARD_JONES
  *_energy = -lj_pair_energy(p1, p2, get_ia_param(p1->p.type, p2->p.type), dx, Utils::veclen(dx));
#endif
  return 0;
}

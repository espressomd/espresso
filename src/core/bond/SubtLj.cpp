#include "SubtLj.hpp"
#include "interaction_data.hpp" // for iaparams LJ
#include "lj.hpp" //add_lj_pair_force(), lj_pair_energy()

int Bond::SubtLj::calc_bonded_pair_force(Particle *p1, Particle *p2, double dx[3], double force[3]) 
  const {

  auto ia_params = get_ia_param(p1->p.type, p2->p.type);

  for (int i = 0; i < 3; i++) {
    dx[i] *= -1;
  }

  add_lj_pair_force(p1, p2, ia_params, dx, Utils::veclen(dx), force);

  return ES_OK;

}

int Bond::SubtLj::calc_bonded_pair_energy(Particle *p1, Particle *p2, double dx[3], double *_energy)
  const {

  auto ia_params = get_ia_param(p1->p.type, p2->p.type);

  *_energy = -lj_pair_energy(p1, p2, ia_params, dx, Utils::veclen(dx));
  return ES_OK;

}

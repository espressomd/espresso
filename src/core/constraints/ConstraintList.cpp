
#include "ConstraintList.hpp"
#include "GeometryConstraint.hpp"

#include "energy.hpp"
#include "errorhandling.hpp"

#include <limits>

namespace Constraints {

ConstraintList list;

/** @TODO: Collect total force on constraint */

void ConstraintList::init_forces() {
  for(iterator it = begin(); it != end(); ++it) {
    for(int i = 0; i < 3; i++)
      (*it)->total_force[i] = 0.0;
  }
}

void ConstraintList::add_forces(Particle *p) {
  double folded_pos[3];
  int img[3];

  memcpy(folded_pos, p->r.p, 3*sizeof(double));
  memcpy(img, p->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);
 
  for(iterator it = begin(); it != end(); ++it) {
    (*it)->add_force(p, folded_pos);
  }
}

void ConstraintList::add_energies(Particle *p) {
  double folded_pos[3];
  int img[3];

  memcpy(folded_pos, p->r.p, 3*sizeof(double));
  memcpy(img, p->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);
 
  for(iterator it = begin(); it != end(); ++it) {
    (*it)->add_energy(p, folded_pos, energy);
  }
}

/** @TODO: Make less ugly. */

double ConstraintList::min_dist(double pos[3]) {
  double dist, vec[3];
  double mind = std::numeric_limits<double>::max();
  for(iterator it = begin(); it != end(); ++it) {
    const Constraints::GeometryConstraint *c = dynamic_cast<Constraints::GeometryConstraint *>(*it);
    if(c) {
      c->calculate_dist(pos, &dist, vec);
      mind = std::min(mind, dist);
    }
  }
  return mind;
}

}

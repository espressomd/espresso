#include "ConstraintList.hpp"
#include "GeometryConstraint.hpp"

#include "energy.hpp"
#include "errorhandling.hpp"

#include <limits>

ConstraintList constraintList;

void ConstraintList::init_forces() {
  for(iterator it = begin(); it != end(); ++it) {
    for(int i = 0; i < 3; i++)
      it->second->total_force[i] = 0.0;
  }
}

int ConstraintList::add_constraint(Constraints::Constraint *c) {
  /* Add c to local list */
  insert(value_type(m_next_id, c));

  /* c is now the last element */
  c->id = m_next_id++;

  /* Check if c is a GeometryConstraint in which case we have
     to do a deep copy of the Shape object. */
  if(dynamic_cast<Constraints::GeometryConstraint *>(c)) {
    ;
  } else {
    ;
  }

  return c->id;
}

void ConstraintList::remove_constraint(int i) {
  /* Remove entry if it exists */
  if(find(i) != end())
    erase(i);
}

void ConstraintList::add_forces(Particle *p) {
  double folded_pos[3];
  int img[3];

  memcpy(folded_pos, p->r.p, 3*sizeof(double));
  memcpy(img, p->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);
 
  for(iterator it = begin(); it != end(); ++it) {
    it->second->add_force(p, folded_pos);
  }
}

void ConstraintList::add_energies(Particle *p) {
  double folded_pos[3];
  int img[3];

  memcpy(folded_pos, p->r.p, 3*sizeof(double));
  memcpy(img, p->l.i, 3*sizeof(int));
  fold_position(folded_pos, img);
 
  for(iterator it = begin(); it != end(); ++it) {
    it->second->add_energy(p, folded_pos, energy);
  }
}

double ConstraintList::min_dist(double pos[3]) {
  double dist, vec[3];
  double mind = std::numeric_limits<double>::max();
  for(iterator it = begin(); it != end(); ++it) {
    Constraints::GeometryConstraint *c = dynamic_cast<Constraints::GeometryConstraint *>(it->second);
    if(c){
      c->m_shape->calculate_dist(pos, &dist, vec);
      mind = std::min(mind, dist);
    }
  }
  return mind;
}

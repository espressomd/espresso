#ifndef __CHARGED_ROD_HPP
#define __CHARGED_ROD_HPP

#include "Constraint.hpp"

namespace Constraints {
  struct ChargedRod : public Constraint {
    ChargedRod(ConstraintType _type, double center_x, double center_y, double _lambda) : Constraint(_type), lambda(_lambda) {
      center[0] = center_x;
      center[1] = center_y;
    }
    void add_energy(const Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    double center[2];
    double lambda;
  };
}

#endif

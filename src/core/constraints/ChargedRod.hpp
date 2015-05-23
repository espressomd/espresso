#ifndef __CHARGED_ROD_HPP
#define __CHARGED_ROD_HPP

#include "Constraint.hpp"

namespace Constraints {
  struct ChargedRod : public Constraint {
    ChargedRod(double center_x, double center_y, double _lambda) : lambda(_lambda) {
      center[0] = center_x;
      center[1] = center_y;
    }
    void add_energy(const Particle *p, const double *folded_pos, Observable_stat &energy);
    void add_force(Particle *p, const double *folded_pos);
    virtual std::string name() { return Constraint::name() + std::string("ChargedRod"); }
    ConstraintType type() { return CONSTRAINT_CHARGED_ROD; }
    double center[2];
    double lambda;
  };
}

#endif

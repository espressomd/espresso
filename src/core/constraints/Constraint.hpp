#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP

#include <string>

#include "energy.hpp"
#include "interaction_data.hpp"

namespace Constraints {

  enum ConstraintType { CONSTRAINT_NONE, CONSTRAINT_GEOMETRY, CONSTRAINT_INTERACTION };

  struct Constraint {
  public:
    Constraint() : id(-1) { };
    virtual ConstraintType type() { return CONSTRAINT_NONE; }
    virtual void add_energy(const Particle *p, const double *folded_pos, Observable_stat &energy) { }
    virtual void add_force(Particle *p, const double *folded_pos) { }
    /* Accumulated force excerted by the constraint */
    double total_force[3];
    /* Numerical id for interface binding */
    int id;
    /* Human readable name */
    virtual std::string name() { return std::string("Constraint::"); }
  };
};

#endif

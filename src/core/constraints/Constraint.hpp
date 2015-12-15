#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP

#include <string>

#include "ScriptObject.hpp"

#include "energy.hpp"
#include "interaction_data.hpp"
#include "utils/Factory.hpp"

namespace Constraints {
struct Constraint : public ScriptObject {
 public:
  virtual void add_energy(const Particle *p, const double *folded_pos, Observable_stat &energy) { }
  virtual void add_force(Particle *p, const double *folded_pos) { }
  /* Accumulated force excerted by the constraint */
  double total_force[3];
  /* Human readable name */
  virtual const std::string name() const { return std::string("Constraint::"); }
  virtual ~Constraint() {}
};

typedef Utils::Factory<Constraint> Factory;
void initialize_factory();

}

#endif

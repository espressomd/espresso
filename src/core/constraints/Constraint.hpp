#ifndef __CONSTRAINT_HPP
#define __CONSTRAINT_HPP

#include <string>

#include "script_interface/ParallelScriptObject.hpp"

#include "energy.hpp"
#include "interaction_data.hpp"
#include "utils/ParallelFactory.hpp"

namespace Constraints {
struct Constraint : public ScriptInterface::ParallelScriptObject {
 public:
  virtual void add_energy(Particle *p, const double *folded_pos, Observable_stat &energy) = 0;
  virtual void add_force(Particle *p, const double *folded_pos) = 0;
  /* Accumulated force excerted by the constraint */
  double total_force[3];
  /* Human readable name */
  virtual const std::string name() const { return std::string("Constraint::"); }
  virtual ~Constraint() {}
};

typedef Utils::ParallelFactory<Constraint> Factory;
void initialize_factory();

}

#endif

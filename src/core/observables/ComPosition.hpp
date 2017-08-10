#ifndef OBSERVABLES_COMPOSITION_HPP
#define OBSERVABLES_COMPOSITION_HPP

#include "PidObservable.hpp"
#include "partCfg.hpp"
#include <vector>

namespace Observables {

class ComPosition : public PidObservable {
public:
  virtual int n_values() const override { return 3; }
  virtual int actual_calculate() override {
    double total_mass = 0;
    for (int i = 0; i < ids().size(); i++) {
      double mass = partCfg[ids()[i]].p.mass;
      last_value[0] += mass * partCfg[ids()[i]].r.p[0];
      last_value[1] += mass * partCfg[ids()[i]].r.p[1];
      last_value[2] += mass * partCfg[ids()[i]].r.p[2];
      total_mass += mass;
    }
    last_value[0] /= total_mass;
    last_value[1] /= total_mass;
    last_value[2] /= total_mass;
    return 0;
  };
};

} // Namespace Observables
#endif

#ifndef OBSERVABLES_ComForce_HPP
#define OBSERVABLES_ComForce_HPP

#include "PidObservable.hpp"
#include "core/partCfg.hpp"
#include "integrate.hpp"

namespace Observables {

class ComForce : public PidObservable {
public:
  virtual int n_values() const override { return 3; }
  virtual int actual_calculate() override {

    double scale = 2 / time_step / time_step;
    for (int i = 0; i < ids().size(); i++) {
      last_value[0] += scale * partCfg[ids()[i]].f.f[0] * partCfg[ids()[i]].p.mass;
      last_value[1] += scale * partCfg[ids()[i]].f.f[1] * partCfg[ids()[i]].p.mass;
      last_value[2] += scale * partCfg[ids()[i]].f.f[2] * partCfg[ids()[i]].p.mass;
    }
    return 0;
  };
};

} // Namespace Observables
#endif

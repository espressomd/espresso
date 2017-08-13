#ifndef OBSERVABLES_CURRENTS_HPP
#define OBSERVABLES_CURRENTS_HPP

#include "PidObservable.hpp"
#include "partCfg.hpp"
#include <vector>

namespace Observables {

class Current : public PidObservable {
public:
  virtual int n_values() const override { return 3; };
  virtual int actual_calculate() override {
    last_value.resize(3);
    for (int i = 0; i < ids().size(); i++) {
#ifdef ELECTROSTATICS
      double charge = partCfg[ids()[i]].p.q;
      last_value[0] += charge * partCfg[ids()[i]].m.v[0] / time_step;
      last_value[1] += charge * partCfg[ids()[i]].m.v[1] / time_step;
      last_value[2] += charge * partCfg[ids()[i]].m.v[2] / time_step;
#endif
    };
    return 0;
  };
};

} // Namespace Observables
#endif

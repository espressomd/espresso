#ifndef OBSERVABLES_PARTICLECURRENTS_HPP
#define OBSERVABLES_PARTICLECURRENTS_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"

#include <vector>

namespace Observables {

class ParticleCurrent : public PidObservable {
public:
  virtual int actual_calculate(PartCfg & partCfg) override {
    last_value.resize(3 * ids().size());
    for (int i = 0; i < ids().size(); i++) {
#ifdef ELECTROSTATICS
      double charge = partCfg[ids()[i]].p.q;
      last_value[3 * i + 0] = charge * partCfg[ids()[i]].m.v[0] / time_step;
      last_value[3 * i + 1] = charge * partCfg[ids()[i]].m.v[1] / time_step;
      last_value[3 * i + 2] = charge * partCfg[ids()[i]].m.v[2] / time_step;
#endif
    };
    return 0;
  };
};

} // Namespace Observables
#endif

#ifndef OBSERVABLES_PARTICLEPOSITIONS_HPP
#define OBSERVABLES_PARTICLEPOSITIONS_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class ParticlePositions : public PidObservable {
public:
  virtual int actual_calculate(PartCfg & partCfg) override {
    for (int i = 0; i < ids().size(); i++) {
      last_value[3 * i + 0] = partCfg[ids()[i]].r.p[0];
      last_value[3 * i + 1] = partCfg[ids()[i]].r.p[1];
      last_value[3 * i + 2] = partCfg[ids()[i]].r.p[2];
    }
    return 0;
  }
};

} // Namespace Observables

#endif

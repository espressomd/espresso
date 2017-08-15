#ifndef OBSERVABLES_PARTICLEBODYANGULARMOMENTUM_HPP
#define OBSERVABLES_PARTICLEBODYANGULARMOMENTUM_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"
#include "partCfg.hpp"
#include <vector>

namespace Observables {

class ParticleBodyAngularMomentum : public PidObservable {
public:
  virtual int actual_calculate() override {
    last_value.resize(3 * ids().size());
    for (int i = 0; i < ids().size(); i++) {
#ifdef ROTATION

      last_value[3 * i + 0] = partCfg[ids()[i]].m.omega[0];
      last_value[3 * i + 1] = partCfg[ids()[i]].m.omega[1];
      last_value[3 * i + 2] = partCfg[ids()[i]].m.omega[2];
#endif
    }
    return 0;
  }
};

} // Namespace Observables
#endif

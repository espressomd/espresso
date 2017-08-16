#ifndef OBSERVABLES_PARTICLEFORCES_HPP
#define OBSERVABLES_PARTICLEFORCES_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include <vector>

namespace Observables {

class ParticleForces : public PidObservable {
public:
  virtual int actual_calculate() override {
    last_value.resize(3 * ids().size());
    double scale = 2. / time_step / time_step;
    for (int i = 0; i < ids().size(); i++) {
      last_value[3 * i + 0] =
          scale * partCfg[ids()[i]].f.f[0] * partCfg[ids()[i]].p.mass;
      last_value[3 * i + 1] =
          scale * partCfg[ids()[i]].f.f[1] * partCfg[ids()[i]].p.mass;
      last_value[3 * i + 2] =
          scale * partCfg[ids()[i]].f.f[2] * partCfg[ids()[i]].p.mass;
    }
    return 0;
  };
};

} // Namespace Observables
#endif

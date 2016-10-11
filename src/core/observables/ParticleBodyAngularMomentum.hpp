#ifndef OBSERVABLES_PARTICLEBODYANGULARMOMENTUM_HPP
#define OBSERVABLES_PARTICLEBODYANGULARMOMENTUM_HPP


#include "PidObservable.hpp"
#include "particle_data.hpp" 
#include <vector>
#include "integrate.hpp"  


namespace Observables {


class ParticleBodyAngularMomentum : public PidObservable {
public:
    virtual int actual_calculate() override {
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  last_value.resize(3*ids.size());
  for (int i = 0; i<ids.size(); i++ ) {
    if (ids[i] >= n_part)
      return 1;
#ifdef ROTATION

    last_value[3*i + 0] = partCfg[ids[i]].m.omega[0];
    last_value[3*i + 1] = partCfg[ids[i]].m.omega[1];
    last_value[3*i + 2] = partCfg[ids[i]].m.omega[2];
#endif


  }
  return 0;
}
};

} // Namespace Observables
#endif


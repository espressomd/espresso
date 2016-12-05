#ifndef OBSERVABLES_ComForce_HPP
#define OBSERVABLES_ComForce_HPP


#include "PidObservable.hpp"
#include "particle_data.hpp" 
#include <vector>


namespace Observables {


class ComForce : public PidObservable {
public:
    virtual int n_values() const override { return 3; }
    virtual int actual_calculate() override {
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  double scale=2/time_step/time_step;
  for (int i = 0; i<ids.size(); i++ ) {
    if (ids[i] >= n_part)
      return 1;
    last_value[0] += scale * partCfg[ids[i]].f.f[0] *partCfg[ids[i]].p.mass;
    last_value[1] += scale * partCfg[ids[i]].f.f[1] *partCfg[ids[i]].p.mass;
    last_value[2] += scale * partCfg[ids[i]].f.f[2] *partCfg[ids[i]].p.mass;
  }
  return 0;
};
};

} // Namespace Observables
#endif


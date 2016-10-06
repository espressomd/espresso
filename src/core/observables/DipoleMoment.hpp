#ifndef OBSERVABLES_DIPOLEMOMENT_HPP
#define OBSERVABLES_DIPOLEMOMENT_HPP

#include "PidObservable.hpp"
#include "particle_data.hpp" 
#include <vector>
#include "integrate.hpp"  


namespace Observables {


class DipoleMoment : public PidObservable {
public:
    virtual int n_values() const override { return 3; };
    virtual int actual_calculate() override {
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  for (int i = 0; i<ids.size(); i++ ) {
    if (ids[i] >= n_part)
      return 1;
    double charge = partCfg[ids[i]].p.q;
   
    last_value[0] += charge * partCfg[ids[i]].r.p[0];
    last_value[1] += charge * partCfg[ids[i]].r.p[1];
    last_value[2] += charge * partCfg[ids[i]].r.p[2];
  }
  return 0;
}
};

} // Namespace Observables

#endif


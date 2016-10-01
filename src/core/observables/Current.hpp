#ifndef OBSERVABLES_CURRENTS_HPP
#define OBSERVABLES_CURRENTS_HPP


#include "PidObservable.hpp"
#include "particle_data.hpp" 
#include <vector>
#include "integrate.hpp"  


namespace Observables {


class ParticleCurrent : public PidObservable {
public:
    virtual int n_values() override { return 3; }; 
    virtual int actual_calculate() override {
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  int count=0;
  last_value.resize(3);
  for (int i = 0; i<ids.size(); i++ ) {
    if (ids[i] >= n_part)
      return 1;
#ifdef ELECTROSTATICS
    charge = partCfg[ids[i]].p.q;
    last_value[3*i + 0] += charge * partCfg[ids[i]].m.v[0]/time_step;
    last_value[3*i + 1] += charge * partCfg[ids[i]].m.v[1]/time_step;
    last_value[3*i + 2] += charge * partCfg[ids[i]].m.v[2]/time_step;
    count++;
 #endif'
};
  if (count==0) {
    runtimeErrorMSg() << "No particles for Current observable.";
 
  last_value/=count;
  return 0;
};
};

} // Namespace Observables
#endif


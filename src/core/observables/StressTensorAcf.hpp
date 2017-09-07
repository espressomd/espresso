#ifndef OBSERVABLES_STRESSTENSORACF_HPP
#define OBSERVABLES_STRESSTENSORACF_HPP

#include "PidObservable.hpp"
#include "particle_data.hpp" 
#include <vector>
#include "pressure.hpp"


namespace Observables {


class StressTensorAcf : public PidObservable {
public:
    virtual int n_values() const override { return 6;};
    virtual int actual_calculate(PartCfg & partCfg) override {
  double stress_tensor[9];
  observable_compute_stress_tensor(1,stress_tensor);
  last_value[0]=stress_tensor[1];
  last_value[1]=stress_tensor[5];
  last_value[2]=stress_tensor[6];
  last_value[3]=stress_tensor[0]-stress_tensor[4];
  last_value[4]=stress_tensor[0]-stress_tensor[8];
  last_value[5]=stress_tensor[4]-stress_tensor[8];
  return 0;
}
};

} // Namespace Observables

#endif


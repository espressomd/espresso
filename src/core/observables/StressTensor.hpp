#ifndef OBSERVABLES_STRESSTENSOR_HPP

#include "PidObservable.hpp"
#include "particle_data.hpp" 
#include <vector>
#include "pressure.hpp"


namespace Observables {


class StressTensor : public PidObservable {
public:
    virtual int n_values() const override { return 9;};
    virtual int actual_calculate(PartCfg & partCfg) override {
  observable_compute_stress_tensor(1,&(last_value[0]));
  return 0;
}
};

} // Namespace Observables

#endif


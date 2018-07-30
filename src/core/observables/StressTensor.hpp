#ifndef OBSERVABLES_STRESSTENSOR_HPP

#include "Observable.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include <vector>

namespace Observables {

class StressTensor : public Observable {
public:
  virtual int n_values() const override { return 9; };
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    observable_compute_stress_tensor(1, res.data());
    return res;
  }
};

} // Namespace Observables

#endif

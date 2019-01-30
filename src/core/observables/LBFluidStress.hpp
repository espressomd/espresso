#ifndef OBSERVABLES_LB_FLUID_STRESS_HPP
#define OBSERVABLES_LB_FLUID_STRESS_HPP

#include "Observable.hpp"
#include "grid_based_algorithms/lb.hpp"

#include <vector>

namespace Observables {

class LBFluidStress : public Observable {
public:
  int n_values() const override {
    return 6;
  }
  std::vector<double> operator()(PartCfg &) const override {
    std::vector<double> res(n_values(), 0.0);

#if defined(LB_GPU) || defined(LB)
    lb_lbfluid_get_pi(res.data());
#endif

    return res;
  }
};

} // Namespace Observables

#endif

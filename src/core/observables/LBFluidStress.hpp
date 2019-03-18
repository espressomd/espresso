#ifndef OBSERVABLES_LB_FLUID_STRESS_HPP
#define OBSERVABLES_LB_FLUID_STRESS_HPP

#include "Observable.hpp"
#include "grid_based_algorithms/lb_interface.hpp"

#include <vector>

namespace Observables {
#if (defined(LB) || defined(LB_GPU))
class LBFluidStress : public Observable {
public:
  int n_values() const override { return 6; }
  std::vector<double> operator()(PartCfg &) const override {

    return lb_lbfluid_get_stress();
  }
};
#endif

} // Namespace Observables

#endif

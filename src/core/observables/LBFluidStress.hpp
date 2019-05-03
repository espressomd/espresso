#ifndef OBSERVABLES_LB_FLUID_STRESS_HPP
#define OBSERVABLES_LB_FLUID_STRESS_HPP

#include "Observable.hpp"
#include "grid_based_algorithms/lb_interface.hpp"

#include <vector>

namespace Observables {
class LBFluidStress : public Observable {
public:
  int n_values() const override { return 6; }
  std::vector<double> operator()(PartCfg &) const override {

    auto const unit_conversion =
        1. / (lb_lbfluid_get_agrid() * pow(lb_lbfluid_get_tau(), 2));
    return lb_lbfluid_get_stress() * unit_conversion;
  }
};

} // Namespace Observables

#endif

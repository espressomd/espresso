#ifndef OBSERVABLES_LB_BULK_STRESS_HPP
#define OBSERVABLES_LB_BULK_STRESS_HPP

#include "Observable.hpp"
#include "lb.hpp"

#include <vector>

namespace Observables {

class LBBulkStress : public Observable {
public:
  virtual int n_values() const override {
    return 6;
  }
  virtual std::vector<double> operator()(PartCfg &) const override {
    std::vector<double> res(n_values(), 0.0);

#if defined(LB_GPU) || defined(LB)
    lb_bulk_get_pi(res.data());
#endif

    return res;
  }
};

} // Namespace Observables

#endif

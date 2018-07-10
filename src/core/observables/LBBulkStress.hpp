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

    lb_bulk_get_pi(res.data());

    return res;
  }
};

} // Namespace Observables

#endif

#ifndef OBSERVABLES_DPDSTRESS_HPP
#define OBSERVABLES_DPDSTRESS_HPP

#include "Observable.hpp"
#include "particle_data.hpp"
#include "pressure.hpp"
#include <vector>

namespace Observables {

class DPDStress : public Observable {
public:
  int n_values() const override { return 9; };
  std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    dpd_stress();
    return res;
  }
};

} // Namespace Observables

#endif

#ifndef OBSERVABLES_COMPOSITION_HPP
#define OBSERVABLES_COMPOSITION_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class ComPosition : public PidObservable {
public:
  int n_values() const override { return 3; }
  std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    double total_mass = 0;
    for (int i = 0; i < ids().size(); i++) {
      double mass = partCfg[ids()[i]].p.mass;
      res[0] += mass * partCfg[ids()[i]].r.p[0];
      res[1] += mass * partCfg[ids()[i]].r.p[1];
      res[2] += mass * partCfg[ids()[i]].r.p[2];
      total_mass += mass;
    }
    res[0] /= total_mass;
    res[1] /= total_mass;
    res[2] /= total_mass;
    return res;
  };
};

} // Namespace Observables
#endif

#ifndef OBSERVABLES_COMVELOCITY_HPP
#define OBSERVABLES_COMVELOCITY_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class ComVelocity : public PidObservable {
public:
  virtual int n_values() const override { return 3; }
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    double total_mass = 0;
    for (int i = 0; i < ids().size(); i++) {
      double mass = partCfg[ids()[i]].p.mass;
      res[0] += mass * partCfg[ids()[i]].m.v[0];
      res[1] += mass * partCfg[ids()[i]].m.v[1];
      res[2] += mass * partCfg[ids()[i]].m.v[2];
      total_mass += mass;
    }
    res[0] /= total_mass * time_step;
    res[1] /= total_mass * time_step;
    res[2] /= total_mass * time_step;
    return res;
  };
};

} // Namespace Observables
#endif

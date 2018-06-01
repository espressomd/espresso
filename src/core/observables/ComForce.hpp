#ifndef OBSERVABLES_ComForce_HPP
#define OBSERVABLES_ComForce_HPP

#include "PidObservable.hpp"

#include "integrate.hpp"

namespace Observables {

class ComForce : public PidObservable {
public:
  virtual int n_values() const override { return 3; }
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    double scale = 2 / time_step / time_step;
    for (int i = 0; i < ids().size(); i++) {
      res[0] += scale * partCfg[ids()[i]].f.f[0] * partCfg[ids()[i]].p.mass;
      res[1] += scale * partCfg[ids()[i]].f.f[1] * partCfg[ids()[i]].p.mass;
      res[2] += scale * partCfg[ids()[i]].f.f[2] * partCfg[ids()[i]].p.mass;
    }
    return res;
  };
};

} // Namespace Observables
#endif

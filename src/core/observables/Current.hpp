#ifndef OBSERVABLES_CURRENTS_HPP
#define OBSERVABLES_CURRENTS_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class Current : public PidObservable {
public:
  virtual int n_values() const override { return 3; };
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
#ifdef ELECTROSTATICS
      double charge = partCfg[ids()[i]].p.q;
      res[0] += charge * partCfg[ids()[i]].m.v[0] / time_step;
      res[1] += charge * partCfg[ids()[i]].m.v[1] / time_step;
      res[2] += charge * partCfg[ids()[i]].m.v[2] / time_step;
#endif
    };
    return res;
  };
};

} // Namespace Observables
#endif

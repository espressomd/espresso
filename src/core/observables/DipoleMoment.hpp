#ifndef OBSERVABLES_DIPOLEMOMENT_HPP
#define OBSERVABLES_DIPOLEMOMENT_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class DipoleMoment : public PidObservable {
public:
  virtual int n_values() const override { return 3; };
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values(), 0.0);
    for (int i = 0; i < ids().size(); i++) {
#ifdef ELECTROSTATICS
      double charge = partCfg[ids()[i]].p.q;

      res[0] += charge * partCfg[ids()[i]].r.p[0];
      res[1] += charge * partCfg[ids()[i]].r.p[1];
      res[2] += charge * partCfg[ids()[i]].r.p[2];
#endif // ELECTROSTATICS
    }
    return res;
  }
};

} // Namespace Observables

#endif

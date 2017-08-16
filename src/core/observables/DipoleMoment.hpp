#ifndef OBSERVABLES_DIPOLEMOMENT_HPP
#define OBSERVABLES_DIPOLEMOMENT_HPP

#include "PidObservable.hpp"
#include "partCfg.hpp"
#include <vector>

namespace Observables {

class DipoleMoment : public PidObservable {
public:
  virtual int n_values() const override { return 3; };
  virtual int actual_calculate() override {
    for (int i = 0; i < ids().size(); i++) {
#ifdef ELECTROSTATICS
      double charge = partCfg[ids()[i]].p.q;

      last_value[0] += charge * partCfg[ids()[i]].r.p[0];
      last_value[1] += charge * partCfg[ids()[i]].r.p[1];
      last_value[2] += charge * partCfg[ids()[i]].r.p[2];
#endif // ELECTROSTATICS
    }
    return 0;
  }
};

} // Namespace Observables

#endif

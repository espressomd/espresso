#ifndef OBSERVABLES_MAGNETICDIPOLEMOMENT_HPP
#define OBSERVABLES_MAGNETICDIPOLEMOMENT_HPP

#include "PidObservable.hpp"
#include "partCfg.hpp"
#include <vector>

namespace Observables {

class MagneticDipoleMoment : public PidObservable {
public:
  virtual int n_values() const override { return 3; };
  virtual int actual_calculate() override {
    for (int i = 0; i < ids().size(); i++) {
#ifdef DIPOLES
      last_value[0] += partCfg[ids()[i]].r.dip[0];
      last_value[1] += partCfg[ids()[i]].r.dip[1];
      last_value[2] += partCfg[ids()[i]].r.dip[2];
#endif
    }
    return 0;
  }
};

} // Namespace Observables

#endif

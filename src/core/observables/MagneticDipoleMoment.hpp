#ifndef OBSERVABLES_MAGNETICDIPOLEMOMENT_HPP
#define OBSERVABLES_MAGNETICDIPOLEMOMENT_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class MagneticDipoleMoment : public PidObservable {
public:
  virtual int n_values() const override { return 3; };
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values(), 0.0);
    for (int i = 0; i < ids().size(); i++) {
#ifdef DIPOLES
      res[0] += partCfg[ids()[i]].r.dip[0];
      res[1] += partCfg[ids()[i]].r.dip[1];
      res[2] += partCfg[ids()[i]].r.dip[2];
#endif
    }
    return res;
  }
};

} // Namespace Observables

#endif

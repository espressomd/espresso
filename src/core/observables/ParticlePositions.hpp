#ifndef OBSERVABLES_PARTICLEPOSITIONS_HPP
#define OBSERVABLES_PARTICLEPOSITIONS_HPP

#include "PidObservable.hpp"

#include <vector>

namespace Observables {

class ParticlePositions : public PidObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
      res[3 * i + 0] = partCfg[ids()[i]].r.p[0];
      res[3 * i + 1] = partCfg[ids()[i]].r.p[1];
      res[3 * i + 2] = partCfg[ids()[i]].r.p[2];
    }
    return res;
  }
  virtual int n_values() const override { return 3 * ids().size(); }
};

} // Namespace Observables

#endif

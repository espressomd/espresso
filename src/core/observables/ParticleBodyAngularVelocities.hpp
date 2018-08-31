#ifndef OBSERVABLES_PARTICLEBODYANGULARVELOCITIES_HPP
#define OBSERVABLES_PARTICLEBODYANGULARVELOCITIES_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"

#include <vector>

namespace Observables {

class ParticleBodyAngularVelocities : public PidObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
#ifdef ROTATION
      res[3 * i + 0] = partCfg[ids()[i]].m.omega[0];
      res[3 * i + 1] = partCfg[ids()[i]].m.omega[1];
      res[3 * i + 2] = partCfg[ids()[i]].m.omega[2];
#endif
    }
    return res;
  }
  virtual int n_values() const override { return 3 * ids().size(); }
};

} // Namespace Observables
#endif

#ifndef OBSERVABLES_PARTICLEFORCES_HPP
#define OBSERVABLES_PARTICLEFORCES_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"
#include "particle_data.hpp"
#include <vector>

namespace Observables {

class ParticleForces : public PidObservable {
public:
  std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
      res[3 * i + 0] = partCfg[ids()[i]].f.f[0];
      res[3 * i + 1] = partCfg[ids()[i]].f.f[1];
      res[3 * i + 2] = partCfg[ids()[i]].f.f[2];
    }
    return res;
  };
  int n_values() const override { return 3 * ids().size(); }
};

} // Namespace Observables
#endif

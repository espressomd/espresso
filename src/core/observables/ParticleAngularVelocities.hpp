#ifndef OBSERVABLES_PARTICLEANGULARVELOCITIES_HPP
#define OBSERVABLES_PARTICLEANGULARVELOCITIES_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"

#include "rotation.hpp"
#include <vector>

namespace Observables {

class ParticleAngularVelocities : public PidObservable {
public:
  std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    for (int i = 0; i < ids().size(); i++) {
#ifdef ROTATION

      double RMat[9];
      double omega[3];
      define_rotation_matrix(partCfg[ids()[i]], RMat);
      omega[0] = RMat[0 + 3 * 0] * partCfg[ids()[i]].m.omega[0] +
                 RMat[1 + 3 * 0] * partCfg[ids()[i]].m.omega[1] +
                 RMat[2 + 3 * 0] * partCfg[ids()[i]].m.omega[2];
      omega[1] = RMat[0 + 3 * 1] * partCfg[ids()[i]].m.omega[0] +
                 RMat[1 + 3 * 1] * partCfg[ids()[i]].m.omega[1] +
                 RMat[2 + 3 * 1] * partCfg[ids()[i]].m.omega[2];
      omega[2] = RMat[0 + 3 * 2] * partCfg[ids()[i]].m.omega[0] +
                 RMat[1 + 3 * 2] * partCfg[ids()[i]].m.omega[1] +
                 RMat[2 + 3 * 2] * partCfg[ids()[i]].m.omega[2];

      res[3 * i + 0] = omega[0];
      res[3 * i + 1] = omega[1];
      res[3 * i + 2] = omega[2];
#endif
    }
    return res;
  }
  int n_values() const override { return 3 * ids().size(); }
};

} // Namespace Observables
#endif

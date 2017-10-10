#ifndef OBSERVABLES_PARTICLEBODYVELOCITIES_HPP
#define OBSERVABLES_PARTICLEBODYVELOCITIES_HPP

#include "PidObservable.hpp"
#include "integrate.hpp"

#include "rotation.hpp"
#include <vector>

namespace Observables {

class ParticleBodyVelocities : public PidObservable {
public:
  virtual int actual_calculate(PartCfg & partCfg) override {
    last_value.resize(3 * ids().size());
    for (int i = 0; i < ids().size(); i++) {
#ifdef ROTATION

      double RMat[9];
      double vel_lab[3];
      double vel_body[3];

      vel_lab[0] = partCfg[ids()[i]].m.v[0] / time_step;
      vel_lab[1] = partCfg[ids()[i]].m.v[1] / time_step;
      vel_lab[2] = partCfg[ids()[i]].m.v[2] / time_step;
      define_rotation_matrix(&partCfg[ids()[i]], RMat);

      vel_body[0] = RMat[0 + 3 * 0] * vel_lab[0] +
                    RMat[0 + 3 * 1] * vel_lab[1] + RMat[0 + 3 * 2] * vel_lab[2];
      vel_body[1] = RMat[1 + 3 * 0] * vel_lab[0] +
                    RMat[1 + 3 * 1] * vel_lab[1] + RMat[1 + 3 * 2] * vel_lab[2];
      vel_body[2] = RMat[2 + 3 * 0] * vel_lab[0] +
                    RMat[2 + 3 * 1] * vel_lab[1] + RMat[2 + 3 * 2] * vel_lab[2];

      last_value[3 * i + 0] = vel_body[0];
      last_value[3 * i + 1] = vel_body[1];
      last_value[3 * i + 2] = vel_body[2];

#endif
    }
    return 0;
  }
};

} // Namespace Observables
#endif

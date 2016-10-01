#ifndef OBSERVABLES_PARTICLEANGULARMOMENTUM_HPP
#define OBSERVABLES_PARTICLEANGULARMOMENTUM_HPP


#include "PidObservable.hpp"
#include "particle_data.hpp" 
#include <vector>
#include "integrate.hpp"  


namespace Observables {


class ParticleAngularMomentum : public PidObservable {
public:
    virtual int actual_calculate() override {
  if (!sortPartCfg()) {
      runtimeErrorMsg() <<"could not sort partCfg";
    return -1;
  }
  last_value.resize(3*ids.size());
  for (int i = 0; i<ids.size(); i++ ) {
    if (ids[i] >= n_part)
      return 1;
#ifdef ROTATION

    double RMat[9];
    double omega[3];
    define_rotation_matrix(&partCfg[ids->e[i]], RMat);
    omega[0] = RMat[0 + 3*0]*partCfg[ids[i]].m.omega[0] + RMat[1 + 3*0]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*0]*partCfg[ids->e[i]].m.omega[2];
    omega[1] = RMat[0 + 3*1]*partCfg[ids[i]].m.omega[0] + RMat[1 + 3*1]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*1]*partCfg[ids->e[i]].m.omega[2];
    omega[2] = RMat[0 + 3*2]*partCfg[ids[i]].m.omega[0] + RMat[1 + 3*2]*partCfg[ids->e[i]].m.omega[1] + RMat[2 + 3*2]*partCfg[ids->e[i]].m.omega[2];

    last_value[3*i + 0] = omega[0];
    last_value[3*i + 1] = omega[1];
    last_value[3*i + 2] = omega[2];
#endif

  }
  return 0;
}
};

} // Namespace Observables
#endif


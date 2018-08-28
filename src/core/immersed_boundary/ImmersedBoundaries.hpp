#ifndef IMMERSED_BOUNDARY_IMMERSED_BOUNDARIES_HPP
#define IMMERSED_BOUNDARY_IMMERSED_BOUNDARIES_HPP

#include "config.hpp"
#include <vector>

#ifdef IMMERSED_BOUNDARY
class ImmersedBoundaries {
public:
  ImmersedBoundaries() : MaxNumIBM(1000), VolumeInitDone(false) {
    VolumesCurrent.resize(MaxNumIBM);
  }
  void init_volume_conservation();
  void volume_conservation();
  int volume_conservation_reset_params(const int bond_type,
                                       const double volRef);
  int volume_conservation_set_params(const int bond_type, const int softID,
                                     const double kappaV);
  void calc_volumes();
  void calc_volume_force();

private:
  const int MaxNumIBM;
  std::vector<double> VolumesCurrent;
  bool VolumeInitDone = false;
};

#endif
#endif

#ifndef OBSERVABLES_DENSITYPROFILE_HPP
#define OBSERVABLES_DENSITYPROFILE_HPP

#include "ProfileObservable.hpp"

#include <vector>

namespace Observables {

class DensityProfile : public ProfileObservable {
public:
  virtual int actual_calculate(PartCfg & partCfg) override {
    double bin_volume =
        (maxx - minx) * (maxy - miny) * (maxz - minz) / xbins / ybins / zbins;

    for (int id : ids) {
      auto const ppos = folded_position(partCfg[id]);

      int binx = (int)floor(xbins * (ppos[0] - minx) / (maxx - minx));
      int biny = (int)floor(ybins * (ppos[1] - miny) / (maxy - miny));
      int binz = (int)floor(zbins * (ppos[2] - minz) / (maxz - minz));
      if (binx >= 0 && binx < xbins && biny >= 0 && biny < ybins && binz >= 0 &&
          binz < zbins) {
        last_value[binx * ybins * zbins + biny * zbins + binz] +=
            1. / bin_volume;
      }
    }
    return 0;
  }
};

} // Namespace Observables

#endif

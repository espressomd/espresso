#ifndef OBSERVABLES_FLUXDENSITYPROFILE_HPP
#define OBSERVABLES_FLUXDENSITYPROFILE_HPP

#include "ProfileObservable.hpp"

#include <vector>

namespace Observables {
class FluxDensityProfile : public ProfileObservable {
public:
  virtual int n_values() const override { return 3 * xbins * ybins * zbins; }
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<double> res(n_values());
    double bin_volume =
        (maxx - minx) * (maxy - miny) * (maxz - minz) / xbins / ybins / zbins;

    for (int id : ids) {
      auto const ppos = folded_position(partCfg[id]);
      int binx = (int)floor(xbins * (ppos[0] - minx) / (maxx - minx));
      int biny = (int)floor(ybins * (ppos[1] - miny) / (maxy - miny));
      int binz = (int)floor(zbins * (ppos[2] - minz) / (maxz - minz));
      if (binx >= 0 && binx < xbins && biny >= 0 && biny < ybins && binz >= 0 &&
          binz < zbins) {
        for (int dim = 0; dim < 3; dim++)
          res[3 * (binx * ybins * zbins + biny * zbins + binz) + dim] +=
              partCfg[id].m.v[dim] / bin_volume;
      }
    }
    return res;
  }
};

} // Namespace Observables

#endif

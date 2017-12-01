#ifndef OBSERVABLES_FLUXDENSITYPROFILE_HPP
#define OBSERVABLES_FLUXDENSITYPROFILE_HPP

#include "ProfileObservable.hpp"

#include <vector>

namespace Observables {
class FluxDensityProfile : public ProfileObservable {
public:
  virtual int n_values() const override { return 3 * xbins * ybins * zbins; }
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<size_t> n_bins{{static_cast<size_t>(xbins),
                                static_cast<size_t>(ybins),
                                static_cast<size_t>(zbins)}};
    std::vector<std::pair<double, double>> limits{
        {std::make_pair(minx, maxx), std::make_pair(miny, maxy),
         std::make_pair(minz, maxz)}};
    Utils::Histogram<double> histogram(n_bins, 3, limits);
    for (int id : ids) {
      auto const ppos = ::Vector<3, double>(folded_position(partCfg[id]));
      histogram.update(ppos, ::Vector<3, double> {{partCfg[id].m.v[0],
                                                   partCfg[id].m.v[1],
                                                   partCfg[id].m.v[2]}});
      }
    histogram.normalize();
    return histogram.get_histogram();
  }
};

} // Namespace Observables

#endif

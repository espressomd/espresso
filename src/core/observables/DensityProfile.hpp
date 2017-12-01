#ifndef OBSERVABLES_DENSITYPROFILE_HPP
#define OBSERVABLES_DENSITYPROFILE_HPP

#include <vector>
#include "ProfileObservable.hpp"
#include "utils/Histogram.hpp"

namespace Observables {

class DensityProfile : public ProfileObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<size_t> n_bins{{static_cast<size_t>(xbins),
                                static_cast<size_t>(ybins),
                                static_cast<size_t>(zbins)}};
    std::vector<std::pair<double, double>> limits{
        {std::make_pair(minx, maxx), std::make_pair(miny, maxy),
         std::make_pair(minz, maxz)}};
    Utils::Histogram<double> histogram(n_bins, 1, limits);
    for (int id : ids) {
      auto const ppos = ::Vector<3, double>(folded_position(partCfg[id]));
      histogram.update(ppos);
    }
    histogram.normalize();
    return histogram.get_histogram();
  }
};

} // Namespace Observables

#endif

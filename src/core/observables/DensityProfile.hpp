#ifndef OBSERVABLES_DENSITYPROFILE_HPP
#define OBSERVABLES_DENSITYPROFILE_HPP

#include "ProfileObservable.hpp"
#include "utils/Histogram.hpp"
#include <vector>

namespace Observables {

class DensityProfile : public ProfileObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::array<size_t, 3> n_bins{{static_cast<size_t>(n_x_bins),
                                static_cast<size_t>(n_y_bins),
                                static_cast<size_t>(n_z_bins)}};
    std::array<std::pair<double, double>, 3> limits{
        {std::make_pair(min_x, max_x), std::make_pair(min_y, max_y),
         std::make_pair(min_z, max_z)}};
    Utils::Histogram<double, 3> histogram(n_bins, 1, limits);
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

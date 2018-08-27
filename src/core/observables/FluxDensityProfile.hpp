#ifndef OBSERVABLES_FLUXDENSITYPROFILE_HPP
#define OBSERVABLES_FLUXDENSITYPROFILE_HPP

#include "PidProfileObservable.hpp"

#include <vector>

namespace Observables {
class FluxDensityProfile : public PidProfileObservable {
public:
  virtual int n_values() const override {
    return 3 * n_x_bins * n_y_bins * n_z_bins;
  }
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::array<size_t, 3> n_bins{{static_cast<size_t>(n_x_bins),
                                  static_cast<size_t>(n_y_bins),
                                  static_cast<size_t>(n_z_bins)}};
    std::array<std::pair<double, double>, 3> limits{
        {std::make_pair(min_x, max_x), std::make_pair(min_y, max_y),
         std::make_pair(min_z, max_z)}};
    Utils::Histogram<double, 3> histogram(n_bins, 3, limits);
    for (auto const &id : ids()) {
      auto const ppos = ::Vector<3, double>(folded_position(partCfg[id]));
      histogram.update(
          ppos, ::Vector<3, double>{{partCfg[id].m.v[0], partCfg[id].m.v[1],
                                     partCfg[id].m.v[2]}});
    }
    histogram.normalize();
    return histogram.get_histogram();
  }
};

} // Namespace Observables

#endif

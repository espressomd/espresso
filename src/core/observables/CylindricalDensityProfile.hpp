#ifndef OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP

#include "CylindricalProfileObservable.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"
#include "utils/coordinate_transformation.hpp"

namespace Observables {
class CylindricalDensityProfile : public CylindricalProfileObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<size_t> n_bins{{static_cast<size_t>(n_r_bins),
                                static_cast<size_t>(n_phi_bins),
                                static_cast<size_t>(n_z_bins)}};
    std::vector<std::pair<double, double>> limits{
        {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
         std::make_pair(min_z, max_z)}};
    Utils::CylindricalHistogram<double> histogram(n_bins, 1, limits);
    auto folded_positions = Utils::get_folded_positions(ids());
    for (auto &p : folded_positions) {
      p -= center;
      histogram.update(Utils::transform_pos_to_cylinder_coordinates(p, axis));
    }
    histogram.normalize();
    return histogram.get_histogram();
  }
  virtual int n_values() const override {
    return n_r_bins * n_phi_bins * n_z_bins;
  }
};

} // Namespace Observables

#endif

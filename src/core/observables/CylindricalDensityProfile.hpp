#ifndef OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP

#include "CylindricalPidProfileObservable.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"
#include "utils/coordinate_transformation.hpp"

namespace Observables {
class CylindricalDensityProfile : public CylindricalPidProfileObservable {
public:
  std::vector<double> operator()(PartCfg &partCfg) const override {
    std::array<size_t, 3> n_bins{{static_cast<size_t>(n_r_bins),
                                  static_cast<size_t>(n_phi_bins),
                                  static_cast<size_t>(n_z_bins)}};
    std::array<std::pair<double, double>, 3> limits{
        {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
         std::make_pair(min_z, max_z)}};
    Utils::CylindricalHistogram<double, 3> histogram(n_bins, 1, limits);
    std::vector<::Vector<3, double>> folded_positions;
    std::transform(ids().begin(), ids().end(),
                   std::back_inserter(folded_positions), [&partCfg](int id) {
                     return ::Vector<3, double>(folded_position(partCfg[id]));
                   });
    for (auto &p : folded_positions) {
      p -= center;
      histogram.update(Utils::transform_pos_to_cylinder_coordinates(p, axis));
    }
    histogram.normalize();
    return histogram.get_histogram();
  }
  int n_values() const override {
    return n_r_bins * n_phi_bins * n_z_bins;
  }
};

} // Namespace Observables

#endif

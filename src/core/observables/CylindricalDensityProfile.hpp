#ifndef OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALDENSITYPROFILE_HPP

#include "CylindricalProfileObservable.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"

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
    for (int id : ids()) {
      auto const ppos = ::Vector<3, double>(folded_position(partCfg[id]));
      ::Vector<3, double> ppos_shifted;
      ppos_shifted = ppos - center;
      if (axis == "x") {
        // x' = -z, y' = y, z'= x
        ppos_shifted = ::Vector<3, double>{-ppos_shifted[2], ppos_shifted[1],
                                           ppos_shifted[0]};
      } else if (axis == "y") {
        // x' = x, y' = -z, z' = y
        ppos_shifted = ::Vector<3, double>{ppos_shifted[0], -ppos_shifted[2],
                                           ppos_shifted[1]};
      }

      auto const ppos_cyl =
          Utils::transform_to_cylinder_coordinates(ppos_shifted);
      histogram.update(ppos_cyl);
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

#ifndef OBSERVABLES_CYLINDRICALVELOCITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALVELOCITYPROFILE_HPP

#include "CylindricalPidProfileObservable.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"

namespace Observables {
class CylindricalVelocityProfile : public CylindricalPidProfileObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::array<size_t, 3> n_bins{{static_cast<size_t>(n_r_bins),
                                  static_cast<size_t>(n_phi_bins),
                                  static_cast<size_t>(n_z_bins)}};
    std::array<std::pair<double, double>, 3> limits{
        {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
         std::make_pair(min_z, max_z)}};
    Utils::CylindricalHistogram<double, 3> histogram(n_bins, 3, limits);
    std::vector<::Vector<3, double>> folded_positions;
    std::transform(ids().begin(), ids().end(),
                   std::back_inserter(folded_positions), [&partCfg](int id) {
                     return ::Vector<3, double>(folded_position(partCfg[id]));
                   });
    std::vector<::Vector<3, double>> velocities;
    std::transform(ids().begin(), ids().end(), std::back_inserter(velocities),
                   [&partCfg](int id) {
                     return ::Vector<3, double>{{partCfg[id].m.v[0],
                                                 partCfg[id].m.v[1],
                                                 partCfg[id].m.v[2]}};
                   });
    for (auto &p : folded_positions)
      p -= center;
    // Write data to the histogram
    for (size_t ind = 0; ind < folded_positions.size(); ++ind) {
      histogram.update(Utils::transform_pos_to_cylinder_coordinates(
                           folded_positions[ind], axis),
                       Utils::transform_vel_to_cylinder_coordinates(
                           velocities[ind], axis, folded_positions[ind]));
    }
    auto hist_tmp = histogram.get_histogram();
    auto tot_count = histogram.get_tot_count();
    for (size_t ind = 0; ind < hist_tmp.size(); ++ind) {
      if (tot_count[ind] > 0) {
        hist_tmp[ind] /= tot_count[ind];
      }
    }
    return hist_tmp;
  }
  virtual int n_values() const override {
    return 3 * n_r_bins * n_phi_bins * n_z_bins;
  }
};

} // Namespace Observables

#endif

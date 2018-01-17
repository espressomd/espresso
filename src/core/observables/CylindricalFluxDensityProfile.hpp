#ifndef OBSERVABLES_CYLINDRICALFLUXDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALFLUXDENSITYPROFILE_HPP

#include "CylindricalProfileObservable.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"

namespace Observables {
class CylindricalFluxDensityProfile : public CylindricalProfileObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<size_t> n_bins{{static_cast<size_t>(n_r_bins),
                                static_cast<size_t>(n_phi_bins),
                                static_cast<size_t>(n_z_bins)}};
    std::vector<std::pair<double, double>> limits{
        {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
         std::make_pair(min_z, max_z)}};
    Utils::CylindricalHistogram<double> histogram(n_bins, 3, limits);
    auto folded_positions = Utils::get_folded_positions(ids());
    auto velocities = Utils::get_velocities(ids());
    for (auto &p : folded_positions)
      p -= center;
    // Write data to the histogram
    for (size_t ind = 0; ind < ids().size(); ++ind) {
      histogram.update(Utils::transform_pos_to_cylinder_coordinates(
                           folded_positions[ind], axis),
                       Utils::transform_vel_to_cylinder_coordinates(
                           velocities[ind], axis, folded_positions[ind]) /
                           time_step);
    }
    histogram.normalize();
    return histogram.get_histogram();
  }
  virtual int n_values() const override {
    return 3 * n_r_bins * n_phi_bins * n_z_bins;
  }

private:
  virtual void do_write() override {
    // We override the implementation to actually write positions not plain
    // indices.
    static const int len_dims[4] = {n_r_bins, n_phi_bins, n_z_bins, 3};
    static const int n_dims = 4;
    static const std::array<double, 3> bin_sizes = {
        {r_bin_size(), phi_bin_size(), z_bin_size()}};
    std::array<double, 3> position;
    int index;
    int unravelled_index[4];
    std::vector<double> tmp = operator()(partCfg());
    for (auto it = tmp.begin(); it != tmp.end(); it += 3) {
      index = std::distance(tmp.begin(), it);
      ::Utils::unravel_index(len_dims, n_dims, index, unravelled_index);
      position = {
          {(static_cast<double>(unravelled_index[0]) + 0.5) * bin_sizes[0],
           (static_cast<double>(unravelled_index[1]) + 0.5) * bin_sizes[1],
           (static_cast<double>(unravelled_index[2]) + 0.5) * bin_sizes[2]}};
      m_ofile << position[0] << " " << position[1] << " " << position[2] << " "
              << *it << " " << *(it + 1) << " " << *(it + 2) << "\n";
    }
    m_ofile << std::endl;
  }
};

} // Namespace Observables

#endif

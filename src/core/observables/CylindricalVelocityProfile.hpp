#ifndef OBSERVABLES_CYLINDRICALVELOCITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALVELOCITYPROFILE_HPP

#include "CylindricalProfileObservable.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"

namespace Observables {
class CylindricalVelocityProfile : public CylindricalProfileObservable {
public:
  virtual std::vector<double> operator()(PartCfg &partCfg) const override {
    std::vector<size_t> n_bins{{static_cast<size_t>(n_r_bins),
                                static_cast<size_t>(n_phi_bins),
                                static_cast<size_t>(n_z_bins)}};
    std::vector<std::pair<double, double>> limits{
        {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
         std::make_pair(min_z, max_z)}};
    Utils::CylindricalHistogram<double> histogram(n_bins, 3, limits);
    for (int id : ids()) {
      auto const ppos = ::Vector<3, double>(folded_position(partCfg[id]));
      ::Vector<3, double> ppos_shifted;
      ppos_shifted = ppos - center;
      ::Vector<3, double> vel = {partCfg[id].m.v[0], partCfg[id].m.v[1],
                                 partCfg[id].m.v[2]};
      if (axis == "x") {
        // x' = -z, y' = y, z'= x
        ppos_shifted = ::Vector<3, double>{-ppos_shifted[2], ppos_shifted[1],
                                           ppos_shifted[0]};
        vel = {-vel[2], vel[1], vel[0]};
      } else if (axis == "y") {
        // x' = x, y' = -z, z' = y
        ppos_shifted = ::Vector<3, double>{ppos_shifted[0], -ppos_shifted[2],
                                           ppos_shifted[1]};
        vel = {vel[0], -vel[2], vel[1]};
      }

      auto const ppos_cyl =
          Utils::transform_to_cylinder_coordinates(ppos_shifted);

      // Coordinate transform the velocities and divide core velocities by
      // time_step to get MD units. v_r = (x * v_x + y * v_y) / sqrt(x^2 +
      // y^2)
      double v_r = (ppos_shifted[0] * vel[0] / time_step +
                    ppos_shifted[1] * vel[1] / time_step) /
                   std::sqrt(ppos_shifted[0] * ppos_shifted[0] +
                             ppos_shifted[1] * ppos_shifted[1]);
      // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
      double v_phi = (ppos_shifted[0] * vel[1] / time_step -
                      ppos_shifted[1] * vel[0] / time_step) /
                     (ppos_shifted[0] * ppos_shifted[0] +
                      ppos_shifted[1] * ppos_shifted[1]);
      // v_z = v_z
      double v_z = vel[2] / time_step;
      // Write data to the histogram.
      histogram.update(ppos_cyl, std::vector<double>{{v_r, v_phi, v_z}});
    }
    auto hist_tmp = histogram.get_histogram();
    auto tot_count = histogram.get_tot_count();
    for (size_t ind=0; ind < hist_tmp.size(); ++ind) {
      if (hist_tmp[ind] > 0.0) {
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

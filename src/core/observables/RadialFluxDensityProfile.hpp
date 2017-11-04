#ifndef OBSERVABLES_RADIALFLUXDENSITYPROFILE_HPP
#define OBSERVABLES_RADIALFLUXDENSITYPROFILE_HPP

#include "RadialProfileObservable.hpp"
#include "utils.hpp"

namespace Observables {
class RadialFluxDensityProfile : public RadialProfileObservable {
public:
  virtual int actual_calculate(PartCfg &partCfg) override {
    double bin_volume;
    int r_bin, phi_bin, z_bin;
    for (int id : ids()) {
      auto const ppos = ::Vector<3, double>(folded_position(partCfg[id]));
      auto const ppos_shifted = ppos - center;
      auto const ppos_cyl = transform_to_cylinder_coordinates(ppos_shifted);
      r_bin = std::floor((ppos_cyl[0] - min_r) / r_bin_size());
      phi_bin = std::floor((ppos_cyl[1] - min_phi) / phi_bin_size());
      z_bin = std::floor((ppos_cyl[2] - min_z) / z_bin_size());
      bin_volume =
          PI *
          ((min_r + (r_bin + 1) * r_bin_size()) *
               (min_r + (r_bin + 1) * r_bin_size()) -
           (min_r + r_bin * r_bin_size()) * (min_r + r_bin * r_bin_size())) *
          z_bin_size() * phi_bin_size() / (2 * PI);
      if (r_bin >= 0 && r_bin < n_r_bins && phi_bin >= 0 &&
          phi_bin < n_phi_bins && z_bin >= 0 && z_bin < n_z_bins) {
        // Coordinate transform the velocities and divide by time_step to get MD
        // units. v_r = (x * v_x + y * v_y) / sqrt(x^2 + y^2)
        double v_r = (ppos_shifted[0] * partCfg[id].m.v[0] / time_step +
                      ppos_shifted[1] * partCfg[id].m.v[1] / time_step) /
                     std::sqrt(ppos_shifted[0] * ppos_shifted[0] +
                               ppos_shifted[1] * ppos_shifted[1]);
        // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
        double v_phi = (ppos_shifted[0] * partCfg[id].m.v[1] / time_step -
                        ppos_shifted[1] * partCfg[id].m.v[0] / time_step) /
                       (ppos_shifted[0] * ppos_shifted[0] +
                        ppos_shifted[1] * ppos_shifted[1]);
        // v_z = v_z
        double v_z = partCfg[id].m.v[2] / time_step;
        // Write a flat histogram.
        // index calculation: using the following formula for N dimensions:
        //   ind = ind_{N-1} + sum_{j=0}^{N-2} (ind_j * prod_{k=j+1}^{N-1} n_k),
        // where ind is the flattened index, ind_i is the ith unflattened index
        // and n_i is the size of the ith dimension.
        int ind =
            3 * (r_bin * n_phi_bins * n_z_bins + phi_bin * n_z_bins + z_bin);
        last_value[ind + 0] += v_r / bin_volume;
        last_value[ind + 1] += v_phi / bin_volume;
        last_value[ind + 2] += v_z / bin_volume;
      }
    }
    return 0;
  }
  virtual int n_values() const override {
    return 3 * n_r_bins * n_phi_bins * n_z_bins;
  }
};

} // Namespace Observables

#endif

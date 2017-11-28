#ifndef OBSERVABLES_CYLINDRICALFLUXDENSITYPROFILE_HPP
#define OBSERVABLES_CYLINDRICALFLUXDENSITYPROFILE_HPP

#include "CylindricalProfileObservable.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"

namespace Observables {
class CylindricalFluxDensityProfile : public CylindricalProfileObservable {
public:
  virtual int actual_calculate(PartCfg &partCfg) override {
    double bin_volume;
    int r_bin, phi_bin, z_bin;
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
        printf("Velocity: %f %f %f\n", v_r, v_phi, v_z);
        // Write a flat histogram.
        // index calculation: using the following formula for N dimensions:
        //   ind = ind_{N-1} + sum_{j=0}^{N-2} (ind_j * prod_{k=j+1}^{N-1} n_k),
        // where ind is the flattened index, ind_i is the ith unflattened index
        // and n_i is the size of the ith dimension.
        int ind =
            3 * (r_bin * n_phi_bins * n_z_bins + phi_bin * n_z_bins + z_bin);
        if (std::isfinite(v_r)) {
          last_value[ind + 0] += v_r / bin_volume;
        }
        if (std::isfinite(v_phi)) {
          last_value[ind + 1] += v_phi / bin_volume;
        }
        if (std::isfinite(v_z)) {
          last_value[ind + 2] += v_z / bin_volume;
        }
      }
    }
    return 0;
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
        r_bin_size(), phi_bin_size(), z_bin_size()};
    std::array<double, 3> position;
    int index;
    int unravelled_index[4];
    for (auto it = last_value.begin(); it != last_value.end(); it += 3) {
      index = std::distance(last_value.begin(), it);
      ::Utils::unravel_index(len_dims, n_dims, index, unravelled_index);
      position = {
          (static_cast<double>(unravelled_index[0]) + 0.5) * bin_sizes[0],
          (static_cast<double>(unravelled_index[1]) + 0.5) * bin_sizes[1],
          (static_cast<double>(unravelled_index[2]) + 0.5) * bin_sizes[2]};
      m_ofile << position[0] << " " << position[1] << " " << position[2] << " "
              << *it << " " << *(it + 1) << " " << *(it + 2) << "\n";
    }
    m_ofile << std::endl;
  }
};

} // Namespace Observables

#endif

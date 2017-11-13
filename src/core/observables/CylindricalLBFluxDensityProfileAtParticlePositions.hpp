#ifndef OBSERVABLES_CYLINDRICALLBFLUXDENSITYPROFILEATPARTICLEPOSITIONS_HPP
#define OBSERVABLES_CYLINDRICALLBFLUXDENSITYPROFILEATPARTICLEPOSITIONS_HPP

#include "CylindricalProfileObservable.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"
#include "lbgpu.hpp"

namespace Observables {
class CylindricalLBFluxDensityProfileAtParticlePositions : public CylindricalProfileObservable {
public:
  virtual int actual_calculate(PartCfg &partCfg) override {
    double bin_volume;
    int r_bin, phi_bin, z_bin;
    // First collect all positions (since we want to call the LB function to
    // get the fluid velocities only once).
    std::vector<double> ppos(3*ids().size());
    std::vector<double> ppos_cyl(3*ids().size());
    std::vector<double> ppos_shifted(3*ids().size());
    for (int index = 0; index < ids().size(); ++index) {
      auto const ppos_tmp = ::Vector<3, double>(folded_position(partCfg[ids()[index]]));
      ppos[3*index + 0] = ppos_tmp[0];
      ppos[3*index + 1] = ppos_tmp[1];
      ppos[3*index + 2] = ppos_tmp[2];
      auto const ppos_shifted_tmp = ppos_tmp - center;
      ppos_shifted[3*index + 0] = ppos_shifted_tmp[0];
      ppos_shifted[3*index + 1] = ppos_shifted_tmp[1];
      ppos_shifted[3*index + 2] = ppos_shifted_tmp[2];
      auto const ppos_cyl_tmp = Utils::transform_to_cylinder_coordinates(ppos_shifted_tmp);
      ppos_cyl[3*index + 0] = ppos_cyl_tmp[0];
      ppos_cyl[3*index + 1] = ppos_cyl_tmp[1];
      ppos_cyl[3*index + 2] = ppos_cyl_tmp[2];
    }
#ifdef LB_GPU
    std::vector<double> velocities(3*ids().size());
    lb_lbfluid_get_interpolated_velocity_at_positions(ppos.data(), velocities.data(), ids().size());
#endif
    for (int index = 0; index < ids().size(); ++index) {
      r_bin = std::floor((ppos_cyl[3*index + 0] - min_r) / r_bin_size());
      phi_bin = std::floor((ppos_cyl[3*index + 1] - min_phi) / phi_bin_size());
      z_bin = std::floor((ppos_cyl[3*index + 2] - min_z) / z_bin_size());
      bin_volume =
          PI *
          ((min_r + (r_bin + 1) * r_bin_size()) *
               (min_r + (r_bin + 1) * r_bin_size()) -
           (min_r + r_bin * r_bin_size()) * (min_r + r_bin * r_bin_size())) *
          z_bin_size() * phi_bin_size() / (2 * PI);
      if (r_bin >= 0 && r_bin < n_r_bins && phi_bin >= 0 &&
          phi_bin < n_phi_bins && z_bin >= 0 && z_bin < n_z_bins) {
        // Coordinate transform the velocities and divide core velocities by time_step to get MD
        // units.
        // v_r = (x * v_x + y * v_y) / sqrt(x^2 + y^2)
        double v_r = (ppos_shifted[3*index + 0] * velocities[3*index + 0] +
                      ppos_shifted[3*index + 1] * velocities[3*index + 1] ) /
                     std::sqrt(ppos_shifted[3*index + 0] * ppos_shifted[3*index + 0] +
                               ppos_shifted[3*index + 1] * ppos_shifted[3*index + 1]);
        // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
        double v_phi = (ppos_shifted[3*index + 0] * velocities[3*index + 1] -
                        ppos_shifted[3*index + 1] * velocities[3*index + 0] ) /
                       (ppos_shifted[3*index + 0] * ppos_shifted[3*index + 0] +
                        ppos_shifted[3*index + 1] * ppos_shifted[3*index + 1]);
        // v_z = v_z
        double v_z = velocities[3*index + 2];
        // Write a flat histogram.
        // index calculation: using the following formula for N dimensions:
        //   ind = ind_{N-1} + sum_{j=0}^{N-2} (ind_j * prod_{k=j+1}^{N-1} n_k),
        // where ind is the flattened index, ind_i is the ith unflattened index
        // and n_i is the size of the ith dimension.
        int ind =
            3 * (r_bin * n_phi_bins * n_z_bins + phi_bin * n_z_bins + z_bin);
        if (std::isfinite(v_r)) {last_value[ind + 0] += v_r / bin_volume;}
        if (std::isfinite(v_phi)) {last_value[ind + 1] += v_phi / bin_volume;}
        if (std::isfinite(v_z)) {last_value[ind + 2] += v_z / bin_volume;}
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

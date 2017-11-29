/*
  Copyright (C) 2016,2017 The ESPResSo project

  This file is part of ESPResSo.


  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include "CylindricalLBFluxDensityProfileAtParticlePositions.hpp"
#include "lbgpu.hpp"
#include "utils.hpp"

namespace Observables {

std::vector<double> CylindricalLBFluxDensityProfileAtParticlePositions::operator() (
    PartCfg &partCfg) const {
  std::vector<double> res(n_values(), 0.0);
#ifdef LB_GPU
  double bin_volume;
  int r_bin, phi_bin, z_bin;
  // First collect all positions (since we want to call the LB function to
  // get the fluid velocities only once).
  std::vector<double> ppos(3 * ids().size());
  std::vector<double> ppos_cyl(3 * ids().size());
  std::vector<double> ppos_shifted(3 * ids().size());
  for (int index = 0; index < ids().size(); ++index) {
    auto const ppos_tmp =
        ::Vector<3, double>(folded_position(partCfg[ids()[index]]));
    ppos[3 * index + 0] = ppos_tmp[0];
    ppos[3 * index + 1] = ppos_tmp[1];
    ppos[3 * index + 2] = ppos_tmp[2];
    auto const ppos_shifted_tmp = ppos_tmp - center;
    if (axis == "z") {
      ppos_shifted[3 * index + 0] = ppos_shifted_tmp[0];
      ppos_shifted[3 * index + 1] = ppos_shifted_tmp[1];
      ppos_shifted[3 * index + 2] = ppos_shifted_tmp[2];
    } else if (axis == "x") {
      ppos_shifted[3 * index + 0] = -ppos_shifted_tmp[2];
      ppos_shifted[3 * index + 1] = ppos_shifted_tmp[1];
      ppos_shifted[3 * index + 2] = ppos_shifted_tmp[0];
    } else if (axis == "y") {
      ppos_shifted[3 * index + 0] = ppos_shifted_tmp[0];
      ppos_shifted[3 * index + 1] = -ppos_shifted_tmp[2];
      ppos_shifted[3 * index + 2] = ppos_shifted_tmp[1];
    }
    auto const ppos_cyl_tmp =
        Utils::transform_to_cylinder_coordinates(ppos_shifted_tmp);
    ppos_cyl[3 * index + 0] = ppos_cyl_tmp[0];
    ppos_cyl[3 * index + 1] = ppos_cyl_tmp[1];
    ppos_cyl[3 * index + 2] = ppos_cyl_tmp[2];
  }
  std::vector<double> velocities(3 * ids().size());
  lb_lbfluid_get_interpolated_velocity_at_positions(
      ppos.data(), velocities.data(), ids().size());
  for (int index = 0; index < ids().size(); ++index) {
    r_bin = std::floor((ppos_cyl[3 * index + 0] - min_r) / r_bin_size());
    phi_bin = std::floor((ppos_cyl[3 * index + 1] - min_phi) / phi_bin_size());
    z_bin = std::floor((ppos_cyl[3 * index + 2] - min_z) / z_bin_size());
    bin_volume =
        PI *
        ((min_r + (r_bin + 1) * r_bin_size()) *
             (min_r + (r_bin + 1) * r_bin_size()) -
         (min_r + r_bin * r_bin_size()) * (min_r + r_bin * r_bin_size())) *
        z_bin_size() * phi_bin_size() / (2 * PI);
    if (r_bin >= 0 && r_bin < n_r_bins && phi_bin >= 0 &&
        phi_bin < n_phi_bins && z_bin >= 0 && z_bin < n_z_bins) {
      // Coordinate transform the velocities and divide core velocities by
      // time_step to get MD units. v_r = (x * v_x + y * v_y) / sqrt(x^2 + y^2)
      double v_r =
          (ppos_shifted[3 * index + 0] * velocities[3 * index + 0] +
           ppos_shifted[3 * index + 1] * velocities[3 * index + 1]) /
          std::sqrt(ppos_shifted[3 * index + 0] * ppos_shifted[3 * index + 0] +
                    ppos_shifted[3 * index + 1] * ppos_shifted[3 * index + 1]);
      // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
      double v_phi =
          (ppos_shifted[3 * index + 0] * velocities[3 * index + 1] -
           ppos_shifted[3 * index + 1] * velocities[3 * index + 0]) /
          (ppos_shifted[3 * index + 0] * ppos_shifted[3 * index + 0] +
           ppos_shifted[3 * index + 1] * ppos_shifted[3 * index + 1]);
      // v_z = v_z
      double v_z = velocities[3 * index + 2];
      // Write a flat histogram.
      // index calculation: using the following formula for N dimensions:
      //   ind = ind_{N-1} + sum_{j=0}^{N-2} (ind_j * prod_{k=j+1}^{N-1} n_k),
      // where ind is the flattened index, ind_i is the ith unflattened index
      // and n_i is the size of the ith dimension.
      int ind =
          3 * (r_bin * n_phi_bins * n_z_bins + phi_bin * n_z_bins + z_bin);
      if (std::isfinite(v_r)) {
        res[ind + 0] += v_r / bin_volume;
      }
      if (std::isfinite(v_phi)) {
        res[ind + 1] += v_phi / bin_volume;
      }
      if (std::isfinite(v_z)) {
        res[ind + 2] += v_z / bin_volume;
      }
    }
  }
#endif // LB_GPU
  return res;
}
} // namespace Observables

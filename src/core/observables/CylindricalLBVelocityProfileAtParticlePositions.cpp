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
#include "CylindricalLBVelocityProfileAtParticlePositions.hpp"
#include "lbgpu.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"

namespace Observables {

std::vector<double> CylindricalLBVelocityProfileAtParticlePositions::
operator()(PartCfg &partCfg) const {
  std::vector<size_t> n_bins{{static_cast<size_t>(n_r_bins),
                              static_cast<size_t>(n_phi_bins),
                              static_cast<size_t>(n_z_bins)}};
  std::vector<std::pair<double, double>> limits{
      {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
       std::make_pair(min_z, max_z)}};
  Utils::CylindricalHistogram<double> histogram(n_bins, 3, limits);
#ifdef LB_GPU
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
    // Coordinate transform the velocities and divide core velocities by
    // time_step to get MD units. v_r = (x * v_x + y * v_y) / sqrt(x^2 + y^2)
    double v_r =
        (ppos_shifted[3 * index + 0] * velocities[3 * index + 0] +
         ppos_shifted[3 * index + 1] * velocities[3 * index + 1]) /
        std::sqrt(ppos_shifted[3 * index + 0] * ppos_shifted[3 * index + 0] +
                  ppos_shifted[3 * index + 1] * ppos_shifted[3 * index + 1]);
    // v_phi = (x * v_y - y * v_x ) / (x^2 + y^2)
    double v_phi = (ppos_shifted[3 * index + 0] * velocities[3 * index + 1] -
                    ppos_shifted[3 * index + 1] * velocities[3 * index + 0]) /
                   (ppos_shifted[3 * index + 0] * ppos_shifted[3 * index + 0] +
                    ppos_shifted[3 * index + 1] * ppos_shifted[3 * index + 1]);
    // v_z = v_z
    double v_z = velocities[3 * index + 2];
    histogram.update(
      std::vector<double>{{ppos_cyl[3 * index + 0], ppos_cyl[3 * index + 1],
                             ppos_cyl[3 * index + 2]}},
      std::vector<double>{{v_r, v_phi, v_z}});
  }
  auto hist_tmp = histogram.get_histogram();
  auto tot_count = histogram.get_tot_count();
  for (size_t ind=0; ind < hist_tmp.size(); ++ind) {
    if (hist_tmp[ind] > 0.0) {
      hist_tmp[ind] /= tot_count[ind];
    }
  }
  return hist_tmp;
#endif // LB_GPU
#ifndef LB_GPU
  return histogram.get_histogram();
#endif
}
} // namespace Observables

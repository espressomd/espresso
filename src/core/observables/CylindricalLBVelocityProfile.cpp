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
#include "CylindricalLBVelocityProfile.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"

namespace Observables {

std::vector<double> CylindricalLBVelocityProfile::
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
  std::vector<double> velocities(m_node_positions.size());
  lb_lbfluid_get_interpolated_velocity_at_positions(
      m_node_positions.data(), velocities.data(), m_node_positions.size() / 3);
  ::Vector<3, double> pos_shifted, pos_cyl, velocity;
  for (size_t ind = 0; ind < m_node_positions.size() / 3; ind+=3) {
    //printf("pos: %f %f %f\n", m_node_positions[ind+0], m_node_positions[ind+1], m_node_positions[ind+2]);
    pos_shifted = {{m_node_positions[ind + 0] - center[0],
                    m_node_positions[ind + 1] - center[1],
                    m_node_positions[ind + 2] - center[2]}};
    pos_cyl = Utils::transform_pos_to_cylinder_coordinates(pos_shifted, axis);
    //printf("pos: %f %f %f\n", pos_shifted[0], pos_shifted[1], pos_shifted[2]);
    //printf("pos_cyl: %f %f %f\n", pos_cyl[0], pos_cyl[1], pos_cyl[2]);
    velocity = {{velocities[ind + 0],
                 velocities[ind + 1],
                 velocities[ind + 2]}};
    histogram.update(Utils::transform_pos_to_cylinder_coordinates(
                         pos_shifted,
                         axis),
                     Utils::transform_vel_to_cylinder_coordinates(
                         velocity,
                         axis, pos_shifted));
  }
  auto hist_tmp = histogram.get_histogram();
  auto tot_count = histogram.get_tot_count();
  for (size_t ind = 0; ind < hist_tmp.size(); ++ind) {
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

void CylindricalLBVelocityProfile::calculate_node_positions() {
  m_node_positions.clear();
  if (lbpar_gpu.agrid == 0.0)
    throw std::runtime_error("LB GPU fluid has to be initialized first.");
  size_t n_nodes_x = box_l[0] / lbpar_gpu.agrid;
  size_t n_nodes_y = box_l[1] / lbpar_gpu.agrid;
  size_t n_nodes_z = box_l[2] / lbpar_gpu.agrid;
  for (size_t x = 0; x < n_nodes_x; ++x) {
    for (size_t y = 0; y < n_nodes_y; ++y) {
      for (size_t z = 0; z < n_nodes_z; ++z) {
        m_node_positions.push_back(0.5 * lbpar_gpu.agrid + (x * lbpar_gpu.agrid));
        m_node_positions.push_back(0.5 * lbpar_gpu.agrid + (y * lbpar_gpu.agrid));
        m_node_positions.push_back(0.5 * lbpar_gpu.agrid + (z * lbpar_gpu.agrid));
      }
    }
  }
}

} // namespace Observables

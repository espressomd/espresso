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
#include "lb.hpp"
#include "CylindricalLBVelocityProfile.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"

namespace Observables {

std::vector<double> CylindricalLBVelocityProfile::
operator()(PartCfg &partCfg) const {
  std::array<size_t, 3> n_bins{{static_cast<size_t>(n_r_bins),
                                static_cast<size_t>(n_phi_bins),
                                static_cast<size_t>(n_z_bins)}};
  std::array<std::pair<double, double>, 3> limits{
      {std::make_pair(min_r, max_r), std::make_pair(min_phi, max_phi),
       std::make_pair(min_z, max_z)}};
  Utils::CylindricalHistogram<double, 3> histogram(n_bins, 3, limits);
  // First collect all positions (since we want to call the LB function to
  // get the fluid velocities only once).
  std::vector<double> velocities(m_sample_positions.size());
  if (lattice_switch & LATTICE_LB_GPU) {
#if defined(LB_GPU)
    lb_lbfluid_get_interpolated_velocity_at_positions(
        m_sample_positions.data(), velocities.data(),
        m_sample_positions.size() / 3);
#endif
  } else if (lattice_switch & LATTICE_LB) {
#if defined(LB)
    for (size_t ind=0; ind < m_sample_positions.size(); ind +=3) {
      Vector3d pos_tmp = {m_sample_positions[ind + 0],
                           m_sample_positions[ind + 1],
                           m_sample_positions[ind + 2]};
      lb_lbfluid_get_interpolated_velocity(pos_tmp, &(velocities[ind + 0]));
    }
#endif
  } else {
    return histogram.get_histogram();
  }
  for (size_t ind = 0; ind < m_sample_positions.size(); ind += 3) {
    const Vector3d pos_shifted = {{m_sample_positions[ind + 0] - center[0],
                                   m_sample_positions[ind + 1] - center[1],
                                   m_sample_positions[ind + 2] - center[2]}};
    const Vector3d pos_cyl =
        Utils::transform_pos_to_cylinder_coordinates(pos_shifted, axis);
    const Vector3d velocity = {
        {velocities[ind + 0], velocities[ind + 1], velocities[ind + 2]}};
    histogram.update(pos_cyl, Utils::transform_vel_to_cylinder_coordinates(
                                  velocity, axis, pos_shifted));
  }
  auto hist_tmp = histogram.get_histogram();
  auto tot_count = histogram.get_tot_count();
  for (size_t ind = 0; ind < hist_tmp.size(); ++ind) {
    if (tot_count[ind] == 0 and not allow_empty_bins)
      throw std::runtime_error("Decrease sampling delta(s), bin without hit "
                               "found!");
    if (tot_count[ind] > 0) {
      hist_tmp[ind] /= tot_count[ind];
    }
  }
  return hist_tmp;
}

} // namespace Observables

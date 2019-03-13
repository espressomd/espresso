/*
  Copyright (C) 2016-2018 The ESPResSo project

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
#include "grid_based_algorithms/lb_interface.hpp"
#include "grid_based_algorithms/lb_interpolation.hpp"
#include "utils.hpp"
#include "utils/Histogram.hpp"
#include "utils/coordinate_transformation.hpp"
#include <boost/range/algorithm/transform.hpp>

namespace Observables {

std::vector<double> CylindricalLBVelocityProfileAtParticlePositions::
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
  std::vector<Vector3d> folded_positions(ids().size());
  boost::transform(ids(), folded_positions.begin(),
                   [&partCfg](int id) { return folded_position(partCfg[id]); });

  std::vector<Vector3d> velocities(ids().size());
#if defined(LB) || defined(LB_GPU)
  boost::transform(
      folded_positions, velocities.begin(), [](const Vector3d &pos) {
        return lb_lbinterpolation_get_interpolated_velocity_global(pos) *
               lb_lbfluid_get_lattice_speed();
      });
#endif
  for (auto &p : folded_positions)
    p -= center;
  for (int ind = 0; ind < ids().size(); ++ind) {
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

} // namespace Observables

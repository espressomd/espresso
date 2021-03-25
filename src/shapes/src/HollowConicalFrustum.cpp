/*
 * Copyright (C) 2010-2019 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <shapes/HollowConicalFrustum.hpp>

#include <utils/Vector.hpp>
#include <utils/math/abs.hpp>
#include <utils/math/coordinate_transformation.hpp>

namespace Shapes {
void HollowConicalFrustum::calculate_dist(const Utils::Vector3d &pos,
                                          double &dist,
                                          Utils::Vector3d &vec) const {

  // transform given position to cylindrical coordinates in the reference frame
  // of the cone
  auto const v = pos - m_cyl_transform_params->center();
  auto const pos_cyl = Utils::transform_coordinate_cartesian_to_cylinder(
      v, m_cyl_transform_params->axis(), m_cyl_transform_params->orientation());

  auto project_on_line = [](auto const vec, auto const line_start, auto const line_director) {
    return line_start + line_director *((vec-line_start)*line_director);
  };

  Utils::Vector3d pos_closest;
  if (Utils::abs(pos_cyl[1])>=m_central_angle/2.){
    // Go to 2d, find the projection onto the cone mantle
    auto const pos_2d = Utils::Vector2d{{pos_cyl[0], pos_cyl[2]}};
    auto const r1_endpoint = Utils::Vector2d{{m_r1,m_length/2.}};
    auto const r2_endpoint = Utils::Vector2d{{m_r2,-m_length/2.}};
    auto const line_director = (r2_endpoint-r1_endpoint).normalized();
    auto closest_point_2d = project_on_line(pos_2d, r1_endpoint, line_director);

    // correct projection if it is outside the frustum
    if (Utils::abs(closest_point_2d[1])>m_length/2.){
      bool at_r1 = closest_point_2d[1]>0;
      closest_point_2d[0] = at_r1 ? m_r1 : m_r2;
      closest_point_2d[1] = at_r1 ? m_length/2. : -m_length/2.;
    }

    // Go back to cartesian coordinates of the box frame
    pos_closest = Utils::transform_coordinate_cylinder_to_cartesian(
        {closest_point_2d[0], pos_cyl[1], closest_point_2d[1]}, m_cyl_transform_params->axis(), m_cyl_transform_params->orientation()) +
                  m_cyl_transform_params->center();

  }
  else{
    // We cannot go to 2d because the central-angle-gap breaks rotational symmetry,
    // so we have to get the line endpoints of the closer edge in 3d cartesian coordinates (but still in the reference frame of the HCF)

    // Cannot use Utils::sgn because of pos_cyl[1]==0 corner case
    auto const endpoint_angle = pos_cyl[1]>=0 ? m_central_angle/2. : -m_central_angle/2.;
    auto const r1_endpoint = Utils::transform_coordinate_cylinder_to_cartesian(Utils::Vector3d{{m_r1,endpoint_angle,m_length/2.}});
    auto const r2_endpoint = Utils::transform_coordinate_cylinder_to_cartesian(Utils::Vector3d{{m_r2,endpoint_angle,-m_length/2.}});
    auto const line_director = (r2_endpoint-r1_endpoint).normalized();
    auto const pos_hcf_frame = Utils::transform_coordinate_cylinder_to_cartesian(pos_cyl);
    auto pos_closest_hcf_frame = project_on_line(pos_hcf_frame, r1_endpoint, line_director);

    // Go back to cylindrical coordinates (HCF reference frame), here we can apply the capping at z = plusminus l/2.
    auto pos_closest_hcf_cyl = Utils::transform_coordinate_cartesian_to_cylinder(pos_closest_hcf_frame);
    if (Utils::abs(pos_closest_hcf_cyl[2])>m_length/2.){
      bool at_r1 = pos_closest_hcf_cyl[2]>0.;
      pos_closest_hcf_cyl[0] = at_r1 ? m_r1 : m_r2;
      pos_closest_hcf_cyl[2] = at_r1 ? m_length/2. : -m_length/2.;
    }

    // Finally, go to cartesian coordinates of the box frame
    pos_closest = Utils::transform_coordinate_cylinder_to_cartesian(
        pos_closest_hcf_cyl, m_cyl_transform_params->axis(), m_cyl_transform_params->orientation()) +
                  m_cyl_transform_params->center();
  }

  auto const u = (pos - pos_closest).normalize();
  auto const d = (pos - pos_closest).norm() - 0.5 * m_thickness;
  dist = d * m_direction;
  vec = d * u;
}
} // namespace Shapes

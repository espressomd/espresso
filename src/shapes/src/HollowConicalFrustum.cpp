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

  // Use the rotational symmetry of the cone: Transformation of pos to the
  // cylindrical coordinates in the frame of the cone.
  auto const v = pos - m_cyl_transform_params->center();
  auto const pos_cyl = Utils::transform_coordinate_cartesian_to_cylinder(
      v, m_cyl_transform_params->axis(), m_cyl_transform_params->orientation());

  auto project_on_line = [](auto const vec, auto const line_start,
                            auto const line_director) {
    return line_start + line_director * ((vec - line_start) * line_director);
  };

  // The point on the frustum closest to pos will be determined, dist and vec
  // follow trivially.
  Utils::Vector3d pos_closest;

  if (Utils::abs(pos_cyl[1]) >= m_central_angle / 2.) {
    /* First case: pos is not in the gap region defined by central_angle.
     * Here, the problem reduces to 2D, because pos_closest lies in the plane
     * defined by axis and pos as shown in the figure below:
     *
     *     r1
     * * ----->X r1_endpoint
     * ^ a      \
     * | x       \
     * | i        \
     * | s         \                       *pos
     * X center     \
     * |             \
     * |              X pos_closest_2d (to be determined)
     * |               \
     * |       r2       \
     * ----------------->X r2_endpoint
     *
     * pos_closest_2d is the projection of pos on the line defined by
     * r1_endpoint and r2_endpoint. The real pos_closest is then calculated
     * using the phi-coordinate of pos_cyl and transformation back to the box
     * frame.
     */
    auto const pos_2d = Utils::Vector2d{{pos_cyl[0], pos_cyl[2]}};
    auto const r1_endpoint = Utils::Vector2d{{m_r1, m_length / 2.}};
    auto const r2_endpoint = Utils::Vector2d{{m_r2, -m_length / 2.}};
    auto const line_director = (r2_endpoint - r1_endpoint).normalized();
    auto pos_closest_2d = project_on_line(pos_2d, r1_endpoint, line_director);

    // If the projection is outside of the frustum, pos_closest_2d is one of the
    // endpoints
    if (Utils::abs(pos_closest_2d[1]) > m_length / 2.) {
      pos_closest_2d = pos_closest_2d[1] > 0 ? r1_endpoint : r2_endpoint;
    }

    pos_closest = Utils::transform_coordinate_cylinder_to_cartesian(
                      {pos_closest_2d[0], pos_cyl[1], pos_closest_2d[1]},
                      m_cyl_transform_params->axis(),
                      m_cyl_transform_params->orientation()) +
                  m_cyl_transform_params->center();

  } else {
    /* If pos is in the gap region, we cannot go to 2d or cylindrical
     * coordinates, because the central-angle-gap breaks rotational symmetry.
     * Instead, we parametrise the closer edge of the frustum cutout by its
     * endpoints in 3D cartesian coordinates, but still in the reference frame
     * of the HCF. The projection onto this edge to find pos_closest is then the
     * same procedure as in the previous case.
     */

    // Cannot use Utils::sgn because of pos_cyl[1]==0 corner case
    auto const endpoint_angle =
        pos_cyl[1] >= 0 ? m_central_angle / 2. : -m_central_angle / 2.;
    auto const r1_endpoint = Utils::transform_coordinate_cylinder_to_cartesian(
        Utils::Vector3d{{m_r1, endpoint_angle, m_length / 2.}});
    auto const r2_endpoint = Utils::transform_coordinate_cylinder_to_cartesian(
        Utils::Vector3d{{m_r2, endpoint_angle, -m_length / 2.}});
    auto const line_director = (r2_endpoint - r1_endpoint).normalized();
    auto const pos_hcf_frame =
        Utils::transform_coordinate_cylinder_to_cartesian(pos_cyl);
    auto pos_closest_hcf_frame =
        project_on_line(pos_hcf_frame, r1_endpoint, line_director);

    if (Utils::abs(pos_closest_hcf_frame[2]) > m_length / 2.) {
      pos_closest_hcf_frame =
          pos_closest_hcf_frame[2] > 0. ? r1_endpoint : r2_endpoint;
    }

    // Now we transform pos_closest_hcf_frame back to the box frame.
    auto const hcf_frame_y_axis = Utils::vector_product(
        m_cyl_transform_params->axis(), m_cyl_transform_params->orientation());
    pos_closest =
        Utils::basis_change(m_cyl_transform_params->orientation(),
                            hcf_frame_y_axis, m_cyl_transform_params->axis(),
                            pos_closest_hcf_frame, true) +
        m_cyl_transform_params->center();
  }

  auto const u = (pos - pos_closest).normalize();
  auto const d = (pos - pos_closest).norm() - 0.5 * m_thickness;
  dist = d * m_direction;
  vec = d * u;
}
} // namespace Shapes

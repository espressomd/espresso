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

  /* Transform pos to the frame of the frustum (origin = center, z = axis, x =
   * orientation, y = cross(z,x) ). Get the angle relative to orientation to
   * determine whether pos is in the gap defined by m_central_angle or not.
   */
  auto const hcf_frame_y_axis = Utils::vector_product(
      m_cyl_transform_params->axis(), m_cyl_transform_params->orientation());
  auto const pos_hcf_frame = Utils::basis_change(
      m_cyl_transform_params->orientation(), hcf_frame_y_axis,
      m_cyl_transform_params->axis(), pos - m_cyl_transform_params->center());
  auto const pos_phi =
      Utils::transform_coordinate_cartesian_to_cylinder(pos_hcf_frame)[1];

  /* The closest point of the frustum to pos will be determined by projection
   * onto a line. Which line, depends on where pos is relative to the frustum
   * gap. If pos is not within the gap region of central_angle, the closest
   * point will be in the plane defined by axis and pos_cyl, as shown in the
   * figure below:
   *
   *     r1
   * * ----->X r1_endpoint
   * ^ a      \
   * | x       \
   * | i        \
   * | s         \                       *pos_hcf_frame
   * X center     \
   * |             \
   * |              X pos_closest (to be determined)
   * |               \
   * |       r2       \
   * ----------------->X r2_endpoint
   *
   * In this case, the line is the intersection between the frustum and the
   * plane of interest. The endpoints of the line are determined by the radii of
   * the frustum, its length and the phi-coordinate of the plane.
   *
   * If pos is in the gap region, the projection must be made onto the closest
   * edge (imagine pos is out-of-plane in the figure). In this case, for the
   * endpoints we do not use the phi-coordinate of pos_cyl, but instead the
   * phi-coordinate of the closer gap edge.
   */

  auto endpoint_angle = pos_phi;
  if (Utils::abs(pos_phi) < m_central_angle / 2.) {
    // Cannot use Utils::sgn because of pos_phi==0 corner case
    endpoint_angle =
        pos_phi > 0. ? m_central_angle / 2. : -m_central_angle / 2.;
  }

  auto const r1_endpoint = Utils::transform_coordinate_cylinder_to_cartesian(
      Utils::Vector3d{{m_r1, endpoint_angle, m_length / 2.}});
  auto const r2_endpoint = Utils::transform_coordinate_cylinder_to_cartesian(
      Utils::Vector3d{{m_r2, endpoint_angle, -m_length / 2.}});
  auto const line_director = (r2_endpoint - r1_endpoint).normalized();

  auto project_on_line = [](auto const vec, auto const line_start,
                            auto const line_director) {
    return line_start + line_director * ((vec - line_start) * line_director);
  };

  auto pos_closest_hcf_frame =
      project_on_line(pos_hcf_frame, r1_endpoint, line_director);

  /* It can be that the projection onto the (infinite) line is outside the
   * frustum. In that case, the closest point is actually one of the endpoints.
   */
  if (Utils::abs(pos_closest_hcf_frame[2]) > m_length / 2.) {
    pos_closest_hcf_frame =
        pos_closest_hcf_frame[2] > 0. ? r1_endpoint : r2_endpoint;
  }

  // calculate distance and distance vector respecting thickness and direction
  auto const u = (pos_hcf_frame - pos_closest_hcf_frame).normalize();
  auto const d =
      (pos_hcf_frame - pos_closest_hcf_frame).norm() - 0.5 * m_thickness;

  // dist does not depend on reference frame, it can be calculated in the hcf
  // frame.
  dist = d * m_direction;

  // vec must be rotated back to the original frame
  auto const vec_hcf_frame = d * u;
  vec = Utils::basis_change(m_cyl_transform_params->orientation(),
                            hcf_frame_y_axis, m_cyl_transform_params->axis(),
                            vec_hcf_frame, true);
}
} // namespace Shapes

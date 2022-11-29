/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include <shapes/Rhomboid.hpp>

#include <utils/Vector.hpp>

#include <functional>

namespace Shapes {
void Rhomboid::calculate_dist(const Utils::Vector3d &pos, double &dist,
                              Utils::Vector3d &vec) const {

  auto const le = std::less_equal<double>();
  auto const ge = std::greater_equal<double>();
  auto const lt = std::less<double>();
  auto const gt = std::greater<double>();

  // calculate vectors and scalars that are going to be used frequently

  auto const axb = vector_product(m_a, m_b);
  auto const bxc = vector_product(m_b, m_c);
  auto const axc = vector_product(m_a, m_c);

  auto const a_dot_bxc = m_a * bxc;
  auto const b_dot_axc = m_b * axc;
  auto const c_dot_axb = m_c * axb;

  auto const dpos = pos - m_pos;

  // compute distance from the rhomboid corners, edges and faces using linear
  // combinations of the rhomboid edge vectors

  auto const corner = [this, &vec, &dist, a = bxc / a_dot_bxc,
                       b = axc / b_dot_axc,
                       c = axb / c_dot_axb](auto op1, auto op2, auto op3,
                                            Utils::Vector3d const &d) {
    /* coefficients A, B, C tell whether ppos lies within a cone defined
     * by pos and the adjacent edges */
    auto const A = a * d;
    auto const B = b * d;
    auto const C = c * d;
    if (op1(A, 0) and op2(B, 0) and op3(C, 0)) {
      vec = d;
      dist = m_direction * vec.norm();
      return true;
    }
    return false;
  };

  if ( // check for cone at m_pos
      corner(le, le, le, dpos) ||
      // check for cone at m_pos+a
      corner(ge, le, le, dpos - m_a) ||
      // check for cone at m_pos+b
      corner(le, ge, le, dpos - m_b) ||
      // check for cone at m_pos+c
      corner(le, le, ge, dpos - m_c) ||
      // check for cone at m_pos+a+b
      corner(ge, ge, le, dpos - m_a - m_b) ||
      // check for cone at m_pos+a+c
      corner(ge, le, ge, dpos - m_a - m_c) ||
      // check for cone at m_pos+b+c
      corner(le, ge, ge, dpos - m_b - m_c) ||
      // check for cone at m_pos+a+b+c
      corner(ge, ge, ge, dpos - m_a - m_b - m_c))
    return;

  auto const edge = [this, &vec, &dist](auto op1, auto op2,
                                        Utils::Vector3d const &d,
                                        Utils::Vector3d const &axis1,
                                        double const dir1_dot_axis1,
                                        Utils::Vector3d const &axis2,
                                        double const dir2_dot_axis2,
                                        Utils::Vector3d const &edge) {
    auto const A = (d * axis1) / dir1_dot_axis1;
    auto const B = (d * axis2) / dir2_dot_axis2;
    if (op1(A, 0) and op2(B, 0)) {
      auto const tmp = (d * edge) / edge.norm2();
      vec = d - edge * tmp;
      dist = m_direction * vec.norm();
      return true;
    }
    return false;
  };

  if ( // check for prism at edge m_pos, a
      edge(le, le, dpos, axc, b_dot_axc, axb, c_dot_axb, m_a) ||
      // check for prism at edge m_pos, b
      edge(le, le, dpos, bxc, a_dot_bxc, axb, c_dot_axb, m_b) ||
      // check for prism at edge m_pos, c
      edge(le, le, dpos, bxc, a_dot_bxc, axc, b_dot_axc, m_c) ||
      // check for prism at edge m_pos+a, b
      edge(ge, le, dpos - m_a, bxc, a_dot_bxc, axb, c_dot_axb, m_b) ||
      // check for prism at edge m_pos+a, c
      edge(ge, le, dpos - m_a, bxc, a_dot_bxc, axc, b_dot_axc, m_c) ||
      // check for prism at edge m_pos+b+c, c
      edge(le, ge, dpos - m_b - m_c, bxc, a_dot_bxc, axc, b_dot_axc, m_c) ||
      // check for prism at edge m_pos+b+c, b
      edge(le, ge, dpos - m_b - m_c, bxc, a_dot_bxc, axb, c_dot_axb, m_b) ||
      // check for prism at edge m_pos+b+c, a
      edge(ge, ge, dpos - m_b - m_c, axc, b_dot_axc, axb, c_dot_axb, m_a) ||
      // check for prism at edge m_pos+a+b, a
      edge(ge, le, dpos - m_a - m_b, axc, b_dot_axc, axb, c_dot_axb, m_a) ||
      // check for prism at edge m_pos+a+b, c
      edge(ge, ge, dpos - m_a - m_b, bxc, a_dot_bxc, axc, b_dot_axc, m_c) ||
      // check for prism at edge m_pos+a+c, a
      edge(le, ge, dpos - m_a - m_c, axc, b_dot_axc, axb, c_dot_axb, m_a) ||
      // check for prism at edge m_pos+a+c, b
      edge(ge, ge, dpos - m_a - m_c, bxc, a_dot_bxc, axb, c_dot_axb, m_b))
    return;

  auto const face_outside =
      [this, &vec, &dist](auto op1, auto op2, Utils::Vector3d const &distance,
                          Utils::Vector3d const &axis,
                          double const dir_dot_axis, int sign) {
        auto d = distance * axis;
        if (op1(dir_dot_axis, 0)) {
          d *= -1;
        }
        if (d >= 0) {
          auto const tmp = axis.norm();
          d /= tmp;
          dist = d * m_direction;
          if (op2(dir_dot_axis, 0)) {
            d *= -1;
          }
          vec = (sign * d / tmp) * axis;
          return true;
        }
        return false;
      };

  if ( // check for face with normal -axb
      face_outside(gt, lt, dpos, axb, c_dot_axb, -1) ||
      // calculate distance to face with normal axc
      face_outside(gt, gt, dpos, axc, b_dot_axc, +1) ||
      // calculate distance to face with normal -bxc
      face_outside(gt, lt, dpos, bxc, a_dot_bxc, -1) ||
      // calculate distance to face with normal axb
      face_outside(lt, lt, dpos - m_a - m_b - m_c, axb, c_dot_axb, +1) ||
      // calculate distance to face with normal -axc
      face_outside(lt, gt, dpos - m_a - m_b - m_c, axc, b_dot_axc, -1) ||
      // calculate distance to face with normal bxc
      face_outside(lt, lt, dpos - m_a - m_b - m_c, bxc, a_dot_bxc, +1))
    return;

  // ppos lies within rhomboid.
  // Find nearest wall for interaction (test all 6 possibilities).

  // calculate distance to face with normal -axb
  {
    auto d = dpos * axb;
    if (c_dot_axb > 0.0) {
      d *= -1;
    }
    auto const tmp = axb.norm();
    d /= tmp;
    dist = d * m_direction;
    if (c_dot_axb < 0.0) {
      d *= -1;
    }
    vec = (-d / tmp) * axb;
  }

  auto const face_inside =
      [this, &vec, &dist](auto op1, auto op2, Utils::Vector3d const &distance,
                          Utils::Vector3d const &axis,
                          double const dir_dot_axis, int sign) {
        auto d = distance * axis;
        if (op1(dir_dot_axis, 0)) {
          d *= -1;
        }
        auto const tmp = axis.norm();
        d /= tmp;
        if (std::abs(d) < std::abs(dist)) {
          dist = d * m_direction;
          if (op2(dir_dot_axis, 0)) {
            d *= -1;
          }
          vec = (sign * d / tmp) * axis;
        }
      };

  // calculate distance to face with normal axc
  face_inside(gt, gt, dpos, axc, b_dot_axc, +1);
  // calculate distance to face with normal -bxc
  face_inside(gt, lt, dpos, bxc, a_dot_bxc, -1);
  // calculate distance to face with normal axb
  face_inside(lt, lt, dpos - m_a - m_b - m_c, axb, c_dot_axb, +1);
  // calculate distance to face with normal -axc
  face_inside(lt, gt, dpos - m_a - m_b - m_c, axc, b_dot_axc, -1);
  // calculate distance to face with normal bxc
  face_inside(lt, lt, dpos - m_a - m_b - m_c, bxc, a_dot_bxc, +1);
}

} // namespace Shapes

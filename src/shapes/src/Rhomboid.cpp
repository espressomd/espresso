/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include <cmath>

using namespace std;

namespace Shapes {
void Rhomboid::calculate_dist(const Utils::Vector3d &pos, double &dist,
                              Utils::Vector3d &vec) const {

  using le = std::less_equal<double>;
  using ge = std::greater_equal<double>;
  using lt = std::less<double>;
  using gt = std::greater<double>;

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

#ifndef DOXYGEN
#define DISTANCE_FROM_CORNER(op1, op2, op3, distance)                          \
  {                                                                            \
    auto const d = distance;                                                   \
    /* coefficients A, B, C tell whether ppos lies within a cone defined */    \
    /* by pos and the adjacent edges */                                        \
    auto const A = (d * bxc) / a_dot_bxc;                                      \
    auto const B = (d * axc) / b_dot_axc;                                      \
    auto const C = (d * axb) / c_dot_axb;                                      \
    if (op1{}(A, 0) & op2{}(B, 0) & op3{}(C, 0)) {                             \
      vec = d;                                                                 \
      dist = m_direction * vec.norm();                                         \
      return;                                                                  \
    }                                                                          \
  }

  // check for cone at pos+a
  DISTANCE_FROM_CORNER(le, le, le, dpos);
  // check for cone at pos+a
  DISTANCE_FROM_CORNER(ge, le, le, dpos - m_a);
  // check for cone at pos+b
  DISTANCE_FROM_CORNER(le, ge, le, dpos - m_b);
  // check for cone at pos+c
  DISTANCE_FROM_CORNER(le, le, ge, dpos - m_c);
  // check for cone at m_pos+a+b
  DISTANCE_FROM_CORNER(ge, ge, le, dpos - m_a - m_b);
  // check for cone at m_pos+a+c
  DISTANCE_FROM_CORNER(ge, le, ge, dpos - m_a - m_c);
  // check for cone at m_pos+b+c
  DISTANCE_FROM_CORNER(le, ge, ge, dpos - m_b - m_c);
  // check for cone at m_pos+a+b+c
  DISTANCE_FROM_CORNER(ge, ge, ge, dpos - m_a - m_b - m_c);

#define DISTANCE_FROM_EDGE(op1, op2, distance, axis1, dir1, axis2, dir2, edge) \
  {                                                                            \
    auto const d = distance;                                                   \
    auto const A = (d * axis1) / dir1##_dot_##axis1;                           \
    auto const B = (d * axis2) / dir2##_dot_##axis2;                           \
    if (op1{}(A, 0) & op2{}(B, 0)) {                                           \
      auto const tmp = (d * edge) / edge.norm2();                              \
      vec = d - edge * tmp;                                                    \
      dist = m_direction * vec.norm();                                         \
      return;                                                                  \
    }                                                                          \
  }

  // check for prism at edge m_pos, a
  DISTANCE_FROM_EDGE(le, le, dpos, axc, b, axb, c, m_a);
  // check for prism at edge m_pos, b
  DISTANCE_FROM_EDGE(le, le, dpos, bxc, a, axb, c, m_b);
  // check for prism at edge m_pos, c
  DISTANCE_FROM_EDGE(le, le, dpos, bxc, a, axc, b, m_c);
  // check for prism at edge m_pos+a, b
  DISTANCE_FROM_EDGE(ge, le, dpos - m_a, bxc, a, axb, c, m_b);
  // check for prism at edge m_pos+a, c
  DISTANCE_FROM_EDGE(ge, le, dpos - m_a, bxc, a, axc, b, m_c);
  // check for prism at edge m_pos+b+c, c
  DISTANCE_FROM_EDGE(le, ge, dpos - m_b - m_c, bxc, a, axc, b, m_c);
  // check for prism at edge m_pos+b+c, b
  DISTANCE_FROM_EDGE(le, ge, dpos - m_b - m_c, bxc, a, axb, c, m_b);
  // check for prism at edge m_pos+b+c, a
  DISTANCE_FROM_EDGE(ge, ge, dpos - m_b - m_c, axc, b, axb, c, m_a);
  // check for prism at edge m_pos+a+b, a
  DISTANCE_FROM_EDGE(ge, le, dpos - m_a - m_b, axc, b, axb, c, m_a);
  // check for prism at edge m_pos+a+b, c
  DISTANCE_FROM_EDGE(ge, ge, dpos - m_a - m_b, bxc, a, axc, b, m_c);
  // check for prism at edge m_pos+a+c, a
  DISTANCE_FROM_EDGE(le, ge, dpos - m_a - m_c, axc, b, axb, c, m_a);
  // check for prism at edge m_pos+a+c, b
  DISTANCE_FROM_EDGE(ge, ge, dpos - m_a - m_c, bxc, a, axb, c, m_b);

#define DISTANCE_FROM_FACE(op1, op2, distance, axis, dir, sign)                \
  {                                                                            \
    auto d = (distance)*axis;                                                  \
    if (op1{}(dir##_dot_##axis, 0)) {                                          \
      d *= -1;                                                                 \
    }                                                                          \
    if (d >= 0) {                                                              \
      auto const tmp = axis.norm();                                            \
      d /= tmp;                                                                \
      dist = d * m_direction;                                                  \
      if (op2{}(dir##_dot_##axis, 0)) {                                        \
        d *= -1;                                                               \
      }                                                                        \
      vec = (sign * d / tmp) * axis;                                           \
      return;                                                                  \
    }                                                                          \
  }

  // check for face with normal -axb
  DISTANCE_FROM_FACE(gt, lt, dpos, axb, c, -1);
  // calculate distance to face with normal axc
  DISTANCE_FROM_FACE(gt, gt, dpos, axc, b, +1);
  // calculate distance to face with normal -bxc
  DISTANCE_FROM_FACE(gt, lt, dpos, bxc, a, -1);
  // calculate distance to face with normal axb
  DISTANCE_FROM_FACE(lt, lt, dpos - m_a - m_b - m_c, axb, c, +1);
  // calculate distance to face with normal -axc
  DISTANCE_FROM_FACE(lt, gt, dpos - m_a - m_b - m_c, axc, b, -1);
  // calculate distance to face with normal bxc
  DISTANCE_FROM_FACE(lt, lt, dpos - m_a - m_b - m_c, bxc, a, +1);

  // ppos lies within rhomboid.
  // Find nearest wall for interaction (test all 6 possibilities).

  // check for face with normal -axb

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

#define DISTANCE_FROM_FACE_INSIDE(op1, op2, distance, axis, dir, sign)         \
  {                                                                            \
    auto d = (distance)*axis;                                                  \
    if (op1{}(dir##_dot_##axis, 0)) {                                          \
      d *= -1;                                                                 \
    }                                                                          \
    auto const tmp = axis.norm();                                              \
    d /= tmp;                                                                  \
    if (abs(d) < abs(dist)) {                                                  \
      dist = d * m_direction;                                                  \
      if (op2{}(dir##_dot_##axis, 0)) {                                        \
        d *= -1;                                                               \
      }                                                                        \
      vec = (sign * d / tmp) * axis;                                           \
    }                                                                          \
  }

  // calculate distance to face with normal axc
  DISTANCE_FROM_FACE_INSIDE(gt, gt, dpos, axc, b, +1);
  // calculate distance to face with normal -bxc
  DISTANCE_FROM_FACE_INSIDE(gt, lt, dpos, bxc, a, -1);
  // calculate distance to face with normal axb
  DISTANCE_FROM_FACE_INSIDE(lt, lt, dpos - m_a - m_b - m_c, axb, c, +1);
  // calculate distance to face with normal -axc
  DISTANCE_FROM_FACE_INSIDE(lt, gt, dpos - m_a - m_b - m_c, axc, b, -1);
  // calculate distance to face with normal bxc
  DISTANCE_FROM_FACE_INSIDE(lt, lt, dpos - m_a - m_b - m_c, bxc, a, +1);
#endif // ifndef DOXYGEN
}

} // namespace Shapes

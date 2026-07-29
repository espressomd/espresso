/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include "algorithm/periodic_fold.hpp"
#include "lees_edwards/LeesEdwardsBC.hpp"

#include <utils/Vector.hpp>
#include <utils/attributes.hpp>

#include <bitset>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace detail {
/**
 * @brief Get the minimum-image distance between two coordinates.
 *
 * Branchless fold: the periodicity is encoded in the masked inverse box
 * length (0 for non-periodic directions, where <tt>rint</tt> then yields a
 * zero image shift). Uses <tt>rint</tt> (round half to even) rather than
 * <tt>round</tt>, so it maps to a single rounding instruction; the two only
 * differ for separations of exactly half a box length, where both images
 * are equidistant.
 *
 * @param a                      Coordinate of the terminal point.
 * @param b                      Coordinate of the initial point.
 * @param box_length             Box length.
 * @param box_length_inv_masked  Inverse box length if periodic, 0 otherwise.
 * @return Shortest distance from @p b to @p a across periodic images,
 *         i.e. <tt>a - b</tt>. Can be negative.
 */
template <typename T>
T get_mi_coord_masked(T a, T b, T box_length, T box_length_inv_masked) {
  auto const dx = a - b;
  return dx - std::rint(dx * box_length_inv_masked) * box_length;
}

/**
 * @brief Get the minimum-image distance between two coordinates.
 * @param a           Coordinate of the terminal point.
 * @param b           Coordinate of the initial point.
 * @param box_length  Box length.
 * @param periodic    Box periodicity.
 * @return Shortest distance from @p b to @p a across periodic images,
 *         i.e. <tt>a - b</tt>. Can be negative.
 */
template <typename T> T get_mi_coord(T a, T b, T box_length, bool periodic) {
  return get_mi_coord_masked(a, b, box_length,
                             periodic ? T{1.} / box_length : T{0.});
}

/** @brief Calculate image box shift vector.
 *  @param image_box  image box offset
 *  @param box        box length
 *  @return Image box coordinates.
 */
inline auto image_shift(Utils::Vector3i const &image_box,
                        Utils::Vector3d const &box) {
  return hadamard_product(image_box, box);
}

/** @brief Unfold particle coordinates to image box.
 *  @param pos        coordinate to unfold
 *  @param image_box  image box offset
 *  @param box        box length
 *  @return Unfolded coordinates.
 */
inline auto unfolded_position(Utils::Vector3d const &pos,
                              Utils::Vector3i const &image_box,
                              Utils::Vector3d const &box) {
  return pos + image_shift(image_box, box);
}
} // namespace detail

enum class BoxType { CUBOID = 0, LEES_EDWARDS = 1 };

/**
 * @brief Cuboid minimum-image fold parameters for hot pair loops.
 *
 * Capture an instance by value in a kernel to hoist the box data out of the
 * pair loop: member loads then come from the kernel's own frame and the
 * compiler can keep them in registers, instead of re-reading them through a
 * @ref BoxGeometry reference for every pair. Only valid for cuboid boxes;
 * Lees-Edwards boxes need the full @ref BoxGeometry::get_mi_vector.
 */
class CuboidMinimumImage {
  Utils::Vector3d m_length;
  Utils::Vector3d m_length_inv_masked;

public:
  CuboidMinimumImage(Utils::Vector3d const &length,
                     Utils::Vector3d const &length_inv_masked)
      : m_length(length), m_length_inv_masked(length_inv_masked) {}

  /** @brief Squared minimum-image distance between two coordinates. */
  ESPRESSO_ATTR_ALWAYS_INLINE inline double
  dist2(Utils::Vector3d const &a, Utils::Vector3d const &b) const {
    double acc = 0.;
    for (auto c = 0u; c < 3u; ++c) {
      auto const dx = detail::get_mi_coord_masked(a[c], b[c], m_length[c],
                                                  m_length_inv_masked[c]);
      acc += dx * dx;
    }
    return acc;
  }

  /**
   * @brief Batched minimum-image vector and squared distance: one point
   * (@p xi, @p yi, @p zi) against @p m others held in the SoA arrays
   * @p sx / @p sy / @p sz.
   *
   * Writes the fold vector components to @p dx0 / @p dx1 / @p dx2 and the
   * squared distance to @p dsq. The loop carries no dependency across the
   * @p m entries and reads contiguous arrays, so it vectorizes. Each entry is
   * computed as three per-axis @c detail::get_mi_coord_masked folds followed
   * by `Utils::Vector::norm2` (accumulating from zero in component order), so
   * the per-pair results are bitwise-identical to the scalar
   * @ref BoxGeometry::get_mi_vector path for cuboid boxes.
   */
  ESPRESSO_ATTR_ALWAYS_INLINE inline void
  batch_vector_dist2(double xi, double yi, double zi, int m, double const *sx,
                     double const *sy, double const *sz, double *dx0,
                     double *dx1, double *dx2, double *dsq) const {
    auto const lx = m_length[0u];
    auto const ly = m_length[1u];
    auto const lz = m_length[2u];
    auto const ix = m_length_inv_masked[0u];
    auto const iy = m_length_inv_masked[1u];
    auto const iz = m_length_inv_masked[2u];
    for (int t = 0; t < m; ++t) {
      auto const a0 = detail::get_mi_coord_masked(xi, sx[t], lx, ix);
      auto const a1 = detail::get_mi_coord_masked(yi, sy[t], ly, iy);
      auto const a2 = detail::get_mi_coord_masked(zi, sz[t], lz, iz);
      dx0[t] = a0;
      dx1[t] = a1;
      dx2[t] = a2;
      auto acc = 0.;
      acc += a0 * a0;
      acc += a1 * a1;
      acc += a2 * a2;
      dsq[t] = acc;
    }
  }
};

class BoxGeometry {
public:
  BoxGeometry() {
    set_length(Utils::Vector3d{1., 1., 1.});
    set_periodic(0u, true);
    set_periodic(1u, true);
    set_periodic(2u, true);
    set_type(BoxType::CUBOID);
  }
  BoxGeometry(BoxGeometry const &rhs) {
    m_type = rhs.type();
    set_length(rhs.length());
    set_periodic(0u, rhs.periodic(0u));
    set_periodic(1u, rhs.periodic(1u));
    set_periodic(2u, rhs.periodic(2u));
    m_lees_edwards_bc = rhs.m_lees_edwards_bc;
  }

private:
  BoxType m_type = BoxType::CUBOID;
  /** Flags for all three dimensions whether pbc are applied (default). */
  std::bitset<3> m_periodic = 0b111;
  /** Side lengths of the box */
  Utils::Vector3d m_length = {1., 1., 1.};
  /** Inverse side lengths of the box */
  Utils::Vector3d m_length_inv = {1., 1., 1.};
  /** Inverse side lengths for periodic directions, 0 for non-periodic ones.
   *  Folding the periodicity into the inverse length makes the cuboid
   *  minimum-image fold branchless (see `detail::get_mi_coord_masked`). */
  Utils::Vector3d m_length_inv_masked = {1., 1., 1.};
  /** Half side lengths of the box */
  Utils::Vector3d m_length_half = {0.5, 0.5, 0.5};

  /** Lees-Edwards boundary conditions */
  LeesEdwardsBC m_lees_edwards_bc;

public:
  /**
   * @brief Set periodicity for direction
   *
   * @param coord The coordinate to set the periodicity for.
   * @param val True if this direction should be periodic.
   */
  void set_periodic(unsigned coord, bool val) {
    m_periodic.set(coord, val);
    m_length_inv_masked[coord] = val ? m_length_inv[coord] : 0.;
  }

  /**
   * @brief Check periodicity in direction.
   *
   * @param coord Direction to check
   * @return true iff periodic in direction.
   */
  constexpr bool periodic(unsigned coord) const {
    assert(coord <= 2u);
    return m_periodic[coord];
  }

  /**
   * @brief Box length
   * @return Return vector of side-lengths of the box.
   */
  Utils::Vector3d const &length() const { return m_length; }

  /**
   * @brief Inverse box length
   * @return Return vector of inverse side-lengths of the box.
   */
  Utils::Vector3d const &length_inv() const { return m_length_inv; }

  /**
   * @brief Half box length
   * @return Return vector of half side-lengths of the box.
   */
  Utils::Vector3d const &length_half() const { return m_length_half; }

  /**
   * @brief Set box side lengths.
   * @param box_l Length that should be set.
   */
  void set_length(Utils::Vector3d const &box_l) {
    assert(box_l > Utils::Vector3d::broadcast(0.));
    m_length = box_l;
    m_length_inv = {1. / box_l[0], 1. / box_l[1], 1. / box_l[2]};
    for (auto c = 0u; c < 3u; ++c) {
      m_length_inv_masked[c] = m_periodic[c] ? m_length_inv[c] : 0.;
    }
    m_length_half = 0.5 * box_l;
  }

  /**
   * @brief Box volume
   * @return Return the volume of the box.
   */
  double volume() const { return Utils::product(m_length); }

  /**
   * @brief Get the minimum-image distance between two coordinates.
   * @param a     Coordinate of the terminal point.
   * @param b     Coordinate of the initial point.
   * @param coord Direction
   * @return Shortest distance from @p b to @p a across periodic images,
   *         i.e. <tt>a - b</tt>. Can be negative.
   */
  template <typename T> T inline get_mi_coord(T a, T b, unsigned coord) const {
    assert(coord <= 2u);

    return detail::get_mi_coord_masked(
        a, b, static_cast<T>(m_length[coord]),
        static_cast<T>(m_length_inv_masked[coord]));
  }

  /** @brief Cuboid minimum-image fold parameters for hoisting into kernels. */
  auto cuboid_minimum_image() const {
    return CuboidMinimumImage{m_length, m_length_inv_masked};
  }

  /**
   * @brief Get the minimum-image vector between two coordinates.
   *
   * @tparam T Floating point type.
   *
   * @param a     Coordinate of the terminal point.
   * @param b     Coordinate of the initial point.
   * @return Vector from @p b to @p a that minimizes the distance across
   *         periodic images, i.e. <tt>a - b</tt>.
   */
  template <typename T>
  ESPRESSO_ATTR_ALWAYS_INLINE inline Utils::Vector3<T>
  get_mi_vector(Utils::Vector3<T> const &a, Utils::Vector3<T> const &b) const {
    if (type() == BoxType::LEES_EDWARDS) {
      auto const shear_plane_normal = lees_edwards_bc().shear_plane_normal;
      auto a_tmp = a;
      auto b_tmp = b;
      a_tmp[shear_plane_normal] = Algorithm::periodic_fold(
          a_tmp[shear_plane_normal], m_length[shear_plane_normal]);
      b_tmp[shear_plane_normal] = Algorithm::periodic_fold(
          b_tmp[shear_plane_normal], m_length[shear_plane_normal]);
      return lees_edwards_bc().distance(a_tmp - b_tmp, m_length, m_length_half,
                                        m_length_inv, m_periodic);
    }
    assert(type() == BoxType::CUBOID);
    return {get_mi_coord(a[0], b[0], 0u), get_mi_coord(a[1], b[1], 1u),
            get_mi_coord(a[2], b[2], 2u)};
  }

  /**
   * @brief Get the squared minimum-image distance between two coordinates.
   *
   * Equivalent to <tt>get_mi_vector(a, b).norm2()</tt>, but for cuboid
   * boxes avoids constructing the intermediate vector.
   *
   * @tparam T Floating point type.
   *
   * @param a Coordinate of the terminal point.
   * @param b Coordinate of the initial point.
   * @return Squared shortest distance from @p b to @p a across periodic
   *         images.
   */
  template <typename T>
  ESPRESSO_ATTR_ALWAYS_INLINE inline T
  get_mi_dist2(Utils::Vector3<T> const &a, Utils::Vector3<T> const &b) const {
    if (type() == BoxType::LEES_EDWARDS) {
      return get_mi_vector(a, b).norm2();
    }
    assert(type() == BoxType::CUBOID);
    auto const d0 = get_mi_coord(a[0], b[0], 0u);
    auto const d1 = get_mi_coord(a[1], b[1], 1u);
    auto const d2 = get_mi_coord(a[2], b[2], 2u);
    return d0 * d0 + d1 * d1 + d2 * d2;
  }

  /**
   * @brief Get the minimum-image vector between two coordinates.
   *
   * @tparam T Floating point type.
   *
   * @param a0     x element of the terminal point.
   * @param a1     y element of the terminal point.
   * @param a2     z element of the terminal point.
   * @param b0     x element of the initial point.
   * @param b1     y element of the initial point.
   * @param b2     z element of the initial point.
   * @return Vector from @p b to @p a that minimizes the distance across
   *         periodic images, i.e. <tt>a - b</tt>.
   */
  template <typename T>
  ESPRESSO_ATTR_ALWAYS_INLINE inline Utils::Vector3<T>
  get_mi_vector(T const &a0, T const &a1, T const &a2, T const &b0, T const &b1,
                T const &b2) const {
    if (type() == BoxType::LEES_EDWARDS) {
      auto const shear_plane_normal = lees_edwards_bc().shear_plane_normal;
      auto a_tmp = Utils::Vector3<T>{a0, a1, a2};
      auto b_tmp = Utils::Vector3<T>{b0, b1, b2};
      a_tmp[shear_plane_normal] = Algorithm::periodic_fold(
          a_tmp[shear_plane_normal], m_length[shear_plane_normal]);
      b_tmp[shear_plane_normal] = Algorithm::periodic_fold(
          b_tmp[shear_plane_normal], m_length[shear_plane_normal]);
      return lees_edwards_bc().distance(a_tmp - b_tmp, m_length, m_length_half,
                                        m_length_inv, m_periodic);
    }
    assert(type() == BoxType::CUBOID);
    return {get_mi_coord(a0, b0, 0u), get_mi_coord(a1, b1, 1u),
            get_mi_coord(a2, b2, 2u)};
  }

  BoxType type() const { return m_type; }
  void set_type(BoxType type) { m_type = type; }

  LeesEdwardsBC const &lees_edwards_bc() const { return m_lees_edwards_bc; }
  void set_lees_edwards_bc(LeesEdwardsBC bc) { m_lees_edwards_bc = bc; }

  /**
   * @brief Update the Lees-Edwards parameters of the box geometry
   * for the current simulation time.
   */
  void lees_edwards_update(double pos_offset, double shear_velocity) {
    assert(type() == BoxType::LEES_EDWARDS);
    m_lees_edwards_bc.pos_offset = pos_offset;
    m_lees_edwards_bc.shear_velocity = shear_velocity;
  }

  /** Calculate the velocity difference including the Lees-Edwards velocity */
  Utils::Vector3d velocity_difference(Utils::Vector3d const &x,
                                      Utils::Vector3d const &y,
                                      Utils::Vector3d const &u,
                                      Utils::Vector3d const &v) const {
    auto ret = u - v;
    if (type() == BoxType::LEES_EDWARDS) {
      auto const &le = m_lees_edwards_bc;
      auto const shear_plane_normal = le.shear_plane_normal;
      auto const shear_direction = le.shear_direction;
      auto const dy = x[shear_plane_normal] - y[shear_plane_normal];
      if (std::fabs(dy) > length_half()[shear_plane_normal]) {
        ret[shear_direction] -= std::copysign(1.0, dy) * le.shear_velocity;
      }
    }
    return ret;
  }

  /** @brief Fold coordinates to primary simulation box in-place.
   *  Lees-Edwards offset is ignored.
   *  @param[in,out] pos        coordinates to fold
   *  @param[in,out] image_box  image box offset
   */
  void fold_position(Utils::Vector3d &pos, Utils::Vector3i &image_box) const {
    for (auto i = 0u; i < 3u; i++) {
      if (m_periodic[i]) {
        auto const result =
            Algorithm::periodic_fold(pos[i], image_box[i], m_length[i]);
        if (result.second == std::numeric_limits<int>::min() or
            result.second == std::numeric_limits<int>::max()) {
          throw std::runtime_error(
              "Overflow in the image box count while folding a particle "
              "coordinate into the primary simulation box. Maybe a particle "
              "experienced a huge force.");
        }
        std::tie(pos[i], image_box[i]) = result;
      }
    }
  }

  /**
   * @brief Calculate coordinates folded to primary simulation box.
   * @param[in] pos    coordinates to fold
   * @return Folded coordinates.
   */
  auto folded_position(Utils::Vector3d const &pos) const {
    auto pos_folded = pos;
    for (auto i = 0u; i < 3u; i++) {
      if (m_periodic[i]) {
        pos_folded[i] = Algorithm::periodic_fold(pos[i], m_length[i]);
      }
    }

    return pos_folded;
  }

  /**
   * @brief Calculate image box of coordinates folded to primary simulation box.
   * @param[in] pos        coordinates
   * @param[in] image_box  image box to fold
   * @return Folded image box.
   */
  auto folded_image_box(Utils::Vector3d const &pos,
                        Utils::Vector3i const &image_box) const {
    auto image_box_folded = image_box;
    for (auto i = 0u; i < 3u; i++) {
      if (m_periodic[i]) {
        image_box_folded[i] =
            Algorithm::periodic_fold(pos[i], image_box[i], m_length[i]).second;
      }
    }

    return image_box_folded;
  }

  /** @brief Calculate image box shift vector */
  auto image_shift(Utils::Vector3i const &image_box) const {
    return detail::image_shift(image_box, m_length);
  }

  /** @brief Unfold particle coordinates to image box. */
  auto unfolded_position(Utils::Vector3d const &pos,
                         Utils::Vector3i const &image_box) const {
    return detail::unfolded_position(pos, image_box, m_length);
  }
};

/*
 * Copyright (C) 2010-2022 The ESPResSo project
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
#include <utils/math/sgn.hpp>

#include <bitset>
#include <cassert>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace detail {
/**
 * @brief Get the minimum-image distance between two coordinates.
 * @param a               Coordinate of the terminal point.
 * @param b               Coordinate of the initial point.
 * @param box_length      Box length.
 * @param box_length_inv  Inverse box length
 * @param box_length_half Half box length
 * @param periodic        Box periodicity.
 * @return Shortest distance from @p b to @p a across periodic images,
 *         i.e. <tt>a - b</tt>. Can be negative.
 */
template <typename T>
T get_mi_coord(T a, T b, T box_length, T box_length_inv, T box_length_half,
               bool periodic) {
  auto const dx = a - b;

  if (periodic && (std::abs(dx) > box_length_half)) {
    return dx - std::round(dx * box_length_inv) * box_length;
  }

  return dx;
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
  return get_mi_coord(a, b, box_length, 1. / box_length, 0.5 * box_length,
                      periodic);
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
  void set_periodic(unsigned coord, bool val) { m_periodic.set(coord, val); }

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
    m_length = box_l;
    m_length_inv = {1. / box_l[0], 1. / box_l[1], 1. / box_l[2]};
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

    return detail::get_mi_coord(a, b, m_length[coord], m_length_inv[coord],
                                m_length_half[coord], m_periodic[coord]);
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
  Utils::Vector<T, 3> get_mi_vector(const Utils::Vector<T, 3> &a,
                                    const Utils::Vector<T, 3> &b) const {
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
    return {get_mi_coord(a[0], b[0], 0), get_mi_coord(a[1], b[1], 1),
            get_mi_coord(a[2], b[2], 2)};
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
      if (fabs(dy) > 0.5 * length_half()[shear_plane_normal]) {
        ret[shear_direction] -= Utils::sgn(dy) * le.shear_velocity;
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
    for (unsigned int i = 0u; i < 3u; i++) {
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

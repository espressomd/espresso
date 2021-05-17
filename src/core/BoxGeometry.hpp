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
#ifndef CORE_BOX_GEOMETRY_HPP
#define CORE_BOX_GEOMETRY_HPP

#include "algorithm/periodic_fold.hpp"

#include <utils/Vector.hpp>

#include <bitset>
#include <cassert>
#include <cmath>

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
} // namespace detail

class BoxGeometry {
private:
  /** Flags for all three dimensions whether pbc are applied (default). */
  std::bitset<3> m_periodic = 0b111;
  /** Side lengths of the box */
  Utils::Vector3d m_length = {1, 1, 1};
  /** Inverse side lengths of the box */
  Utils::Vector3d m_length_inv = {1, 1, 1};
  /** Half side lengths of the box */
  Utils::Vector3d m_length_half = {0.5, 0.5, 0.5};

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
    assert(coord <= 2);
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
  double volume() const { return m_length[0] * m_length[1] * m_length[2]; }

  /**
   * @brief Get the minimum-image distance between two coordinates.
   * @param a     Coordinate of the terminal point.
   * @param b     Coordinate of the initial point.
   * @param coord Direction
   * @return Shortest distance from @p b to @p a across periodic images,
   *         i.e. <tt>a - b</tt>. Can be negative.
   */
  template <typename T> T inline get_mi_coord(T a, T b, unsigned coord) const {
    assert(coord <= 2);

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
    return {get_mi_coord(a[0], b[0], 0), get_mi_coord(a[1], b[1], 1),
            get_mi_coord(a[2], b[2], 2)};
  }
};

/** @brief Fold a coordinate to primary simulation box.
 *  @param pos        coordinate to fold
 *  @param image_box  image box offset
 *  @param length     box length
 */
inline std::pair<double, int> fold_coordinate(double pos, int image_box,
                                              double const &length) {
  std::tie(pos, image_box) = Algorithm::periodic_fold(pos, image_box, length);

  if ((image_box == std::numeric_limits<int>::min()) ||
      (image_box == std::numeric_limits<int>::max())) {
    throw std::runtime_error(
        "Overflow in the image box count while folding a particle coordinate "
        "into the primary simulation box. Maybe a particle experienced a "
        "huge force.");
  }

  return {pos, image_box};
}

/** @brief Fold particle coordinates to primary simulation box.
 *  @param[in,out] pos        coordinate to fold
 *  @param[in,out] image_box  image box offset
 *  @param[in] box            box parameters (side lengths, periodicity)
 */
inline void fold_position(Utils::Vector3d &pos, Utils::Vector3i &image_box,
                          const BoxGeometry &box) {
  for (int i = 0; i < 3; i++) {
    if (box.periodic(i)) {
      std::tie(pos[i], image_box[i]) =
          fold_coordinate(pos[i], image_box[i], box.length()[i]);
    }
  }
}

/** @brief Fold particle coordinates to primary simulation box.
 *  @param p    coordinate to fold
 *  @param box  box parameters (side lengths, periodicity)
 *  @return Folded coordinates.
 */
inline Utils::Vector3d folded_position(const Utils::Vector3d &p,
                                       const BoxGeometry &box) {
  Utils::Vector3d p_folded;
  for (int i = 0; i < 3; i++) {
    if (box.periodic(i)) {
      p_folded[i] = Algorithm::periodic_fold(p[i], box.length()[i]);
    } else {
      p_folded[i] = p[i];
    }
  }

  return p_folded;
}

#endif // CORE_BOX_GEOMETRY_HPP

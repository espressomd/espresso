/*
 * Copyright (C) 2026 The ESPResSo project
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

#include "Shape.hpp"

#include <utils/Vector.hpp>

#include <array>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <utility>
#include <vector>

namespace Shapes {

/**
 * @brief Solid/pore mask on a regular voxel grid, stored per MPI rank.
 *
 * Only the voxels of the local region plus a halo are stored, one bit per
 * voxel. The unfolded voxel index @c j along an axis refers to the voxel
 * center @c origin + (j + 0.5) * voxel_size. Voxels outside of the mask
 * wrap periodically by the box length, or are pore.
 */
class VoxelMask : public Shape {
public:
  struct Geometry {
    /** Dataset shape, in voxels. */
    Utils::Vector3i mask_shape{};
    double voxel_size = 1.;
    /** Position of the lower corner of the mask. */
    Utils::Vector3d origin{};
    Utils::Vector3d box_l{};
    std::array<bool, 3> periodic{{true, true, true}};
    /** Halo width, in voxels (at least 1). */
    int halo = 1;
  };

  /** Range of unfolded voxel indices: [lower, lower + extent). */
  struct Region {
    Utils::Vector3i lower{};
    Utils::Vector3i extent{};
  };

  /** @brief Voxels to store for a rank that owns [local_lower, local_upper). */
  static Region stored_region(Geometry const &g,
                              Utils::Vector3d const &local_lower,
                              Utils::Vector3d const &local_upper);

  /**
   * @brief Map an unfolded voxel index to a mask index along one axis.
   * @return Mask index, or -1 if the voxel is pore by definition.
   */
  static int mask_index(Geometry const &g, int axis, int j);

  /** @brief Set the geometry. Discards the local data. */
  void set_geometry(Geometry const &g);
  Geometry const &geometry() const { return m_geometry; }
  Region const &region() const { return m_region; }

  /**
   * @brief Set the voxels of a region.
   * @param region Unfolded index range covered by @p voxels.
   * @param voxels Row-major over the extent, x slowest; nonzero is solid.
   */
  void set_local_data(Region const &region,
                      std::vector<std::uint8_t> const &voxels);

  /** @brief Set the region, with all voxels pore. */
  void allocate_local_data(Region const &region);

  /**
   * @brief Set the voxels of a box inside the region.
   * @param lower  First voxel of the box, relative to the region.
   * @param extent Extent of the box.
   * @param voxels Row-major over @p extent, x slowest; nonzero is solid.
   */
  void set_local_block(Utils::Vector3i const &lower,
                       Utils::Vector3i const &extent,
                       std::vector<std::uint8_t> const &voxels);

  /** @brief Set the callback that (re)loads the local data. */
  void set_loader(std::function<void(VoxelMask &)> loader) {
    m_loader = std::move(loader);
  }
  /** @brief Call the loader, if set. */
  void on_local_geometry_change() override;

  /**
   * @brief Distance reported if no interface is within the halo.
   * It lies beyond any interaction cutoff, and its square is finite.
   */
  static constexpr double far_distance = 1e100;

  /**
   * @brief Signed distance to the nearest interface within the halo.
   * Positive in pore, negative in solid. @c vec is the position minus the
   * closest point on the interface. Distances are exact up to
   * <tt>halo * voxel_size</tt> only. If no interface is within the halo,
   * the result is only a "far" flag: the magnitude is @ref far_distance
   * along the x axis.
   * @throw std::runtime_error if @p pos is outside the local region.
   */
  void calculate_dist(Utils::Vector3d const &pos, double &dist,
                      Utils::Vector3d &vec) const override;
  bool is_inside(Utils::Vector3d const &pos) const override;
  /** @brief <tt>halo * voxel_size</tt>. */
  double local_distance_reach() const override {
    return m_geometry.halo * m_geometry.voxel_size;
  }

  std::size_t local_voxel_count() const;
  std::size_t local_memory_bytes() const;

private:
  /** Shift @p pos by box vectors into the stored region. */
  Utils::Vector3d image_in_region(Utils::Vector3d const &pos) const;
  /** Solid state of an unfolded voxel inside the stored region. */
  bool solid(Utils::Vector3i const &index) const;

  Geometry m_geometry{};
  Region m_region{};
  std::vector<std::uint64_t> m_bits;
  std::function<void(VoxelMask &)> m_loader;
};

} /* namespace Shapes */

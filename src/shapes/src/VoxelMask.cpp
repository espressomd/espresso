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

#include <shapes/VoxelMask.hpp>

#include <utils/Vector.hpp>

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <vector>

namespace Shapes {

VoxelMask::Region VoxelMask::stored_region(Geometry const &g,
                                           Utils::Vector3d const &local_lower,
                                           Utils::Vector3d const &local_upper) {
  Region region;
  for (int a = 0; a < 3; ++a) {
    auto const lower = static_cast<int>(std::floor(
                           (local_lower[a] - g.origin[a]) / g.voxel_size)) -
                       g.halo;
    auto const upper = static_cast<int>(std::ceil(
                           (local_upper[a] - g.origin[a]) / g.voxel_size)) +
                       g.halo;
    region.lower[a] = lower;
    region.extent[a] = upper - lower;
  }
  return region;
}

int VoxelMask::mask_index(Geometry const &g, int axis, int j) {
  auto const h = g.voxel_size;
  auto c = g.origin[axis] + (j + 0.5) * h;
  if (g.periodic[axis]) {
    c -= std::floor(c / g.box_l[axis]) * g.box_l[axis];
  } else if (c < 0. or c >= g.box_l[axis]) {
    return -1;
  }
  auto const m = static_cast<int>(std::floor((c - g.origin[axis]) / h));
  return (m >= 0 and m < g.mask_shape[axis]) ? m : -1;
}

void VoxelMask::set_geometry(Geometry const &g) {
  m_geometry = g;
  m_region = {};
  m_bits.clear();
  m_bits.shrink_to_fit();
}

static std::size_t voxel_count(VoxelMask::Region const &region) {
  return static_cast<std::size_t>(region.extent[0]) *
         static_cast<std::size_t>(region.extent[1]) *
         static_cast<std::size_t>(region.extent[2]);
}

void VoxelMask::set_local_data(Region const &region,
                               std::vector<std::uint8_t> const &voxels) {
  if (voxels.size() != voxel_count(region)) {
    throw std::runtime_error("VoxelMask: data size does not match the region");
  }
  allocate_local_data(region);
  set_local_block({0, 0, 0}, region.extent, voxels);
}

void VoxelMask::allocate_local_data(Region const &region) {
  m_region = region;
  m_bits.assign((voxel_count(region) + 63) / 64, std::uint64_t{0});
  m_bits.shrink_to_fit();
}

void VoxelMask::set_local_block(Utils::Vector3i const &lower,
                                Utils::Vector3i const &extent,
                                std::vector<std::uint8_t> const &voxels) {
  auto const &e = m_region.extent;
  for (int a = 0; a < 3; ++a) {
    if (lower[a] < 0 or extent[a] < 0 or lower[a] + extent[a] > e[a]) {
      throw std::runtime_error("VoxelMask: block lies outside of the region");
    }
  }
  if (voxels.size() != voxel_count({lower, extent})) {
    throw std::runtime_error("VoxelMask: data size does not match the block");
  }
  std::size_t n = 0;
  for (int i = 0; i < extent[0]; ++i) {
    for (int j = 0; j < extent[1]; ++j) {
      auto bit = (static_cast<std::size_t>(lower[0] + i) *
                      static_cast<std::size_t>(e[1]) +
                  static_cast<std::size_t>(lower[1] + j)) *
                     static_cast<std::size_t>(e[2]) +
                 static_cast<std::size_t>(lower[2]);
      for (int k = 0; k < extent[2]; ++k, ++n, ++bit) {
        auto const mask = std::uint64_t{1} << (bit % 64);
        if (voxels[n]) {
          m_bits[bit / 64] |= mask;
        } else {
          m_bits[bit / 64] &= ~mask;
        }
      }
    }
  }
}

void VoxelMask::on_local_geometry_change() {
  if (m_loader) {
    m_loader(*this);
  }
}

std::size_t VoxelMask::local_voxel_count() const {
  return voxel_count(m_region);
}

std::size_t VoxelMask::local_memory_bytes() const {
  return m_bits.capacity() * sizeof(std::uint64_t);
}

bool VoxelMask::solid(Utils::Vector3i const &index) const {
  auto const i = static_cast<std::size_t>(index[0] - m_region.lower[0]);
  auto const j = static_cast<std::size_t>(index[1] - m_region.lower[1]);
  auto const k = static_cast<std::size_t>(index[2] - m_region.lower[2]);
  auto const n = (i * m_region.extent[1] + j) * m_region.extent[2] + k;
  return ((m_bits[n / 64] >> (n % 64)) & 1u) != 0u;
}

static Utils::Vector3i voxel_of(VoxelMask::Geometry const &g,
                                Utils::Vector3d const &pos) {
  Utils::Vector3i index;
  for (int a = 0; a < 3; ++a) {
    index[a] =
        static_cast<int>(std::floor((pos[a] - g.origin[a]) / g.voxel_size));
  }
  return index;
}

Utils::Vector3d VoxelMask::image_in_region(Utils::Vector3d const &pos) const {
  if (m_bits.empty()) {
    throw std::runtime_error("VoxelMask: no data loaded on this rank");
  }
  auto const &g = m_geometry;
  Utils::Vector3d image = pos;
  for (int a = 0; a < 3; ++a) {
    bool found = false;
    for (int n : {0, -1, 1}) {
      if (n != 0 and not g.periodic[a]) {
        break;
      }
      auto const x = pos[a] + n * g.box_l[a];
      auto const j =
          static_cast<int>(std::floor((x - g.origin[a]) / g.voxel_size));
      if (j >= m_region.lower[a] and
          j < m_region.lower[a] + m_region.extent[a]) {
        image[a] = x;
        found = true;
        break;
      }
    }
    if (not found) {
      throw std::runtime_error(
          "VoxelMask: position outside the local region of this rank");
    }
  }
  return image;
}

bool VoxelMask::is_inside(Utils::Vector3d const &pos) const {
  return solid(voxel_of(m_geometry, image_in_region(pos)));
}

void VoxelMask::calculate_dist(Utils::Vector3d const &pos, double &dist,
                               Utils::Vector3d &vec) const {
  auto const &g = m_geometry;
  auto const h = g.voxel_size;
  auto const p = image_in_region(pos);
  auto const center = voxel_of(g, p);
  auto const query_solid = solid(center);

  // inclusive search range
  Utils::Vector3i lo, hi;
  for (int a = 0; a < 3; ++a) {
    lo[a] = std::max(center[a] - g.halo, m_region.lower[a]);
    hi[a] = std::min(center[a] + g.halo,
                     m_region.lower[a] + m_region.extent[a] - 1);
  }

  auto best = std::numeric_limits<double>::max();
  Utils::Vector3d best_vec{};
  Utils::Vector3i index;
  for (index[0] = lo[0]; index[0] <= hi[0]; ++index[0]) {
    for (index[1] = lo[1]; index[1] <= hi[1]; ++index[1]) {
      for (index[2] = lo[2]; index[2] <= hi[2]; ++index[2]) {
        // look for solid cubes from a pore, and for pore cubes from a solid
        if (solid(index) == query_solid) {
          continue;
        }
        Utils::Vector3d v;
        for (int a = 0; a < 3; ++a) {
          auto const cube_lower = g.origin[a] + index[a] * h;
          v[a] = p[a] - std::clamp(p[a], cube_lower, cube_lower + h);
        }
        auto const d2 = v.norm2();
        if (d2 < best) {
          best = d2;
          best_vec = v;
        }
      }
    }
  }

  if (best == std::numeric_limits<double>::max()) {
    dist = far_distance;
    if (query_solid) {
      dist = -dist;
    }
    vec = dist * Utils::Vector3d{1., 0., 0.};
    return;
  }
  dist = query_solid ? -std::sqrt(best) : std::sqrt(best);
  vec = best_vec;
}

} /* namespace Shapes */

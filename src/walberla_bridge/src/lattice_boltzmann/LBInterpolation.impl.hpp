/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

/**
 * @file
 * Out-of-class position-based interpolation definitions for
 * @ref walberla::LBWalberlaImpl.
 */

#include <utils/Vector.hpp>
#include <utils/interpolation/bspline_3d.hpp>

#include <algorithm>
#include <array>
#include <cstddef>
#include <iostream>
#include <iterator>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace walberla {

/**
 * @brief Exception for accessing a lattice node outside the local domain
 *  and ghost layers during B-spline interpolation.
 */
class interpolation_illegal_access : public std::runtime_error {
public:
  interpolation_illegal_access(std::string const &field,
                               Utils::Vector3d const &pos,
                               std::array<int, 3> const &node, double weight)
      : std::runtime_error("Access to LB " + field + " field failed") {
    std::cerr << "pos [" << pos << "], node [" << Utils::Vector3i(node)
              << "], weight " << weight << "\n";
  }
};

void interpolate_bspline_at_pos(Utils::Vector3d const &pos, auto const &&f) {
  Utils::Interpolation::bspline_3d<2>(
      pos, f, Utils::Vector3d::broadcast(1.), // grid spacing
      Utils::Vector3d::broadcast(.5));        // offset
}

template <typename FloatType, lbmpy::Arch Architecture>
std::function<bool(Utils::Vector3d const &)>
LBWalberlaImpl<FloatType, Architecture>::make_lattice_position_checker(
    bool consider_points_in_halo) const {
  auto const &lat = *m_lattice;
  if (consider_points_in_halo) {
    return [&](Utils::Vector3d const &p) { return lat.pos_in_local_halo(p); };
  }
  return [&](Utils::Vector3d const &p) { return lat.pos_in_local_domain(p); };
}

/**
 * @brief Distribute forces to the lattice at given positions.
 * Uses B-spline interpolation to spread each force over the surrounding
 * lattice nodes. On GPU, positions are transformed to block-local
 * coordinates and the operation is performed in a single kernel launch.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::add_forces_at_pos(
    std::vector<Utils::Vector3d> const &pos,
    std::vector<Utils::Vector3d> const &forces) {
  assert(pos.size() == forces.size());
  if (pos.empty()) {
    return;
  }
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_force_interpolation_kernel();
    for (std::size_t i = 0ul; i < pos.size(); ++i) {
      kernel(pos[i], forces[i]);
    }
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    auto const &lattice = get_lattice();
    auto const &block = *(lattice.get_blocks()->begin());
    auto const origin = block.getAABB().min();
    std::vector<FloatType> host_pos;
    std::vector<FloatType> host_force;
    host_pos.reserve(3ul * pos.size());
    host_force.reserve(3ul * forces.size());
    assert(lattice.get_blocks()->getNumberOfBlocks() == 1u);
    for (auto const &vec : pos) {
#pragma unroll
      for (std::size_t i : {0ul, 1ul, 2ul}) {
        host_pos.emplace_back(static_cast<FloatType>(vec[i] - origin[i]));
      }
    }
    for (auto const &vec : forces) {
#pragma unroll
      for (std::size_t i : {0ul, 1ul, 2ul}) {
        host_force.emplace_back(static_cast<FloatType>(vec[i]));
      }
    }
    zero_centered_to_lb_in_place(host_force);
    auto const gl = lattice.get_ghost_layers();
    auto field = block.template uncheckedFastGetData<VectorField>(
        m_force_to_be_applied_id);
    lbm::accessor::Interpolation::add_force(field, host_pos, host_force, gl);
  }
#endif
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImpl<FloatType, Architecture>::make_force_interpolation_kernel()
    const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos, Utils::Vector3d const &force) {
    if (not get_block_extended(lattice, pos, 1u)) {
      return;
    }
    interpolate_bspline_at_pos(
        pos, [&, conv = m_zc_to_lb, field_id = m_force_to_be_applied_id](
                 std::array<int, 3> const node, double weight) {
          auto block = get_block_extended(lattice, node, 0u);
          if (!block)
            block = get_block_extended(lattice, node, 1u);
          if (block) {
            auto cell = to_cell(node);
            blocks.transformGlobalToBlockLocalCell(cell, *block);
            weight *= conv;
            auto const weighted_force = to_vector3<FloatType>(weight * force);
            auto field =
                block->template uncheckedFastGetData<VectorField>(field_id);
            lbm::accessor::Vector::add(field, weighted_force, cell);
          }
        });
  };
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImpl<FloatType,
                    Architecture>::make_velocity_interpolation_kernel() const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos) {
    Utils::Vector3d acc{0., 0., 0.};
    interpolate_bspline_at_pos(
        pos, [&, field_id = m_velocity_field_id](std::array<int, 3> const node,
                                                 double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0.) {
            auto block = get_block_extended(lattice, node, 1u);
            if (!block)
              throw interpolation_illegal_access("velocity", pos, node, weight);
            Vector3<FloatType> vel;
            if (m_has_boundaries and m_boundary->node_is_boundary(node)) {
              vel = m_boundary->get_node_value_at_boundary(node);
            } else {
              auto cell = to_cell(node);
              blocks.transformGlobalToBlockLocalCell(cell, *block);
              auto field =
                  block->template uncheckedFastGetData<VectorField>(field_id);
              vel = lbm::accessor::Vector::get(field, cell);
            }
            acc += to_vector3d(vel) * weight;
          }
        });
    return acc;
  };
}

template <typename FloatType, lbmpy::Arch Architecture>
auto LBWalberlaImpl<FloatType,
                    Architecture>::make_density_interpolation_kernel() const {
  auto const &lattice = *m_lattice;
  auto const &blocks = *lattice.get_blocks();
  assert(lattice.get_ghost_layers() == 1u);
  return [&](Utils::Vector3d const &pos) {
    double acc = 0.;
    interpolate_bspline_at_pos(
        pos, [&, density = m_density, field_id = m_pdf_field_id](
                 std::array<int, 3> const node, double weight) {
          // Nodes with zero weight might not be accessible, because they can be
          // outside ghost layers
          if (weight != 0.) {
            auto block = get_block_extended(lattice, node, 1u);
            if (!block)
              throw interpolation_illegal_access("density", pos, node, weight);
            auto cell = to_cell(node);
            blocks.transformGlobalToBlockLocalCell(cell, *block);
            auto field =
                block->template uncheckedFastGetData<PdfField>(field_id);
            auto const rho = lbm::accessor::Density::get(field, density, cell);
            acc += rho * weight;
          }
        });
    return acc;
  };
}

/**
 * @brief Interpolate velocities at given positions (batch version).
 * On GPU, boundary slip velocities are written into the velocity field
 * before interpolation, since the field has indeterminate values inside
 * boundary regions.
 */
template <typename FloatType, lbmpy::Arch Architecture>
std::vector<Utils::Vector3d>
LBWalberlaImpl<FloatType, Architecture>::get_velocities_at_pos(
    std::vector<Utils::Vector3d> const &pos) {
  if (pos.empty()) {
    return {};
  }
  std::vector<Utils::Vector3d> vel{};
  vel.reserve(pos.size());
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_velocity_interpolation_kernel();
    std::ranges::transform(pos, std::back_inserter(vel), kernel);
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    auto const &lattice = get_lattice();
    auto const &block = *(lattice.get_blocks()->begin());
    auto const origin = block.getAABB().min();
    std::vector<FloatType> host_pos;
    host_pos.reserve(3ul * pos.size());
    assert(lattice.get_blocks()->getNumberOfBlocks() == 1u);
    for (auto const &vec : pos) {
#pragma unroll
      for (std::size_t i : {0ul, 1ul, 2ul}) {
        host_pos.emplace_back(static_cast<FloatType>(vec[i] - origin[i]));
      }
    }
    auto const gl = lattice.get_ghost_layers();
    auto field =
        block.template uncheckedFastGetData<VectorField>(m_velocity_field_id);
    // the velocity field has indeterminate values inside boundary regions;
    // we overwrite them with boundary slip velocities before interpolation
    auto const [dev_idx, dev_vel] = m_boundary->get_flattened_map_device();
    if (not dev_idx->empty()) {
      lbm::accessor::Vector::set_from_list(field, *dev_idx, *dev_vel, gl);
    }
    auto const res = lbm::accessor::Interpolation::get_vel(field, host_pos, gl);
    for (auto it = res.begin(); it != res.end(); it += 3) {
      vel.emplace_back(Utils::Vector3d{static_cast<double>(*(it + 0)),
                                       static_cast<double>(*(it + 1)),
                                       static_cast<double>(*(it + 2))});
    }
  }
#endif
  return vel;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImpl<FloatType, Architecture>::get_densities_at_pos(
    std::vector<Utils::Vector3d> const &pos) {
  if (pos.empty()) {
    return {};
  }
  std::vector<double> rho{};
  rho.reserve(pos.size());
  if constexpr (Architecture == lbmpy::Arch::CPU) {
    auto const kernel = make_density_interpolation_kernel();
    std::ranges::transform(pos, std::back_inserter(rho), kernel);
  }
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    auto const &lattice = get_lattice();
    auto const &block = *(lattice.get_blocks()->begin());
    auto const origin = block.getAABB().min();
    std::vector<FloatType> host_pos;
    host_pos.reserve(3ul * pos.size());
    assert(lattice.get_blocks()->getNumberOfBlocks() == 1u);
    for (auto const &vec : pos) {
#pragma unroll
      for (std::size_t i : {0ul, 1ul, 2ul}) {
        host_pos.emplace_back(static_cast<FloatType>(vec[i] - origin[i]));
      }
    }
    auto const gl = lattice.get_ghost_layers();
    auto field = block.template uncheckedFastGetData<PdfField>(m_pdf_field_id);
    auto res =
        lbm::accessor::Interpolation::get_rho(field, host_pos, m_density, gl);
    if constexpr (std::is_same_v<FloatType, double>) {
      std::swap(rho, res);
    } else {
      for (auto const &v : res) {
        rho.emplace_back(static_cast<double>(v));
      }
    }
  }
#endif
  return rho;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<Utils::Vector3d>
LBWalberlaImpl<FloatType, Architecture>::get_velocity_at_pos(
    Utils::Vector3d const &pos, bool consider_points_in_halo) const {
  assert(not m_pending_ghost_comm.test(GhostComm::VEL));
  assert(not m_pending_ghost_comm.test(GhostComm::UBB));
  if (!consider_points_in_halo and !m_lattice->pos_in_local_domain(pos))
    return std::nullopt;
  if (consider_points_in_halo and !m_lattice->pos_in_local_halo(pos))
    return std::nullopt;
  auto const kernel = make_velocity_interpolation_kernel();
  return {kernel(pos)};
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<double>
LBWalberlaImpl<FloatType, Architecture>::get_density_at_pos(
    Utils::Vector3d const &pos, bool consider_points_in_halo) const {
  assert(not m_pending_ghost_comm.test(GhostComm::PDF));
  if (!consider_points_in_halo and !m_lattice->pos_in_local_domain(pos))
    return std::nullopt;
  if (consider_points_in_halo and !m_lattice->pos_in_local_halo(pos))
    return std::nullopt;
  auto const kernel = make_density_interpolation_kernel();
  return {kernel(pos)};
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImpl<FloatType, Architecture>::add_force_at_pos(
    Utils::Vector3d const &pos, Utils::Vector3d const &force) {
  if (!m_lattice->pos_in_local_halo(pos))
    return false;
  auto const kernel = make_force_interpolation_kernel();
  kernel(pos, force);
  return true;
}

} // namespace walberla

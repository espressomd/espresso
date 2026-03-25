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
 * Out-of-class collision model setup definitions for
 * @ref walberla::LBWalberlaImpl.
 */

#include <memory>
#include <stdexcept>
#include <utility>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
FloatType
LBWalberlaImpl<FloatType, Architecture>::shear_mode_relaxation_rate() const {
  return FloatType{2} / (FloatType{6} * m_viscosity + FloatType{1});
}

template <typename FloatType, lbmpy::Arch Architecture>
FloatType LBWalberlaImpl<FloatType, Architecture>::odd_mode_relaxation_rate(
    FloatType shear_relaxation, FloatType magic_number) const {
  return (FloatType{4} - FloatType{2} * shear_relaxation) /
         (FloatType{4} * magic_number * shear_relaxation + FloatType{2} -
          shear_relaxation);
}

/**
 * @brief Set up the thermalized collision model.
 * Configures MRT relaxation rates from the viscosity and initializes the
 * random number generator for fluctuating LB.
 * When @p kT is zero, fluctuations are suppressed but the thermalized
 * kernel is still used (required for RNG state tracking).
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_collision_model(
    double kT, unsigned int seed) {
  auto const omega = shear_mode_relaxation_rate();
  auto const omega_odd = odd_mode_relaxation_rate(omega);
  auto const blocks = get_lattice().get_blocks();
  m_kT = FloatType_c(kT);
  m_seed = seed;
  auto obj = typename Kernels::StreamCollisionModelThermalized(
      m_last_applied_force_field_id, m_pdf_field_id, zero_centered_to_lb(m_kT),
      omega, omega, omega_odd, omega, seed, uint32_t{0u});
  m_collision_model = std::make_shared<CollisionModel>(std::move(obj));
  m_run_stream_collide_sweep = StreamCollideSweepVisitor(blocks);
  setup_streaming_communicator();
}

/**
 * @brief Set up the Lees-Edwards collision model.
 * Configures a modified collision step that applies the shear velocity
 * at the Lees-Edwards boundary planes, and sets up the interpolation
 * sweeps for PDFs, velocities, and forces at those boundaries.
 * Currently restricted to shear_plane_normal="y" and no domain
 * decomposition along the shear direction.
 */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_collision_model(
    std::unique_ptr<LeesEdwardsPack> &&lees_edwards_pack) {
  assert(m_kT == 0.);
#if defined(__CUDACC__) and defined(WALBERLA_BUILD_WITH_CUDA)
  if constexpr (Architecture == lbmpy::Arch::GPU) {
    throw std::runtime_error("Lees-Edwards LB doesn't support GPU yet");
  }
#endif
  auto const shear_direction = lees_edwards_pack->shear_direction;
  auto const shear_plane_normal = lees_edwards_pack->shear_plane_normal;
  auto const shear_vel = FloatType_c(lees_edwards_pack->get_shear_velocity());
  auto const omega = shear_mode_relaxation_rate();
  auto const omega_odd = odd_mode_relaxation_rate(omega);
  if (shear_plane_normal != 1u) {
    throw std::domain_error(
        "Lees-Edwards LB only supports shear_plane_normal=\"y\"");
  }
  auto const &lattice = get_lattice();
  auto const n_ghost_layers = lattice.get_ghost_layers();
  auto const blocks = lattice.get_blocks();
  if (lattice.get_node_grid()[shear_direction] != 1 or
      blocks->getSize(shear_direction) != 1ul or
      blocks->getSize(shear_plane_normal) !=
          lattice.get_node_grid()[shear_plane_normal]) {
    throw std::domain_error("LB LEbc doesn't support domain decomposition "
                            "along the shear direction, nor multiple blocks "
                            "along the normal direction");
  }
  auto const &grid_dimensions = lattice.get_grid_dimensions();
  auto const block_origin = lattice.get_local_grid_range(false).first;
  auto const lebc_slab_origin = block_origin[shear_plane_normal];
  auto const lebc_slab_total_thickness = grid_dimensions[shear_plane_normal];
  auto const lebc_bot_index = 0 - lebc_slab_origin;
  auto const lebc_top_index = lebc_slab_total_thickness - lebc_slab_origin;
  m_collision_model = std::make_shared<CollisionModel>(
      typename Kernels::StreamCollisionModelLeesEdwards(
          m_last_applied_force_field_id, m_pdf_field_id, lebc_bot_index,
          lebc_top_index, omega, omega, omega_odd, omega, shear_vel));
  m_lees_edwards_callbacks = std::move(lees_edwards_pack);
  m_run_stream_collide_sweep =
      StreamCollideSweepVisitor(blocks, m_lees_edwards_callbacks);
  m_lees_edwards_pdf_interpol_sweep =
      std::make_shared<InterpolateAndShiftAtBoundary<_PdfField, FloatType>>(
          blocks, m_pdf_field_id, m_pdf_tmp_field_id, n_ghost_layers,
          shear_direction, shear_plane_normal,
          m_lees_edwards_callbacks->get_pos_offset);
  m_lees_edwards_vel_interpol_sweep =
      std::make_shared<InterpolateAndShiftAtBoundary<_VectorField, FloatType>>(
          blocks, m_velocity_field_id, m_vel_tmp_field_id, n_ghost_layers,
          shear_direction, shear_plane_normal,
          m_lees_edwards_callbacks->get_pos_offset,
          m_lees_edwards_callbacks->get_shear_velocity);
  m_lees_edwards_last_applied_force_interpol_sweep =
      std::make_shared<InterpolateAndShiftAtBoundary<_VectorField, FloatType>>(
          blocks, m_last_applied_force_field_id, m_vel_tmp_field_id,
          n_ghost_layers, shear_direction, shear_plane_normal,
          m_lees_edwards_callbacks->get_pos_offset);
  setup_streaming_communicator();
}

/** @brief Verify that MD and LB Lees-Edwards parameters are consistent. */
template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::check_lebc(
    unsigned int shear_direction, unsigned int shear_plane_normal) const {
  if (m_lees_edwards_callbacks) {
    if (m_lees_edwards_callbacks->shear_direction != shear_direction or
        m_lees_edwards_callbacks->shear_plane_normal != shear_plane_normal) {
      throw std::runtime_error(
          "MD and LB Lees-Edwards boundary conditions disagree");
    }
  }
}

} // namespace walberla

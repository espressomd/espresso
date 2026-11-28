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
 * Out-of-class node access definitions for
 * @ref walberla::LBWalberlaImpl.
 */

#include <utils/Vector.hpp>

#include <array>
#include <optional>
#include <utility>
#include <vector>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<Utils::Vector3d>
LBWalberlaImpl<FloatType, Architecture>::get_node_velocity(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::VEL)));
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::UBB)));
  if (m_has_boundaries) {
    auto const is_boundary = get_node_is_boundary(node, consider_ghosts);
    if (is_boundary and *is_boundary) {
      return get_node_velocity_at_boundary(node, consider_ghosts);
    }
  }
  auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto field = bc->block->template uncheckedFastGetData<VectorField>(
      m_velocity_field_id);
  auto const vec = lbm::accessor::Vector::get(field, bc->cell);
  return to_vector3d(vec);
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImpl<FloatType, Architecture>::set_node_velocity(
    Utils::Vector3i const &node, Utils::Vector3d const &v) {
  m_pending_ghost_comm.set(GhostComm::PDF);
  m_pending_ghost_comm.set(GhostComm::VEL);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  // We have to set both, the pdf and the stored velocity field
  auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
  auto vel_field =
      bc->block->template getData<VectorField>(m_velocity_field_id);
  auto force_field =
      bc->block->template getData<VectorField>(m_last_applied_force_field_id);
  auto vel = to_vector3<FloatType>(v);
  lbm::accessor::Velocity::set(pdf_field, vel_field, force_field, vel,
                               bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<double> LBWalberlaImpl<FloatType, Architecture>::get_node_density(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::PDF)));
  auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto pdf_field =
      bc->block->template uncheckedFastGetData<PdfField>(m_pdf_field_id);
  auto const density =
      lbm::accessor::Density::get(pdf_field, m_density, bc->cell);
  return {double_c(density)};
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImpl<FloatType, Architecture>::set_node_density(
    Utils::Vector3i const &node, double density) {
  m_pending_ghost_comm.set(GhostComm::PDF);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
  lbm::accessor::Density::set(pdf_field, FloatType_c(density), m_density,
                              bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<std::vector<double>>
LBWalberlaImpl<FloatType, Architecture>::get_node_population(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::PDF)));
  auto bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
  auto const pop = lbm::accessor::Population::get(pdf_field, bc->cell);
  std::vector<double> population(Stencil::Size);
  for (uint_t f = 0u; f < Stencil::Size; ++f) {
    population[f] = double_c(pop[f]);
  }

  return {std::move(population)};
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImpl<FloatType, Architecture>::set_node_population(
    Utils::Vector3i const &node, std::vector<double> const &population) {
  m_pending_ghost_comm.set(GhostComm::PDF);
  m_pending_ghost_comm.set(GhostComm::VEL);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
  auto force_field =
      bc->block->template getData<VectorField>(m_last_applied_force_field_id);
  auto vel_field =
      bc->block->template getData<VectorField>(m_velocity_field_id);
  std::array<FloatType, Stencil::Size> pop;
  for (uint_t f = 0u; f < Stencil::Size; ++f) {
    pop[f] = FloatType_c(population[f]);
  }
  lbm::accessor::Population::set(pdf_field, vel_field, force_field, pop,
                                 bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<Utils::Vector3d>
LBWalberlaImpl<FloatType, Architecture>::get_node_force_to_be_applied(
    Utils::Vector3i const &node) const {
  auto const bc = get_block_and_cell(get_lattice(), node, true);
  if (!bc)
    return std::nullopt;

  auto field =
      bc->block->template getData<VectorField>(m_force_to_be_applied_id);
  auto const vec = lbm::accessor::Vector::get(field, bc->cell);
  return zero_centered_to_md(to_vector3d(vec));
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<Utils::Vector3d>
LBWalberlaImpl<FloatType, Architecture>::get_node_last_applied_force(
    Utils::Vector3i const &node, bool consider_ghosts) const {
  assert(not(consider_ghosts and m_pending_ghost_comm.test(GhostComm::LAF)));
  auto const bc = get_block_and_cell(get_lattice(), node, consider_ghosts);
  if (!bc)
    return std::nullopt;

  auto const field =
      bc->block->template getData<VectorField>(m_last_applied_force_field_id);
  auto const vec = lbm::accessor::Vector::get(field, bc->cell);
  return zero_centered_to_md(to_vector3d(vec));
}

template <typename FloatType, lbmpy::Arch Architecture>
bool LBWalberlaImpl<FloatType, Architecture>::set_node_last_applied_force(
    Utils::Vector3i const &node, Utils::Vector3d const &force) {
  m_pending_ghost_comm.set(GhostComm::VEL);
  m_pending_ghost_comm.set(GhostComm::LAF);
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return false;

  auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
  auto force_field =
      bc->block->template getData<VectorField>(m_last_applied_force_field_id);
  auto vel_field =
      bc->block->template getData<VectorField>(m_velocity_field_id);
  auto const vec = to_vector3<FloatType>(force);
  lbm::accessor::Force::set(pdf_field, vel_field, force_field, vec, m_density,
                            bc->cell);

  return true;
}

template <typename FloatType, lbmpy::Arch Architecture>
std::optional<Utils::VectorXd<9>>
LBWalberlaImpl<FloatType, Architecture>::get_node_pressure_tensor(
    Utils::Vector3i const &node) const {
  auto bc = get_block_and_cell(get_lattice(), node, false);
  if (!bc)
    return std::nullopt;

  auto pdf_field = bc->block->template getData<PdfField>(m_pdf_field_id);
  auto tensor =
      lbm::accessor::PressureTensor::get(pdf_field, m_density, bc->cell);
  pressure_tensor_correction(tensor);
  return to_vector9d(tensor);
}

} // namespace walberla

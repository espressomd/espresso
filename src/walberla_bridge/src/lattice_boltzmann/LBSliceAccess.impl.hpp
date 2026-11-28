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
 * Out-of-class slice access definitions for
 * @ref walberla::LBWalberlaImpl.
 */

#include <utils/Vector.hpp>

#include <span>
#include <vector>

namespace walberla {

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double> LBWalberlaImpl<FloatType, Architecture>::get_slice_velocity(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(3u * ci.numCells());
        auto const field =
            block.template getData<VectorField>(m_velocity_field_id);
        auto values = lbm::accessor::Vector::get(field, bci);

        auto kernel = [&values, &out, this](unsigned const block_index,
                                            unsigned const local_index,
                                            Utils::Vector3i const &node) {
          if (m_boundary->node_is_boundary(node)) {
            auto const &vec = m_boundary->get_node_value_at_boundary(node);
            for (uint_t f = 0u; f < 3u; ++f) {
              out[3u * local_index + f] = vec[f];
            }
          } else {
            for (uint_t f = 0u; f < 3u; ++f) {
              out[3u * local_index + f] =
                  double_c(values[3u * block_index + f]);
            }
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_slice_velocity(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &velocity) {
  m_pending_ghost_comm.set(GhostComm::PDF);
  m_pending_ghost_comm.set(GhostComm::VEL);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(velocity.size() == 3u * ci.numCells());
        auto pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        auto force_field =
            block.template getData<VectorField>(m_last_applied_force_field_id);
        auto vel_field =
            block.template getData<VectorField>(m_velocity_field_id);
        std::vector<FloatType> values(3u * bci.numCells());

        auto kernel = [&values, &velocity](unsigned const block_index,
                                           unsigned const local_index,
                                           Utils::Vector3i const &) {
          for (uint_t f = 0u; f < 3u; ++f) {
            values[3u * block_index + f] =
                numeric_cast<FloatType>(velocity[3u * local_index + f]);
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        lbm::accessor::Velocity::set(pdf_field, vel_field, force_field, values,
                                     bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImpl<FloatType, Architecture>::get_slice_last_applied_force(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(3u * ci.numCells());
        auto const field =
            block.template getData<VectorField>(m_last_applied_force_field_id);
        auto const values = lbm::accessor::Vector::get(field, bci);

        auto kernel = [&values, &out](unsigned const block_index,
                                      unsigned const local_index,
                                      Utils::Vector3i const &) {
          for (uint_t f = 0u; f < 3u; ++f) {
            out[3u * local_index + f] = values[3u * block_index + f];
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  zero_centered_to_md_in_place(out);
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_slice_last_applied_force(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &force) {
  m_pending_ghost_comm.set(GhostComm::VEL);
  m_pending_ghost_comm.set(GhostComm::LAF);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(force.size() == 3u * ci.numCells());
        auto pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        auto force_field =
            block.template getData<VectorField>(m_last_applied_force_field_id);
        auto vel_field =
            block.template getData<VectorField>(m_velocity_field_id);
        std::vector<FloatType> values(3u * bci.numCells());

        auto kernel = [&values, &force](unsigned const block_index,
                                        unsigned const local_index,
                                        Utils::Vector3i const &) {
          for (uint_t f = 0u; f < 3u; ++f) {
            values[3u * block_index + f] =
                numeric_cast<FloatType>(force[3u * local_index + f]);
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        lbm::accessor::Force::set(pdf_field, vel_field, force_field, values,
                                  m_density, bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImpl<FloatType, Architecture>::get_slice_population(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(stencil_size() * ci.numCells());
        auto const pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        auto const values = lbm::accessor::Population::get(pdf_field, bci);

        auto kernel = [&values, &out, this](unsigned const block_index,
                                            unsigned const local_index,
                                            Utils::Vector3i const &) {
          for (uint_t f = 0u; f < stencil_size(); ++f) {
            out[stencil_size() * local_index + f] =
                values[stencil_size() * block_index + f];
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_slice_population(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &population) {
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(population.size() == stencil_size() * ci.numCells());
        auto pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        auto force_field =
            block.template getData<VectorField>(m_last_applied_force_field_id);
        auto vel_field =
            block.template getData<VectorField>(m_velocity_field_id);
        std::vector<FloatType> values(stencil_size() * bci.numCells());

        auto kernel = [&values, &population, this](unsigned const block_index,
                                                   unsigned const local_index,
                                                   Utils::Vector3i const &) {
          for (uint_t f = 0u; f < stencil_size(); ++f) {
            values[stencil_size() * block_index + f] = numeric_cast<FloatType>(
                population[stencil_size() * local_index + f]);
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        lbm::accessor::Population::set(pdf_field, vel_field, force_field,
                                       values, bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double> LBWalberlaImpl<FloatType, Architecture>::get_slice_density(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(ci.numCells());
        auto const pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        auto const values =
            lbm::accessor::Density::get(pdf_field, m_density, bci);

        auto kernel = [&values, &out](unsigned const block_index,
                                      unsigned const local_index,
                                      Utils::Vector3i const &) {
          out[local_index] = values[block_index];
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

template <typename FloatType, lbmpy::Arch Architecture>
void LBWalberlaImpl<FloatType, Architecture>::set_slice_density(
    Utils::Vector3i const &lower_corner, Utils::Vector3i const &upper_corner,
    std::vector<double> const &density) {
  m_pending_ghost_comm.set(GhostComm::PDF);
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        assert(density.size() == ci.numCells());
        auto pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        std::vector<FloatType> values(bci.numCells());

        auto kernel = [&values, &density](unsigned const block_index,
                                          unsigned const local_index,
                                          Utils::Vector3i const &) {
          values[block_index] = numeric_cast<FloatType>(density[local_index]);
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
        lbm::accessor::Density::set(pdf_field, values, m_density, bci);
      });
}

template <typename FloatType, lbmpy::Arch Architecture>
std::vector<double>
LBWalberlaImpl<FloatType, Architecture>::get_slice_pressure_tensor(
    Utils::Vector3i const &lower_corner,
    Utils::Vector3i const &upper_corner) const {
  std::vector<double> out;
  for_each_block_in_slice(
      get_lattice(), lower_corner, upper_corner,
      [&](auto const &block, auto const &bci, auto const &ci,
          auto const &block_offset) {
        if (out.empty())
          out.resize(9u * ci.numCells());
        auto const pdf_field = block.template getData<PdfField>(m_pdf_field_id);
        auto values =
            lbm::accessor::PressureTensor::get(pdf_field, m_density, bci);

        auto kernel = [&values, &out, this](unsigned const block_index,
                                            unsigned const local_index,
                                            Utils::Vector3i const &) {
          pressure_tensor_correction(
              std::span<FloatType, 9ul>(&values[9u * block_index], 9ul));
          for (uint_t f = 0u; f < 9u; ++f) {
            out[9u * local_index + f] = values[9u * block_index + f];
          }
        };

        copy_block_buffer(bci, ci, block_offset, lower_corner, kernel);
      });
  return out;
}

} // namespace walberla

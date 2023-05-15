/*
 * Copyright (C) 2022-2023 The ESPResSo project
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

#include "PoissonSolver.hpp"

#include <blockforest/communication/UniformBufferedScheme.h>
#include <domain_decomposition/BlockDataID.h>
#include <fft/Fft.h>
#include <field/AddToStorage.h>
#include <field/GhostLayerField.h>
#include <field/communication/PackInfo.h>
#include <stencil/D3Q27.h>

#include <utils/constants.hpp>

#include <cmath>
#include <cstddef>
#include <memory>
#include <utility>

namespace walberla {

template <typename FloatType> class FFT : public PoissonSolver {
private:
  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  domain_decomposition::BlockDataID m_potential_field_id;

  using PotentialField = GhostLayerField<FloatType, 1>;

  std::shared_ptr<fft::FourierTransform<PotentialField>> m_ft;
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

  using FullCommunicator = blockforest::communication::UniformBufferedScheme<
      typename stencil::D3Q27>;
  std::shared_ptr<FullCommunicator> m_full_communication;

public:
  FFT(std::shared_ptr<LatticeWalberla> lattice, double permittivity)
      : PoissonSolver(std::move(lattice), permittivity) {
    m_blocks = get_lattice().get_blocks();

    Vector3<uint_t> dim(m_blocks->getNumberOfXCells(),
                        m_blocks->getNumberOfYCells(),
                        m_blocks->getNumberOfZCells());
    auto const greens = [dim](uint_t x, uint_t y, uint_t z) -> real_t {
      if (x == 0u && y == 0u && z == 0u)
        return 0.;
      return -0.5 /
             (std::cos(2. * Utils::pi() * real_c(x) / real_c(dim[0])) +
              std::cos(2. * Utils::pi() * real_c(y) / real_c(dim[1])) +
              std::cos(2. * Utils::pi() * real_c(z) / real_c(dim[2])) - 3.) /
             real_c(dim[0] * dim[1] * dim[2]);
    };

    m_potential_field_id = field::addToStorage<PotentialField>(
        get_lattice().get_blocks(), "potential field", 0.0, field::fzyx,
        get_lattice().get_ghost_layers());

    m_ft = std::make_shared<fft::FourierTransform<PotentialField>>(
        m_blocks, m_potential_field_id, greens);

    m_full_communication =
        std::make_shared<FullCommunicator>(get_lattice().get_blocks());
    m_full_communication->addPackInfo(
        std::make_shared<field::communication::PackInfo<PotentialField>>(
            m_potential_field_id));
  }
  ~FFT() override = default;

  void reset_charge_field() override {
    // the FFT-solver re-uses the potential field for the charge
    auto const potential_id = walberla::BlockDataID(get_potential_field_id());

    for (auto &block : *get_lattice().get_blocks()) {
      auto field = block.template getData<PotentialField>(potential_id);
      WALBERLA_FOR_ALL_CELLS_XYZ(field, field->get(x, y, z) = 0.;)
    }
  }

  void add_charge_to_field(std::size_t id, double valency,
                           bool is_double_precision) override {
    auto const factor = FloatType_c(valency) / FloatType_c(get_permittivity());
    // the FFT-solver re-uses the potential field for the charge
    const auto charge_id = walberla::BlockDataID(get_potential_field_id());
    const auto density_id = walberla::BlockDataID(id);
    for (auto &block : *get_lattice().get_blocks()) {
      auto charge_field = block.template getData<PotentialField>(charge_id);
      if (is_double_precision) {
        auto density_field =
            block.template getData<walberla::GhostLayerField<double, 1>>(
                density_id);
        WALBERLA_FOR_ALL_CELLS_XYZ(
            charge_field, charge_field->get(x, y, z) +=
                          factor * FloatType_c(density_field->get(x, y, z));)
      } else {
        auto density_field =
            block.template getData<walberla::GhostLayerField<float, 1>>(
                density_id);
        WALBERLA_FOR_ALL_CELLS_XYZ(
            charge_field, charge_field->get(x, y, z) +=
                          factor * FloatType_c(density_field->get(x, y, z));)
      }
    }
  }

  [[nodiscard]] std::size_t get_potential_field_id() const noexcept override {
    return static_cast<std::size_t>(m_potential_field_id);
  }

  void solve() override {
    (*m_ft)();
    ghost_communication();
  }

private:
  void ghost_communication() { (*m_full_communication)(); }
};

} // namespace walberla

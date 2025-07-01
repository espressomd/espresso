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
#include <FFT_CUDA.cuh>

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <utility>

namespace walberla {

template <typename FloatType> class FFT_GPU : public PoissonSolver {
private:
  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  std::shared_ptr<FFT_CUDA<FloatType>> fft_cuda;


public:
  FFT_GPU() = default;
  FFT_GPU(std::shared_ptr<LatticeWalberla> lattice, double permittivity)
      : PoissonSolver(lattice, permittivity) {
    fft_cuda = std::make_shared<FFT_CUDA<FloatType>>(lattice, permittivity);
  }
  ~FFT_GPU() override = default;

  void reset_charge_field() override {
    fft_cuda->reset_charge_field();
  }

  void add_charge_to_field(std::size_t id, double valency,
                           bool is_double_precision) override {
    fft_cuda->add_charge_to_field(id, valency, is_double_precision);

  }

  std::size_t get_potential_field_id() const noexcept override {
    return fft_cuda->get_potential_field_id();
  }

  void solve() {
    fft_cuda->solve();
  }

  void set_permittivity(double permittivity) noexcept override{
    fft_cuda->set_permittivity(permittivity);
  }

  [[nodiscard]] double get_permittivity() const noexcept override{
    return fft_cuda->get_permittivity();
  }

  [[nodiscard]] auto const &get_lattice() const noexcept{
    return fft_cuda->get_lattice();
  }

private:
  void ghost_communication() {
    fft_cuda->ghost_communication();
  }
};

} // namespace walberla

/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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

#include <config/config.hpp>

#ifdef ESPRESSO_P3M

#include "electrostatics/p3m_heffte.impl.hpp"

#include <memory>
#include <utility>

template <typename FloatType, Arch Architecture>
std::shared_ptr<CoulombP3M>
new_coulomb_p3m_impl(P3MParameters &&p3m, TuningParameters const &tuning_params,
                     double prefactor) {
  auto constexpr r_space_order = Utils::MemoryOrder::ROW_MAJOR;
  auto constexpr k_space_order = Utils::MemoryOrder::ROW_MAJOR;
  using FFTConfig = P3MFFTConfig<r_space_order, k_space_order, true, 2>;
  auto state_ptr =
      std::make_unique<CoulombP3MState<FloatType, Architecture, FFTConfig>>(
          std::move(p3m));
  auto obj =
      std::make_shared<CoulombP3MHeffte<FloatType, Architecture, FFTConfig>>(
          std::move(state_ptr), tuning_params, prefactor);
  return obj;
}

std::shared_ptr<CoulombP3M>
new_coulomb_p3m_heffte(P3MParameters &&p3m_params,
                       TuningParameters const &tuning_params, double prefactor,
                       bool single_precision, Arch arch) {
  auto fptr = &new_coulomb_p3m_impl<float, Arch::CPU>;
  if (single_precision) {
    if (arch == Arch::CPU) {
      fptr = new_coulomb_p3m_impl<float, Arch::CPU>;
    } else {
      fptr = new_coulomb_p3m_impl<float, Arch::CUDA>;
    }
  } else {
    if (arch == Arch::CPU) {
      fptr = new_coulomb_p3m_impl<double, Arch::CPU>;
    } else {
      throw std::invalid_argument(
          "P3M GPU only implemented in single-precision mode");
    }
  }
  return fptr(std::move(p3m_params), tuning_params, prefactor);
}
#endif // ESPRESSO_P3M

/*
 * Copyright (C) 2010-2025 The ESPResSo project
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

#ifdef ESPRESSO_DP3M

#include "magnetostatics/dp3m_heffte.impl.hpp"

#include "p3m/FFTBackendLegacy.hpp"
#include "p3m/FFTBuffersLegacy.hpp"

#include <memory>
#include <utility>

template <typename FloatType, Arch Architecture>
std::shared_ptr<DipolarP3M>
new_dipolar_p3m_impl(P3MParameters &&p3m, TuningParameters const &tuning_params,
                     double prefactor) {
  using DefaultFFTConfig =
      P3MFFTConfig<Utils::MemoryOrder::ROW_MAJOR, Utils::MemoryOrder::ROW_MAJOR,
                   false, 2>;
  auto obj = std::make_shared<
      DipolarP3MHeffte<FloatType, Architecture, DefaultFFTConfig>>(
      std::make_unique<DipolarP3MState<FloatType, DefaultFFTConfig>>(
          std::move(p3m)),
      tuning_params, prefactor);
  obj->dp3m.template make_mesh_instance<FFTBuffersLegacy<FloatType>>();
  obj->dp3m.template make_fft_instance<FFTBackendLegacy<FloatType>>();
  return obj;
}

std::shared_ptr<DipolarP3M>
new_dipolar_p3m_heffte(P3MParameters &&p3m_params,
                       TuningParameters const &tuning_params, double prefactor,
                       bool single_precision, Arch arch) {
  auto fptr = &new_dipolar_p3m_impl<float, Arch::CPU>;
  if (single_precision) {
    fptr = new_dipolar_p3m_impl<float, Arch::CPU>;
  } else {
    fptr = new_dipolar_p3m_impl<double, Arch::CPU>;
  }
  return fptr(std::move(p3m_params), tuning_params, prefactor);
}

#endif // ESPRESSO_DP3M

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

#include "generated_kernels/AdvectiveFluxKernel_double_precision.h"
#include "generated_kernels/AdvectiveFluxKernel_single_precision.h"
#include "generated_kernels/ContinuityKernel_double_precision.h"
#include "generated_kernels/ContinuityKernel_single_precision.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostatic_double_precision.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostatic_single_precision.h"
#include "generated_kernels/DiffusiveFluxKernel_double_precision.h"
#include "generated_kernels/DiffusiveFluxKernel_single_precision.h"
#include "generated_kernels/FrictionCouplingKernel_double_precision.h"
#include "generated_kernels/FrictionCouplingKernel_single_precision.h"

#include "generated_kernels/Dirichlet_double_precision.h"
#include "generated_kernels/Dirichlet_single_precision.h"
#include "generated_kernels/FixedFlux_double_precision.h"
#include "generated_kernels/FixedFlux_single_precision.h"

namespace walberla {
namespace detail {
template <typename FloatType = double> struct KernelTrait {
  using ContinuityKernel = pystencils::ContinuityKernel_double_precision;
  using DiffusiveFluxKernel = pystencils::DiffusiveFluxKernel_double_precision;
  using AdvectiveFluxKernel = pystencils::AdvectiveFluxKernel_double_precision;
  using FrictionCouplingKernel =
      pystencils::FrictionCouplingKernel_double_precision;
  using DiffusiveFluxKernelElectrostatic =
      pystencils::DiffusiveFluxKernelWithElectrostatic_double_precision;

  using Dirichlet = pystencils::Dirichlet_double_precision;
  using FixedFlux = pystencils::FixedFlux_double_precision;
};
template <> struct KernelTrait<float> {
  using ContinuityKernel = pystencils::ContinuityKernel_single_precision;
  using DiffusiveFluxKernel = pystencils::DiffusiveFluxKernel_single_precision;
  using AdvectiveFluxKernel = pystencils::AdvectiveFluxKernel_single_precision;
  using FrictionCouplingKernel =
      pystencils::FrictionCouplingKernel_single_precision;
  using DiffusiveFluxKernelElectrostatic =
      pystencils::DiffusiveFluxKernelWithElectrostatic_single_precision;

  using Dirichlet = pystencils::Dirichlet_single_precision;
  using FixedFlux = pystencils::FixedFlux_single_precision;
};
} // namespace detail
} // namespace walberla

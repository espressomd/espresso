/*
 * Copyright (C) 2025-2026 The ESPResSo project
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

#include "generated_kernels/AdvectiveFluxKernelDoublePrecisionCUDA.h"
#include "generated_kernels/AdvectiveFluxKernelSinglePrecisionCUDA.h"
#include "generated_kernels/ContinuityKernelDoublePrecisionCUDA.h"
#include "generated_kernels/ContinuityKernelSinglePrecisionCUDA.h"
#include "generated_kernels/DiffusiveFluxKernelDoublePrecisionCUDA.h"
#include "generated_kernels/DiffusiveFluxKernelSinglePrecisionCUDA.h"
#include "generated_kernels/DiffusiveFluxKernelThermalizedDoublePrecisionCUDA.h"
#include "generated_kernels/DiffusiveFluxKernelThermalizedSinglePrecisionCUDA.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticDoublePrecisionCUDA.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticSinglePrecisionCUDA.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticThermalizedDoublePrecisionCUDA.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticThermalizedSinglePrecisionCUDA.h"
#include "generated_kernels/DirichletDoublePrecisionCUDA.h"
#include "generated_kernels/DirichletSinglePrecisionCUDA.h"
#include "generated_kernels/EKFieldAccessorsDoublePrecisionCUDA.cuh"
#include "generated_kernels/EKFieldAccessorsSinglePrecisionCUDA.cuh"
#include "generated_kernels/FixedFluxDoublePrecisionCUDA.h"
#include "generated_kernels/FixedFluxSinglePrecisionCUDA.h"
#include "generated_kernels/FrictionCouplingKernelDoublePrecisionCUDA.h"
#include "generated_kernels/FrictionCouplingKernelSinglePrecisionCUDA.h"

namespace walberla {
namespace detail {

using lbmpy::Arch;

template <> struct KernelTrait<double, Arch::GPU> {
  using ContinuityKernel = pystencils::ContinuityKernelDoublePrecisionCUDA;
  using DiffusiveFluxKernel =
      pystencils::DiffusiveFluxKernelDoublePrecisionCUDA;
  using DiffusiveFluxKernelThermalized =
      pystencils::DiffusiveFluxKernelThermalizedDoublePrecisionCUDA;
  using AdvectiveFluxKernel =
      pystencils::AdvectiveFluxKernelDoublePrecisionCUDA;
  using FrictionCouplingKernel =
      pystencils::FrictionCouplingKernelDoublePrecisionCUDA;
  using DiffusiveFluxKernelElectrostatic =
      pystencils::DiffusiveFluxKernelWithElectrostaticDoublePrecisionCUDA;
  using DiffusiveFluxKernelElectrostaticThermalized = pystencils::
      DiffusiveFluxKernelWithElectrostaticThermalizedDoublePrecisionCUDA;

  using Dirichlet = ek::DirichletDoublePrecisionCUDA;
  using FixedFlux = ek::FixedFluxDoublePrecisionCUDA;
};
template <> struct KernelTrait<float, Arch::GPU> {
  using ContinuityKernel = pystencils::ContinuityKernelSinglePrecisionCUDA;
  using DiffusiveFluxKernel =
      pystencils::DiffusiveFluxKernelSinglePrecisionCUDA;
  using DiffusiveFluxKernelThermalized =
      pystencils::DiffusiveFluxKernelThermalizedSinglePrecisionCUDA;
  using AdvectiveFluxKernel =
      pystencils::AdvectiveFluxKernelSinglePrecisionCUDA;
  using FrictionCouplingKernel =
      pystencils::FrictionCouplingKernelSinglePrecisionCUDA;
  using DiffusiveFluxKernelElectrostatic =
      pystencils::DiffusiveFluxKernelWithElectrostaticSinglePrecisionCUDA;
  using DiffusiveFluxKernelElectrostaticThermalized = pystencils::
      DiffusiveFluxKernelWithElectrostaticThermalizedSinglePrecisionCUDA;

  using Dirichlet = ek::DirichletSinglePrecisionCUDA;
  using FixedFlux = ek::FixedFluxSinglePrecisionCUDA;
};
} // namespace detail
} // namespace walberla

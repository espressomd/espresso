/*
 * Copyright (C) 2025 The ESPResSo project
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

#include "generated_kernels/AdvectiveFluxKernel_double_precision_CUDA.h"
#include "generated_kernels/AdvectiveFluxKernel_single_precision_CUDA.h"
#include "generated_kernels/ContinuityKernel_double_precision_CUDA.h"
#include "generated_kernels/ContinuityKernel_single_precision_CUDA.h"
#include "generated_kernels/DiffusiveFluxKernelThermalized_double_precision_CUDA.h"
#include "generated_kernels/DiffusiveFluxKernelThermalized_single_precision_CUDA.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticThermalized_double_precision_CUDA.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticThermalized_single_precision_CUDA.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostatic_double_precision_CUDA.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostatic_single_precision_CUDA.h"
#include "generated_kernels/DiffusiveFluxKernel_double_precision_CUDA.h"
#include "generated_kernels/DiffusiveFluxKernel_single_precision_CUDA.h"
#include "generated_kernels/EK_FieldAccessors_double_precision_CUDA.cuh"
#include "generated_kernels/EK_FieldAccessors_single_precision_CUDA.cuh"
#include "generated_kernels/FrictionCouplingKernel_double_precision_CUDA.h"
#include "generated_kernels/FrictionCouplingKernel_single_precision_CUDA.h"

#include "generated_kernels/Dirichlet_double_precision_CUDA.h"
#include "generated_kernels/Dirichlet_single_precision_CUDA.h"
#include "generated_kernels/FixedFlux_double_precision_CUDA.h"
#include "generated_kernels/FixedFlux_single_precision_CUDA.h"

namespace walberla {
namespace detail {

using lbmpy::Arch;

template <> struct KernelTrait<double, Arch::GPU> {
  using ContinuityKernel = pystencils::ContinuityKernel_double_precision_CUDA;
  using DiffusiveFluxKernel =
      pystencils::DiffusiveFluxKernel_double_precision_CUDA;
  using DiffusiveFluxKernelThermalized =
      pystencils::DiffusiveFluxKernelThermalized_double_precision_CUDA;
  using AdvectiveFluxKernel =
      pystencils::AdvectiveFluxKernel_double_precision_CUDA;
  using FrictionCouplingKernel =
      pystencils::FrictionCouplingKernel_double_precision_CUDA;
  using DiffusiveFluxKernelElectrostatic =
      pystencils::DiffusiveFluxKernelWithElectrostatic_double_precision_CUDA;
  using DiffusiveFluxKernelElectrostaticThermalized = pystencils::
      DiffusiveFluxKernelWithElectrostaticThermalized_double_precision_CUDA;

  using Dirichlet = pystencils::Dirichlet_double_precision_CUDA;
  using FixedFlux = pystencils::FixedFlux_double_precision_CUDA;
};
template <> struct KernelTrait<float, Arch::GPU> {
  using ContinuityKernel = pystencils::ContinuityKernel_single_precision_CUDA;
  using DiffusiveFluxKernel =
      pystencils::DiffusiveFluxKernel_single_precision_CUDA;
  using DiffusiveFluxKernelThermalized =
      pystencils::DiffusiveFluxKernelThermalized_single_precision_CUDA;
  using AdvectiveFluxKernel =
      pystencils::AdvectiveFluxKernel_single_precision_CUDA;
  using FrictionCouplingKernel =
      pystencils::FrictionCouplingKernel_single_precision_CUDA;
  using DiffusiveFluxKernelElectrostatic =
      pystencils::DiffusiveFluxKernelWithElectrostatic_single_precision_CUDA;
  using DiffusiveFluxKernelElectrostaticThermalized = pystencils::
      DiffusiveFluxKernelWithElectrostaticThermalized_single_precision_CUDA;

  using Dirichlet = pystencils::Dirichlet_single_precision_CUDA;
  using FixedFlux = pystencils::FixedFlux_single_precision_CUDA;
};
} // namespace detail
} // namespace walberla

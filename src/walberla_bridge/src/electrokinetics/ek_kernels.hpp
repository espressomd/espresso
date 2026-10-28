/*
 * Copyright (C) 2022-2026 The ESPResSo project
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

#include <walberla_bridge/Architecture.hpp>

#include "generated_kernels/AdvectiveFluxKernelDoublePrecision.h"
#include "generated_kernels/AdvectiveFluxKernelSinglePrecision.h"
#include "generated_kernels/ContinuityKernelDoublePrecision.h"
#include "generated_kernels/ContinuityKernelSinglePrecision.h"
#include "generated_kernels/DiffusiveFluxKernelDoublePrecision.h"
#include "generated_kernels/DiffusiveFluxKernelSinglePrecision.h"
#include "generated_kernels/DiffusiveFluxKernelThermalizedDoublePrecision.h"
#include "generated_kernels/DiffusiveFluxKernelThermalizedSinglePrecision.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticDoublePrecision.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticSinglePrecision.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticThermalizedDoublePrecision.h"
#include "generated_kernels/DiffusiveFluxKernelWithElectrostaticThermalizedSinglePrecision.h"
#include "generated_kernels/DirichletDoublePrecision.h"
#include "generated_kernels/DirichletSinglePrecision.h"
#include "generated_kernels/EKFieldAccessorsDoublePrecision.h"
#include "generated_kernels/EKFieldAccessorsSinglePrecision.h"
#include "generated_kernels/FixedFluxDoublePrecision.h"
#include "generated_kernels/FixedFluxSinglePrecision.h"
#include "generated_kernels/FrictionCouplingKernelDoublePrecision.h"
#include "generated_kernels/FrictionCouplingKernelSinglePrecision.h"

namespace walberla {
namespace detail {

using lbmpy::Arch;

template <typename FT = double, Arch AT = Arch::CPU> struct KernelTrait {
  using ContinuityKernel = pystencils::ContinuityKernelDoublePrecision;
  using DiffusiveFluxKernel = pystencils::DiffusiveFluxKernelDoublePrecision;
  using DiffusiveFluxKernelThermalized =
      pystencils::DiffusiveFluxKernelThermalizedDoublePrecision;
  using AdvectiveFluxKernel = pystencils::AdvectiveFluxKernelDoublePrecision;
  using FrictionCouplingKernel =
      pystencils::FrictionCouplingKernelDoublePrecision;
  using DiffusiveFluxKernelElectrostatic =
      pystencils::DiffusiveFluxKernelWithElectrostaticDoublePrecision;
  using DiffusiveFluxKernelElectrostaticThermalized = pystencils::
      DiffusiveFluxKernelWithElectrostaticThermalizedDoublePrecision;

  using Dirichlet = ek::DirichletDoublePrecision;
  using FixedFlux = ek::FixedFluxDoublePrecision;
};
template <> struct KernelTrait<float, Arch::CPU> {
  using ContinuityKernel = pystencils::ContinuityKernelSinglePrecision;
  using DiffusiveFluxKernel = pystencils::DiffusiveFluxKernelSinglePrecision;
  using DiffusiveFluxKernelThermalized =
      pystencils::DiffusiveFluxKernelThermalizedSinglePrecision;
  using AdvectiveFluxKernel = pystencils::AdvectiveFluxKernelSinglePrecision;
  using FrictionCouplingKernel =
      pystencils::FrictionCouplingKernelSinglePrecision;
  using DiffusiveFluxKernelElectrostatic =
      pystencils::DiffusiveFluxKernelWithElectrostaticSinglePrecision;
  using DiffusiveFluxKernelElectrostaticThermalized = pystencils::
      DiffusiveFluxKernelWithElectrostaticThermalizedSinglePrecision;

  using Dirichlet = ek::DirichletSinglePrecision;
  using FixedFlux = ek::FixedFluxSinglePrecision;
};
} // namespace detail
} // namespace walberla

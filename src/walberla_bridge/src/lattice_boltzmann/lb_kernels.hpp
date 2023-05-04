/*
 * Copyright (C) 2021-2023 The ESPResSo project
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

#include "generated_kernels/Dynamic_UBB_double_precision.h"
#include "generated_kernels/Dynamic_UBB_single_precision.h"
#include "generated_kernels/FieldAccessorsDoublePrecision.h"
#include "generated_kernels/FieldAccessorsSinglePrecision.h"
#include "generated_kernels/InitialPDFsSetterDoublePrecision.h"
#include "generated_kernels/InitialPDFsSetterSinglePrecision.h"
#include "generated_kernels/StreamSweepDoublePrecision.h"
#include "generated_kernels/StreamSweepSinglePrecision.h"

#ifdef __AVX2__
#include "generated_kernels/CollideSweepDoublePrecisionLeesEdwardsAVX.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalizedAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionLeesEdwardsAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalizedAVX.h"
#else
#include "generated_kernels/CollideSweepDoublePrecisionLeesEdwards.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalized.h"
#include "generated_kernels/CollideSweepSinglePrecisionLeesEdwards.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalized.h"
#endif

namespace walberla {
namespace detail {

using lbmpy::Arch;

template <typename FT = double, Arch AT = Arch::CPU> struct KernelTrait {
#ifdef __AVX2__
  using CollisionModelThermalized =
      pystencils::CollideSweepDoublePrecisionThermalizedAVX;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepDoublePrecisionLeesEdwardsAVX;
#else
  using CollisionModelThermalized =
      pystencils::CollideSweepDoublePrecisionThermalized;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepDoublePrecisionLeesEdwards;
#endif
  using StreamSweep = pystencils::StreamSweepDoublePrecision;
  using InitialPDFsSetter = pystencils::InitialPDFsSetterDoublePrecision;
};

template <> struct KernelTrait<float, Arch::CPU> {
#ifdef __AVX2__
  using CollisionModelThermalized =
      pystencils::CollideSweepSinglePrecisionThermalizedAVX;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepSinglePrecisionLeesEdwardsAVX;
#else
  using CollisionModelThermalized =
      pystencils::CollideSweepSinglePrecisionThermalized;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepSinglePrecisionLeesEdwards;
#endif
  using StreamSweep = pystencils::StreamSweepSinglePrecision;
  using InitialPDFsSetter = pystencils::InitialPDFsSetterSinglePrecision;
};

template <typename FT = double, Arch AT = Arch::CPU>
struct BoundaryHandlingTrait {
  using Dynamic_UBB = lbm::Dynamic_UBB_double_precision;
};

template <> struct BoundaryHandlingTrait<float, Arch::CPU> {
  using Dynamic_UBB = lbm::Dynamic_UBB_single_precision;
};

} // namespace detail
} // namespace walberla

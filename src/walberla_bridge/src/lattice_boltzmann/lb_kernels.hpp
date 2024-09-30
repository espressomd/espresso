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
#include "generated_kernels/PackInfoPdfDoublePrecision.h"
#include "generated_kernels/PackInfoPdfSinglePrecision.h"
#include "generated_kernels/PackInfoVecDoublePrecision.h"
#include "generated_kernels/PackInfoVecSinglePrecision.h"

#ifdef __AVX2__
#include "generated_kernels/CollideSweepDoublePrecisionLeesEdwardsAVX.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalizedAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionLeesEdwardsAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalizedAVX.h"
#include "generated_kernels/StreamSweepDoublePrecisionAVX.h"
#include "generated_kernels/StreamSweepSinglePrecisionAVX.h"
#else
#include "generated_kernels/CollideSweepDoublePrecisionLeesEdwards.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalized.h"
#include "generated_kernels/CollideSweepSinglePrecisionLeesEdwards.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalized.h"
#include "generated_kernels/StreamSweepDoublePrecision.h"
#include "generated_kernels/StreamSweepSinglePrecision.h"
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
  using StreamSweep = pystencils::StreamSweepDoublePrecisionAVX;
#else
  using CollisionModelThermalized =
      pystencils::CollideSweepDoublePrecisionThermalized;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepDoublePrecisionLeesEdwards;
  using StreamSweep = pystencils::StreamSweepDoublePrecision;
#endif
  using InitialPDFsSetter = pystencils::InitialPDFsSetterDoublePrecision;
  using PackInfoPdf = pystencils::PackInfoPdfDoublePrecision;
  using PackInfoVec = pystencils::PackInfoVecDoublePrecision;
};

template <> struct KernelTrait<float, Arch::CPU> {
#ifdef __AVX2__
  using CollisionModelThermalized =
      pystencils::CollideSweepSinglePrecisionThermalizedAVX;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepSinglePrecisionLeesEdwardsAVX;
  using StreamSweep = pystencils::StreamSweepSinglePrecisionAVX;
#else
  using CollisionModelThermalized =
      pystencils::CollideSweepSinglePrecisionThermalized;
  using CollisionModelLeesEdwards =
      pystencils::CollideSweepSinglePrecisionLeesEdwards;
  using StreamSweep = pystencils::StreamSweepSinglePrecision;
#endif
  using InitialPDFsSetter = pystencils::InitialPDFsSetterSinglePrecision;
  using PackInfoPdf = pystencils::PackInfoPdfSinglePrecision;
  using PackInfoVec = pystencils::PackInfoVecSinglePrecision;
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

/*
 * Copyright (C) 2021-2022 The ESPResSo project
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

#include "src/generated_kernels/Dynamic_UBB_double_precision.h"
#include "src/generated_kernels/Dynamic_UBB_single_precision.h"
#include "src/generated_kernels/InitialPDFsSetterDoublePrecision.h"
#include "src/generated_kernels/InitialPDFsSetterSinglePrecision.h"
#include "src/generated_kernels/StreamSweepDoublePrecision.h"
#include "src/generated_kernels/StreamSweepSinglePrecision.h"
#include "src/generated_kernels/macroscopic_values_accessors_double_precision.h"
#include "src/generated_kernels/macroscopic_values_accessors_single_precision.h"

#ifdef __AVX2__
#include "src/generated_kernels/CollideSweepDoublePrecisionLeesEdwardsAVX.h"
#include "src/generated_kernels/CollideSweepDoublePrecisionThermalizedAVX.h"
#include "src/generated_kernels/CollideSweepSinglePrecisionLeesEdwardsAVX.h"
#include "src/generated_kernels/CollideSweepSinglePrecisionThermalizedAVX.h"
#else
#include "src/generated_kernels/CollideSweepDoublePrecisionLeesEdwards.h"
#include "src/generated_kernels/CollideSweepDoublePrecisionThermalized.h"
#include "src/generated_kernels/CollideSweepSinglePrecisionLeesEdwards.h"
#include "src/generated_kernels/CollideSweepSinglePrecisionThermalized.h"
#endif

namespace walberla {

namespace detail {
template <typename FT = double> struct KernelTrait {
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
template <> struct KernelTrait<float> {
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

template <typename FT> struct BoundaryHandlingTrait {
  using Dynamic_UBB = lbm::Dynamic_UBB_double_precision;
};
template <> struct BoundaryHandlingTrait<float> {
  using Dynamic_UBB = lbm::Dynamic_UBB_single_precision;
};
} // namespace detail
} // namespace walberla

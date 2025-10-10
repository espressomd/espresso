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

#include "generated_kernels/DynamicUBBDoublePrecision.h"
#include "generated_kernels/DynamicUBBSinglePrecision.h"
#include "generated_kernels/FieldAccessorsDoublePrecision.h"
#include "generated_kernels/FieldAccessorsSinglePrecision.h"
#include "generated_kernels/InitialPDFsSetterDoublePrecision.h"
#include "generated_kernels/InitialPDFsSetterSinglePrecision.h"
#include "generated_kernels/PackInfoPdfDoublePrecision.h"
#include "generated_kernels/PackInfoPdfSinglePrecision.h"
#include "generated_kernels/PackInfoVecDoublePrecision.h"
#include "generated_kernels/PackInfoVecSinglePrecision.h"
#include "generated_kernels/UpdateVelFromPDFDoublePrecision.h"
#include "generated_kernels/UpdateVelFromPDFSinglePrecision.h"

#ifdef __AVX2__
#include "generated_kernels/StreamCollideSweepLeesEdwardsDoublePrecisionAVX.h"
#include "generated_kernels/StreamCollideSweepLeesEdwardsSinglePrecisionAVX.h"
#include "generated_kernels/StreamCollideSweepThermalizedDoublePrecisionAVX.h"
#include "generated_kernels/StreamCollideSweepThermalizedSinglePrecisionAVX.h"
#else
#include "generated_kernels/StreamCollideSweepLeesEdwardsDoublePrecision.h"
#include "generated_kernels/StreamCollideSweepLeesEdwardsSinglePrecision.h"
#include "generated_kernels/StreamCollideSweepThermalizedDoublePrecision.h"
#include "generated_kernels/StreamCollideSweepThermalizedSinglePrecision.h"
#endif

namespace walberla {
namespace detail {

using lbmpy::Arch;

template <typename FT = double, Arch AT = Arch::CPU> struct KernelTrait {
#ifdef __AVX2__
  using StreamCollisionModelThermalized =
      pystencils::StreamCollideSweepThermalizedDoublePrecisionAVX;
  using StreamCollisionModelLeesEdwards =
      pystencils::StreamCollideSweepLeesEdwardsDoublePrecisionAVX;
#else
  using StreamCollisionModelThermalized =
      pystencils::StreamCollideSweepThermalizedDoublePrecision;
  using StreamCollisionModelLeesEdwards =
      pystencils::StreamCollideSweepLeesEdwardsDoublePrecision;
#endif
  using InitialPDFsSetter = pystencils::InitialPDFsSetterDoublePrecision;
  using UpdateVelFromPDF = pystencils::UpdateVelFromPDFDoublePrecision;
  using PackInfoPdf = pystencils::PackInfoPdfDoublePrecision;
  using PackInfoVec = pystencils::PackInfoVecDoublePrecision;
};

template <> struct KernelTrait<float, Arch::CPU> {
#ifdef __AVX2__
  using StreamCollisionModelThermalized =
      pystencils::StreamCollideSweepThermalizedSinglePrecisionAVX;
  using StreamCollisionModelLeesEdwards =
      pystencils::StreamCollideSweepLeesEdwardsSinglePrecisionAVX;
#else
  using StreamCollisionModelThermalized =
      pystencils::StreamCollideSweepThermalizedSinglePrecision;
  using StreamCollisionModelLeesEdwards =
      pystencils::StreamCollideSweepLeesEdwardsSinglePrecision;
#endif
  using InitialPDFsSetter = pystencils::InitialPDFsSetterSinglePrecision;
  using UpdateVelFromPDF = pystencils::UpdateVelFromPDFSinglePrecision;
  using PackInfoPdf = pystencils::PackInfoPdfSinglePrecision;
  using PackInfoVec = pystencils::PackInfoVecSinglePrecision;
};

template <typename FT = double, Arch AT = Arch::CPU>
struct BoundaryHandlingTrait {
  using DynamicUBB = lbm::DynamicUBBDoublePrecision;
};

template <> struct BoundaryHandlingTrait<float, Arch::CPU> {
  using DynamicUBB = lbm::DynamicUBBSinglePrecision;
};

} // namespace detail
} // namespace walberla

#pragma once

#include "generated_kernels/InitialPDFsSetterDoublePrecision.h"
#include "generated_kernels/InitialPDFsSetterSinglePrecision.h"
#include "generated_kernels/StreamSweepDoublePrecision.h"
#include "generated_kernels/StreamSweepSinglePrecision.h"
#include "generated_kernels/macroscopic_values_accessors_double_precision.h"
#include "generated_kernels/macroscopic_values_accessors_single_precision.h"
#include "generated_kernels/CollideSweepDoublePrecision.h"
#include "generated_kernels/CollideSweepDoublePrecisionLeesEdwards.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalized.h"
#include "generated_kernels/CollideSweepSinglePrecision.h"
#include "generated_kernels/CollideSweepSinglePrecisionLeesEdwards.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalized.h"


#ifdef __AVX2__
#include "generated_kernels/CollideSweepDoublePrecisionAVX.h"
#include "generated_kernels/CollideSweepDoublePrecisionLeesEdwardsAVX.h"
#include "generated_kernels/CollideSweepDoublePrecisionThermalizedAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionLeesEdwardsAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionThermalizedAVX.h"
#include "generated_kernels/CollideSweepSinglePrecisionAVX.h"
#endif


namespace walberla {

namespace detail {
template <typename FT = double> struct KernelTrait {
#ifdef __AVX2__
  using ThermalizedCollisionModel =
      pystencils::CollideSweepDoublePrecisionThermalizedAVX;
  using UnthermalizedCollisionModel =
      pystencils::CollideSweepDoublePrecisionAVX;
  using LeesEdwardsCollisionModel =
      pystencils::CollideSweepDoublePrecisionLeesEdwardsAVX;
#else
  using ThermalizedCollisionModel =
      pystencils::CollideSweepDoublePrecisionThermalized;
  using UnthermalizedCollisionModel = pystencils::CollideSweepDoublePrecision;
  using LeesEdwardsCollisionModel =
      pystencils::CollideSweepDoublePrecisionLeesEdwards;
#endif
  using StreamSweep = pystencils::StreamSweepDoublePrecision;
  using InitialPDFsSetter = pystencils::InitialPDFsSetterDoublePrecision;
};
template <> struct KernelTrait<float> {
#ifdef __AVX2__
  using ThermalizedCollisionModel =
      pystencils::CollideSweepSinglePrecisionThermalizedAVX;
  using UnthermalizedCollisionModel =
      pystencils::CollideSweepSinglePrecisionAVX;
  using LeesEdwardsCollisionModel =
      pystencils::CollideSweepSinglePrecisionLeesEdwardsAVX;
#else
  using ThermalizedCollisionModel =
      pystencils::CollideSweepSinglePrecisionThermalized;
  using UnthermalizedCollisionModel = pystencils::CollideSweepSinglePrecision;
  using LeesEdwardsCollisionModel =
      pystencils::CollideSweepSinglePrecisionLeesEdwards;
#endif
  using StreamSweep = pystencils::StreamSweepSinglePrecision;
  using InitialPDFsSetter = pystencils::InitialPDFsSetterSinglePrecision;
};
} // namespace detail
}

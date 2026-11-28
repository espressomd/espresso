/*
 * Copyright (C) 2019-2026 The ESPResSo project
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

#include "lb_kernels.cuh"

#include <walberla_bridge/Architecture.hpp>

#include <blockforest/communication/UniformBufferedScheme.h>
#include <gpu/GPUField.h>
#include <gpu/communication/MemcpyPackInfo.h>
#include <gpu/communication/UniformGPUScheme.h>

namespace walberla {
template <typename FT, class PdfStencil>
struct FieldTrait<FT, PdfStencil, lbmpy::Arch::GPU> {
private:
  static auto constexpr AT = lbmpy::Arch::GPU;
  template <class Field>
  using MemcpyPackInfo = gpu::communication::MemcpyPackInfo<Field>;

public:
  template <typename Stencil>
  class UniformGPUScheme
      : public gpu::communication::UniformGPUScheme<Stencil> {
  public:
    explicit UniformGPUScheme(auto const &bf)
        : gpu::communication::UniformGPUScheme<Stencil>(
              bf, /* sendDirectlyFromGPU */ false,
              /* useLocalCommunication */ false) {}
  };
  using PdfField = gpu::GPUField<FT>;
  using VectorField = gpu::GPUField<FT>;
  template <class Field> using PackInfo = MemcpyPackInfo<Field>;
  using PackInfoStreamingPdf = detail::KernelTrait<FT, AT>::PackInfoPdf;
  using PackInfoStreamingVec = detail::KernelTrait<FT, AT>::PackInfoVec;
  template <class Stencil> using RegularCommScheme = UniformGPUScheme<Stencil>;
  template <class Stencil>
  using BoundaryCommScheme =
      blockforest::communication::UniformBufferedScheme<Stencil>;
};
} // namespace walberla

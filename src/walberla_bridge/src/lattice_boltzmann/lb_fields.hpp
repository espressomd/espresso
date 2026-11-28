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

#include <blockforest/communication/UniformBufferedScheme.h>
#include <field/GhostLayerField.h>
#include <field/communication/PackInfo.h>

#include <walberla_bridge/Architecture.hpp>

#include "lb_kernels.hpp"

namespace walberla {
template <typename FT, class PdfStencil, lbmpy::Arch AT = lbmpy::Arch::CPU>
struct FieldTrait {
  using PdfField = field::GhostLayerField<FT, PdfStencil::Size>;
  using VectorField = field::GhostLayerField<FT, uint_t{3u}>;
  template <class Field> using PackInfo = field::communication::PackInfo<Field>;
  using PackInfoStreamingPdf = detail::KernelTrait<FT, AT>::PackInfoPdf;
  using PackInfoStreamingVec = detail::KernelTrait<FT, AT>::PackInfoVec;
  template <class Stencil>
  using RegularCommScheme =
      blockforest::communication::UniformBufferedScheme<Stencil>;
  template <class Stencil>
  using BoundaryCommScheme =
      blockforest::communication::UniformBufferedScheme<Stencil>;
};
} // namespace walberla

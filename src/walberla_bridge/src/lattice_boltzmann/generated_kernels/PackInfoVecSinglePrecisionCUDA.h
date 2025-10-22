//======================================================================================================================
//
//  This file is part of waLBerla. waLBerla is free software: you can
//  redistribute it and/or modify it under the terms of the GNU General Public
//  License as published by the Free Software Foundation, either version 3 of
//  the License, or (at your option) any later version.
//
//  waLBerla is distributed in the hope that it will be useful, but WITHOUT
//  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
//  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
//  for more details.
//
//  You should have received a copy of the GNU General Public License along
//  with waLBerla (see COPYING.txt). If not, see <http://www.gnu.org/licenses/>.
//
//! \\file PackInfoVecSinglePrecisionCUDA.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7, sympy v1.12.1,
// lbmpy_walberla/pystencils_walberla from waLBerla commit
// 0aab9c0af2335b1f6fec75deae06e514ccb233ab

#pragma once

#include "core/DataTypes.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"
#include "gpu/communication/GeneratedGPUPackInfo.h"

namespace walberla {
namespace pystencils {

class PackInfoVecSinglePrecisionCUDA
    : public ::walberla::gpu::GeneratedGPUPackInfo {
public:
  PackInfoVecSinglePrecisionCUDA(BlockDataID fieldID_) : fieldID(fieldID_) {}
  ~PackInfoVecSinglePrecisionCUDA() override = default;

  void pack(stencil::Direction dir, unsigned char *buffer, IBlock *block,
            gpuStream_t stream) override;
  void communicateLocal(stencil::Direction /*dir*/, const IBlock * /* sender */,
                        IBlock * /* receiver */,
                        gpuStream_t /* stream */) override {
    WALBERLA_ABORT(
        "Local Communication not implemented yet for standard PackInfos. To "
        "run your application, turn off local communication in the "
        "communication class, e.g. with useLocalCommunication=false")
  }
  void unpack(stencil::Direction dir, unsigned char *buffer, IBlock *block,
              gpuStream_t stream) override;
  uint_t size(stencil::Direction dir, IBlock *block) override;

private:
  BlockDataID fieldID;
};

} // namespace pystencils
} // namespace walberla

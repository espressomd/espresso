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
//! \\file PackInfoPdfDoublePrecisionCUDA.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34,
// sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit
// 007e77e077ad9d22b5eed6f3d3118240993e553c

#pragma once

#include "core/DataTypes.h"

#include "domain_decomposition/IBlock.h"

#include "stencil/Directions.h"

#include "gpu/GPUField.h"
#include "gpu/GPUWrapper.h"
#include "gpu/communication/GeneratedGPUPackInfo.h"

namespace walberla {
namespace pystencils {

class PackInfoPdfDoublePrecisionCUDA
    : public ::walberla::gpu::GeneratedGPUPackInfo {
public:
  PackInfoPdfDoublePrecisionCUDA(BlockDataID pdfsID_) : pdfsID(pdfsID_) {}
  ~PackInfoPdfDoublePrecisionCUDA() override = default;

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
  BlockDataID pdfsID;
};

} // namespace pystencils
} // namespace walberla

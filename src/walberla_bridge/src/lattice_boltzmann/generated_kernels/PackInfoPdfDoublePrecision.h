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
//! \\file PackInfoPdfDoublePrecision.h
//! \\author pystencils
//======================================================================================================================

// kernel generated with pystencils v1.4+1.ge851f4e, lbmpy v1.4+1.ge9efe34,
// sympy v1.12.1, lbmpy_walberla/pystencils_walberla from waLBerla commit
// 007e77e077ad9d22b5eed6f3d3118240993e553c

#pragma once
#include "communication/UniformPackInfo.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "domain_decomposition/IBlock.h"
#include "field/GhostLayerField.h"
#include "stencil/Directions.h"

#include <memory>

#define FUNC_PREFIX

#ifdef __GNUC__
#define RESTRICT __restrict__
#else
#define RESTRICT
#endif

namespace walberla {
namespace pystencils {

class PackInfoPdfDoublePrecision
    : public ::walberla::communication::UniformPackInfo {
public:
  PackInfoPdfDoublePrecision(BlockDataID pdfsID_) : pdfsID(pdfsID_) {}
  ~PackInfoPdfDoublePrecision() override = default;

  bool constantDataExchange() const override { return true; }
  bool threadsafeReceiving() const override { return true; }

  void unpackData(IBlock *receiver, stencil::Direction dir,
                  mpi::RecvBuffer &buffer) override {
    const auto dataSize = size(dir, receiver);
    auto bufferSize = dataSize + sizeof(double);
    auto bufferPtr = reinterpret_cast<void *>(buffer.skip(bufferSize));
    std::align(alignof(double), dataSize, bufferPtr, bufferSize);
    unpack(dir, reinterpret_cast<unsigned char *>(bufferPtr), receiver);
  }

  void communicateLocal(const IBlock *sender, IBlock *receiver,
                        stencil::Direction dir) override {
    mpi::SendBuffer sBuffer;
    packData(sender, dir, sBuffer);
    mpi::RecvBuffer rBuffer(sBuffer);
    unpackData(receiver, stencil::inverseDir[dir], rBuffer);
  }

  void packDataImpl(const IBlock *sender, stencil::Direction dir,
                    mpi::SendBuffer &outBuffer) const override {
    const auto dataSize = size(dir, sender);
    auto bufferSize = dataSize + sizeof(double);
    auto bufferPtr = reinterpret_cast<void *>(outBuffer.forward(bufferSize));
    std::align(alignof(double), dataSize, bufferPtr, bufferSize);
    pack(dir, reinterpret_cast<unsigned char *>(bufferPtr),
         const_cast<IBlock *>(sender));
  }

  void pack(stencil::Direction dir, unsigned char *buffer, IBlock *block) const;
  void unpack(stencil::Direction dir, unsigned char *buffer,
              IBlock *block) const;
  uint_t size(stencil::Direction dir, const IBlock *block) const;

private:
  BlockDataID pdfsID;
};

} // namespace pystencils
} // namespace walberla

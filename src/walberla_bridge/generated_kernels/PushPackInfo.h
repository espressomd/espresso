// kernel generated with pystencils v0.3.4+4.g4fecf0c, lbmpy v0.3.4+6.g2faceda,
// lbmpy_walberla/pystencils_walberla from commit
// b17ca5caf00db7d19f86c5f85c6f67fec6c16aff

#pragma once
#include "communication/UniformPackInfo.h"
#include "core/DataTypes.h"
#include "core/cell/CellInterval.h"
#include "domain_decomposition/IBlock.h"
#include "field/GhostLayerField.h"
#include "stencil/Directions.h"

#define FUNC_PREFIX

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

namespace walberla {
namespace lbm {

class PushPackInfo : public ::walberla::communication::UniformPackInfo {
public:
  PushPackInfo(BlockDataID pdfsID_) : pdfsID(pdfsID_){};
  virtual ~PushPackInfo() {}

  bool constantDataExchange() const { return true; }
  bool threadsafeReceiving() const { return true; }

  void unpackData(IBlock *receiver, stencil::Direction dir,
                  mpi::RecvBuffer &buffer) {
    const auto dataSize = size(dir, receiver);
    unpack(dir, buffer.skip(dataSize), receiver);
  }

  void communicateLocal(const IBlock *sender, IBlock *receiver,
                        stencil::Direction dir) {
    // TODO: optimize by generating kernel for this case
    mpi::SendBuffer sBuffer;
    packData(sender, dir, sBuffer);
    mpi::RecvBuffer rBuffer(sBuffer);
    unpackData(receiver, stencil::inverseDir[dir], rBuffer);
  }

private:
  void packDataImpl(const IBlock *sender, stencil::Direction dir,
                    mpi::SendBuffer &outBuffer) const {
    const auto dataSize = size(dir, sender);
    pack(dir, outBuffer.forward(dataSize), const_cast<IBlock *>(sender));
  }

  void pack(stencil::Direction dir, unsigned char *buffer, IBlock *block) const;
  void unpack(stencil::Direction dir, unsigned char *buffer,
              IBlock *block) const;
  uint_t size(stencil::Direction dir, const IBlock *block) const;

  BlockDataID pdfsID;
};

} // namespace lbm
} // namespace walberla
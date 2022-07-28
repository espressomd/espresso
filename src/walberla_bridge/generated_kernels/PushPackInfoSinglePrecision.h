// kernel generated with pystencils v1.0+12.g54b91e2, lbmpy v1.0+8.gac750b5,
// lbmpy_walberla/pystencils_walberla from commit
// e1fe2ad1dcbe8f31ea79d95e8a5a5cc0ee3691f3

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

class PushPackInfoSinglePrecision
    : public ::walberla::communication::UniformPackInfo {
public:
  PushPackInfoSinglePrecision(BlockDataID pdfsID_) : pdfsID(pdfsID_){};
  virtual ~PushPackInfoSinglePrecision() {}

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

  void packDataImpl(const IBlock *sender, stencil::Direction dir,
                    mpi::SendBuffer &outBuffer) const {
    const auto dataSize = size(dir, sender);
    pack(dir, outBuffer.forward(dataSize), const_cast<IBlock *>(sender));
  }

  void pack(stencil::Direction dir, unsigned char *buffer, IBlock *block) const;
  void unpack(stencil::Direction dir, unsigned char *buffer,
              IBlock *block) const;
  uint_t size(stencil::Direction dir, const IBlock *block) const;

private:
  BlockDataID pdfsID;
};

} // namespace lbm
} // namespace walberla

// kernel generated with pystencils v1.2, lbmpy v1.2,
// lbmpy_walberla/pystencils_walberla from waLBerla commit ref:
// a839fac6ef7d0c58e7710e4d50490e9dd7146b4a

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
namespace pystencils {

class DensityPackInfo_single_precision
    : public ::walberla::communication::UniformPackInfo {
public:
  DensityPackInfo_single_precision(BlockDataID jID_) : jID(jID_){};
  virtual ~DensityPackInfo_single_precision() {}

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
  BlockDataID jID;
};

} // namespace pystencils
} // namespace walberla
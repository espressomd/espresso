// kernel generated with pystencils v0.4.4, lbmpy v0.4.4, lbmpy_walberla/pystencils_walberla from commit 08f04ef64f95609b47838db85862033a1600afa1

#pragma once
#include "stencil/Directions.h"
#include "core/cell/CellInterval.h"
#include "core/DataTypes.h"
#include "field/GhostLayerField.h"
#include "domain_decomposition/IBlock.h"
#include "communication/UniformPackInfo.h"

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


class DensityPackInfo_double_precision : public ::walberla::communication::UniformPackInfo
{
public:
    DensityPackInfo_double_precision( BlockDataID jID_ )
        : jID(jID_)
    {};
    virtual ~DensityPackInfo_double_precision() {}

   bool constantDataExchange() const { return true; }
   bool threadsafeReceiving()  const { return true; }

   void unpackData(IBlock * receiver, stencil::Direction dir, mpi::RecvBuffer & buffer) {
        const auto dataSize = size(dir, receiver);
        unpack(dir, buffer.skip(dataSize), receiver);
   }

   void communicateLocal(const IBlock * sender, IBlock * receiver, stencil::Direction dir) {
       //TODO: optimize by generating kernel for this case
       mpi::SendBuffer sBuffer;
       packData( sender, dir, sBuffer );
       mpi::RecvBuffer rBuffer( sBuffer );
       unpackData( receiver, stencil::inverseDir[dir], rBuffer );
   }

private:
   void packDataImpl(const IBlock * sender, stencil::Direction dir, mpi::SendBuffer & outBuffer) const {
        const auto dataSize = size(dir, sender);
        pack(dir, outBuffer.forward(dataSize), const_cast<IBlock*>(sender));
   }

   void pack  (stencil::Direction dir, unsigned char * buffer, IBlock * block) const;
   void unpack(stencil::Direction dir, unsigned char * buffer, IBlock * block) const;
   uint_t size  (stencil::Direction dir, const IBlock * block) const;

    BlockDataID jID;
};


} // namespace pystencils
} // namespace walberla
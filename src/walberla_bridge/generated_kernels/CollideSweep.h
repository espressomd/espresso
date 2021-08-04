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
//! \\file CollideSweep.h
//! \\author pystencils
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"

#include "blockforest/StructuredBlockForest.h"
#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "domain_decomposition/StructuredBlockStorage.h"
#include <memory>
#include <set>



#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic push
#   pragma GCC diagnostic ignored "-Wunused-parameter"
#   pragma GCC diagnostic ignored "-Wreorder"
#endif

namespace walberla {
namespace pystencils {


class CollideSweep
{
public:
    CollideSweep( IBlock * const block, BlockDataID forceID_, BlockDataID pdfsID_, double omega_bulk, double omega_even, double omega_odd, double omega_shear)
        : block_(block), forceID(forceID_), pdfsID(pdfsID_), omega_bulk_(omega_bulk), omega_even_(omega_even), omega_odd_(omega_odd), omega_shear_(omega_shear)
    {};

    

    void run(IBlock * block);
    
    void runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block);

    
    void operator() ()
    {
        run(block_);
    }
    

    static std::function<void (IBlock *)> getSweep(const shared_ptr<CollideSweep> & kernel)
    {
        return [kernel] 
               (IBlock * b) 
               { kernel->run(b); };
    }

    static std::function<void (IBlock*)> getSweepOnCellInterval(const shared_ptr<CollideSweep> & kernel, const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1)
    {
        return [kernel, blocks, globalCellInterval, ghostLayers]
               (IBlock * b) 
               { kernel->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b); };
    }

    std::function<void (IBlock *)> getSweep()
    {
        return [this] 
               (IBlock * b) 
               { this->run(b); };
    }

    std::function<void (IBlock *)> getSweepOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1)
    {
        return [this, blocks, globalCellInterval, ghostLayers]
               (IBlock * b) 
               { this->runOnCellInterval(blocks, globalCellInterval, ghostLayers, b); };
    }


    IBlock * block_;
    BlockDataID forceID;
    BlockDataID pdfsID;
    double omega_bulk_;
    double omega_even_;
    double omega_odd_;
    double omega_shear_;

};

struct CollideSweepFactory {
   CollideSweepFactory( std::shared_ptr<blockforest::StructuredBlockForest>, BlockDataID forceID, BlockDataID pdfsID, double, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t, uint32_t )
        : m_forceID(forceID), m_pdfsID(pdfsID), m_omega_bulk(omega_bulk), m_omega_even(omega_even), m_omega_odd(omega_odd), m_omega_shear(omega_shear) {}

   CollideSweep * operator() ( IBlock * const block ) {
      return new CollideSweep( block, m_forceID, m_pdfsID,
        m_omega_bulk, m_omega_even, m_omega_odd, m_omega_shear);
   }

   private:
    BlockDataID m_forceID;
    BlockDataID m_pdfsID;
    double m_omega_bulk;
    double m_omega_even;
    double m_omega_odd;
    double m_omega_shear;
};

} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif
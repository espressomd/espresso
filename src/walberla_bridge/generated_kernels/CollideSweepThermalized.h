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
//! \\file CollideSweepThermalized.h
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


class CollideSweepThermalized
{
public:
    CollideSweepThermalized( IBlock * const block, BlockDataID forceID_, BlockDataID pdfsID_, uint32_t block_offset_0, uint32_t block_offset_1, uint32_t block_offset_2, double kT, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t seed, uint32_t time_step )
        : block_(block), forceID(forceID_), pdfsID(pdfsID_), block_offset_0_(block_offset_0), block_offset_1_(block_offset_1), block_offset_2_(block_offset_2), kT_(kT), omega_bulk_(omega_bulk), omega_even_(omega_even), omega_odd_(omega_odd), omega_shear_(omega_shear), seed_(seed), time_step_(time_step)
    {};

    

    void run(IBlock * block);
    
    void runOnCellInterval(const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers, IBlock * block);

    
    void operator() ()
    {
        run(block_);
    }
    

    static std::function<void (IBlock *)> getSweep(const shared_ptr<CollideSweepThermalized> & kernel)
    {
        return [kernel] 
               (IBlock * b) 
               { kernel->run(b); };
    }

    static std::function<void (IBlock*)> getSweepOnCellInterval(const shared_ptr<CollideSweepThermalized> & kernel, const shared_ptr<StructuredBlockStorage> & blocks, const CellInterval & globalCellInterval, cell_idx_t ghostLayers=1)
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
    uint32_t block_offset_0_;
    uint32_t block_offset_1_;
    uint32_t block_offset_2_;
    double kT_;
    double omega_bulk_;
    double omega_even_;
    double omega_odd_;
    double omega_shear_;
    uint32_t seed_;
    uint32_t time_step_;

};

struct CollideSweepThermalizedFactory {
   CollideSweepThermalizedFactory( std::shared_ptr<blockforest::StructuredBlockForest> blockforest, BlockDataID forceID, BlockDataID pdfsID, double kT, double omega_bulk, double omega_even, double omega_odd, double omega_shear, uint32_t seed, uint32_t time_step )
        : m_blockforest(blockforest), m_forceID(forceID), m_pdfsID(pdfsID), m_kT(kT), m_omega_bulk(omega_bulk), m_omega_even(omega_even), m_omega_odd(omega_odd), m_omega_shear(omega_shear), m_seed(seed), m_time_step(time_step) {}

   CollideSweepThermalized * operator() ( IBlock * const block ) {
      return new CollideSweepThermalized( block, m_forceID, m_pdfsID,
        uint32_c(m_blockforest->getBlockCellBB(*block).xMin()),
        uint32_c(m_blockforest->getBlockCellBB(*block).yMin()),
        uint32_c(m_blockforest->getBlockCellBB(*block).zMin()),
        m_kT, m_omega_bulk, m_omega_even, m_omega_odd, m_omega_shear, m_seed, m_time_step );
   }

   private:
    std::shared_ptr<blockforest::StructuredBlockForest> m_blockforest;
    BlockDataID m_forceID;
    BlockDataID m_pdfsID;
    double m_kT;
    double m_omega_bulk;
    double m_omega_even;
    double m_omega_odd;
    double m_omega_shear;
    uint32_t m_seed;
    uint32_t m_time_step;
};

} // namespace pystencils
} // namespace walberla


#if ( defined WALBERLA_CXX_COMPILER_IS_GNU ) || ( defined WALBERLA_CXX_COMPILER_IS_CLANG )
#   pragma GCC diagnostic pop
#endif
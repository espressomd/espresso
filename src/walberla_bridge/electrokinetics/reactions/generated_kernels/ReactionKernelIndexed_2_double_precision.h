// kernel generated with pystencils v1.0, lbmpy v1.0, lbmpy_walberla/pystencils_walberla from commit 01a28162ae1aacf7b96152c9f886ce54cc7f53ff

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
//! \\file ReactionKernelIndexed_2_double_precision.h
//! \\author pystencils
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"

#include "field/GhostLayerField.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "blockforest/StructuredBlockForest.h"
#include "field/FlagField.h"
#include "core/debug/Debug.h"

#include <set>
#include <vector>



#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

namespace walberla {
namespace pystencils {


class ReactionKernelIndexed_2_double_precision
{
public:
    struct IndexInfo { 
        int32_t x;
        int32_t y;
        int32_t z;
        IndexInfo(int32_t x_, int32_t y_, int32_t z_) : x(x_), y(y_), z(z_) {}
        bool operator==(const IndexInfo & o) const {
            return x == o.x && y == o.y && z == o.z;
        }
    };



    class IndexVectors
    {
    public:
        using CpuIndexVector = std::vector<IndexInfo>;

        enum Type {
            ALL = 0,
            INNER = 1,
            OUTER = 2,
            NUM_TYPES = 3
        };

        IndexVectors() = default;
        bool operator==(IndexVectors & other) { return other.cpuVectors_ == cpuVectors_; }

        CpuIndexVector & indexVector(Type t) { return cpuVectors_[t]; }
        IndexInfo * pointerCpu(Type t)  { return &(cpuVectors_[t][0]); }

        void syncGPU()
        {
            
        }

    private:
        std::vector<CpuIndexVector> cpuVectors_{NUM_TYPES};

        
    };

    ReactionKernelIndexed_2_double_precision( const shared_ptr<StructuredBlockForest> & blocks,
                   BlockDataID rho_0ID_, BlockDataID rho_1ID_, double order_0, double order_1, double rate_coefficient, double stoech_0, double stoech_1)
        : rho_0ID(rho_0ID_), rho_1ID(rho_1ID_), order_0_(order_0), order_1_(order_1), rate_coefficient_(rate_coefficient), stoech_0_(stoech_0), stoech_1_(stoech_1)
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_ReactionKernelIndexed_2_double_precision");
    };
    
    ReactionKernelIndexed_2_double_precision(BlockDataID indexVectorID_, BlockDataID rho_0ID_, BlockDataID rho_1ID_, double order_0, double order_1, double rate_coefficient, double stoech_0, double stoech_1)
        :  indexVectorID(indexVectorID_), rho_0ID(rho_0ID_), rho_1ID(rho_1ID_), order_0_(order_0), order_1_(order_1), rate_coefficient_(rate_coefficient), stoech_0_(stoech_0), stoech_1_(stoech_1)
    {};

    void run (IBlock * block);

    void operator() (IBlock * block)
    {
        run(block);
    }

    void inner (IBlock * block);

    void outer (IBlock * block);

    std::function<void (IBlock *)> getSweep()
    {
        return [this]
               (IBlock * b)
               { this->run(b); };
    }

    std::function<void (IBlock *)> getInnerSweep()
    {
        return [this]
               (IBlock * b)
               { this->inner(b); };
    }

    std::function<void (IBlock *)> getOuterSweep()
    {
        return [this]
               (IBlock * b)
               { this->outer(b); };
    }

    template<typename FlagField_T>
    void fillFromFlagField( const shared_ptr<StructuredBlockForest> & blocks, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID)
    {
        for( auto blockIt = blocks->begin(); blockIt != blocks->end(); ++blockIt )
            fillFromFlagField<FlagField_T>(&*blockIt, flagFieldID, boundaryFlagUID, domainFlagUID );
    }


    template<typename FlagField_T>
    void fillFromFlagField(IBlock * block, ConstBlockDataID flagFieldID,
                            FlagUID boundaryFlagUID, FlagUID domainFlagUID )
    {
        auto * indexVectors = block->getData< IndexVectors > ( indexVectorID );
        auto & indexVectorAll = indexVectors->indexVector(IndexVectors::ALL);
        auto & indexVectorInner = indexVectors->indexVector(IndexVectors::INNER);
        auto & indexVectorOuter = indexVectors->indexVector(IndexVectors::OUTER);

        auto * flagField = block->getData< FlagField_T > ( flagFieldID );
        

        if( !(flagField->flagExists(boundaryFlagUID) && flagField->flagExists(domainFlagUID) ))
            return;

        auto boundaryFlag = flagField->getFlag(boundaryFlagUID);
        auto domainFlag = flagField->getFlag(domainFlagUID);

        auto inner = flagField->xyzSize();
        inner.expand( cell_idx_t(-1) );

        indexVectorAll.clear();
        indexVectorInner.clear();
        indexVectorOuter.clear();

        
        auto flagWithGLayers = flagField->xyzSizeWithGhostLayer();
        real_t dot = 0.0; real_t maxn = 0.0;
        cell_idx_t calculated_idx = 0;
        cell_idx_t dx = 0; cell_idx_t dy = 0;   cell_idx_t dz = 0; 
        cell_idx_t sum_x = 0; cell_idx_t sum_y = 0;  cell_idx_t sum_z = 0; 
        for( auto it = flagField->beginWithGhostLayerXYZ(); it != flagField->end(); ++it )
        {
            
            if( ! isFlagSet(it, boundaryFlag) )
                continue;
            if ( flagWithGLayers.contains(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(0)) && isFlagSet( it.neighbor(0, 0, 0 , 0 ), domainFlag ) )
            {
                
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  0 );
                
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
                
            }
            

        }
        

        indexVectors->syncGPU();
    }

private:
    void run_impl(IBlock * block, IndexVectors::Type type);

    BlockDataID indexVectorID;
    
public:
    BlockDataID rho_0ID;
    BlockDataID rho_1ID;
    double order_0_;
    double order_1_;
    double rate_coefficient_;
    double stoech_0_;
    double stoech_1_;
};



} // namespace pystencils
} // namespace walberla
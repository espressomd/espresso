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
//! \\file Dirichlet.h
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


class Dirichlet
{
public:
    struct IndexInfo { 
        int32_t x;
        int32_t y;
        int32_t z;
        int32_t dir;
        double value;
        IndexInfo(int32_t x_, int32_t y_, int32_t z_, int32_t dir_) : x(x_), y(y_), z(z_), dir(dir_), value() {}
        bool operator==(const IndexInfo & o) const {
            return x == o.x && y == o.y && z == o.z && dir == o.dir && floatIsEqual(value, o.value);
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

    Dirichlet( const shared_ptr<StructuredBlockForest> & blocks,
                   BlockDataID fieldID_, std::function<real_t(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)>& dirichletCallback )
        :elementInitaliser(dirichletCallback), fieldID(fieldID_)
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_Dirichlet");
    };

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
            fillFromFlagField<FlagField_T>(blocks, &*blockIt, flagFieldID, boundaryFlagUID, domainFlagUID );
    }


    template<typename FlagField_T>
    void fillFromFlagField( const shared_ptr<StructuredBlockForest> &blocks, IBlock * block, ConstBlockDataID flagFieldID,
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
            sum_x = 0; sum_y = 0;  sum_z = 0; 
            
            if( ! isFlagSet(it, boundaryFlag) )
                continue;
            if ( flagWithGLayers.contains(it.x() + cell_idx_c(-1), it.y() + cell_idx_c(0), it.z() + cell_idx_c(0)) && isFlagSet( it.neighbor(-1, 0, 0 , 0 ), domainFlag ) )
            {
                sum_x += cell_idx_c(-1); sum_y += cell_idx_c(0);  sum_z += cell_idx_c(0); 
                
            }
            
            if ( flagWithGLayers.contains(it.x() + cell_idx_c(0), it.y() + cell_idx_c(-1), it.z() + cell_idx_c(0)) && isFlagSet( it.neighbor(0, -1, 0 , 0 ), domainFlag ) )
            {
                sum_x += cell_idx_c(0); sum_y += cell_idx_c(-1);  sum_z += cell_idx_c(0); 
                
            }
            
            if ( flagWithGLayers.contains(it.x() + cell_idx_c(0), it.y() + cell_idx_c(0), it.z() + cell_idx_c(-1)) && isFlagSet( it.neighbor(0, 0, -1 , 0 ), domainFlag ) )
            {
                sum_x += cell_idx_c(0); sum_y += cell_idx_c(0);  sum_z += cell_idx_c(-1); 
                
            }
            

        
            dot = 0.0; maxn = 0.0; calculated_idx = 0;
            if(sum_x != 0 or sum_y !=0  or sum_z !=0 )
            {
                dx = -1; dy = 0;  dz = 0; 
                dot = real_c( dx*sum_x + dy*sum_y  + dz*sum_z );
                if (dot > maxn)
                {
                    maxn = dot;
                    calculated_idx = 0;
                }
            
                dx = 0; dy = -1;  dz = 0; 
                dot = real_c( dx*sum_x + dy*sum_y  + dz*sum_z );
                if (dot > maxn)
                {
                    maxn = dot;
                    calculated_idx = 1;
                }
            
                dx = 0; dy = 0;  dz = -1; 
                dot = real_c( dx*sum_x + dy*sum_y  + dz*sum_z );
                if (dot > maxn)
                {
                    maxn = dot;
                    calculated_idx = 2;
                }
            
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  calculated_idx );
                real_t InitialisatonAdditionalData = elementInitaliser(Cell(it.x(), it.y(), it.z()), blocks, *block);
                element.value = InitialisatonAdditionalData;
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
    std::function<real_t(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)> elementInitaliser; 
public:
    BlockDataID fieldID;
};



} // namespace pystencils
} // namespace walberla
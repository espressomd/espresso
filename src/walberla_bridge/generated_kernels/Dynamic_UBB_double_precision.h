// kernel generated with pystencils v0.4.3+4.g30da657, lbmpy v0.4.3+2.g0e17e61, lbmpy_walberla/pystencils_walberla from commit 88f85eb7a979f81d68e76009811aeed53ec3014e

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
//! \\file Dynamic_UBB_double_precision.h
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
namespace lbm {


class Dynamic_UBB_double_precision
{
public:
    struct IndexInfo { 
        int32_t x;
        int32_t y;
        int32_t z;
        int32_t dir;
        double vel_0;
        double vel_1;
        double vel_2;
        IndexInfo(int32_t x_, int32_t y_, int32_t z_, int32_t dir_) : x(x_), y(y_), z(z_), dir(dir_), vel_0(), vel_1(), vel_2() {}
        bool operator==(const IndexInfo & o) const {
            return x == o.x && y == o.y && z == o.z && dir == o.dir && floatIsEqual(vel_0, o.vel_0) && floatIsEqual(vel_1, o.vel_1) && floatIsEqual(vel_2, o.vel_2);
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
        bool operator==(IndexVectors const &other) const { return other.cpuVectors_ == cpuVectors_; }

        CpuIndexVector & indexVector(Type t) { return cpuVectors_[t]; }
        IndexInfo * pointerCpu(Type t)  { return cpuVectors_[t].data(); }

        void syncGPU()
        {
            
        }

    private:
        std::vector<CpuIndexVector> cpuVectors_{NUM_TYPES};

        
    };

    Dynamic_UBB_double_precision( const shared_ptr<StructuredBlockForest> & blocks,
                   BlockDataID pdfsID_, std::function<Vector3<double>(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)>& velocityCallback )
        :elementInitaliser(velocityCallback), pdfsID(pdfsID_)
    {
        auto createIdxVector = []( IBlock * const , StructuredBlockStorage * const ) { return new IndexVectors(); };
        indexVectorID = blocks->addStructuredBlockData< IndexVectors >( createIdxVector, "IndexField_Dynamic_UBB_double_precision");
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

        for( auto it = flagField->begin(); it != flagField->end(); ++it )
        {
            if( ! isFlagSet(it, domainFlag) )
                continue;
            if ( isFlagSet( it.neighbor(0, 0, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  0 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()+0, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  1 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()+1, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, -1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  2 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()-1, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, 0, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  3 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()-1, it.y()+0, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, 0, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  4 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+1, it.y()+0, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 0, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  5 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()+0, it.z()+1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 0, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  6 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()+0, it.z()-1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, 1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  7 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()-1, it.y()+1, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, 1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  8 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+1, it.y()+1, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, -1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  9 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()-1, it.y()-1, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, -1, 0 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  10 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+1, it.y()-1, it.z()+0), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 1, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  11 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()+1, it.z()+1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, -1, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  12 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()-1, it.z()+1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, 0, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  13 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()-1, it.y()+0, it.z()+1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, 0, 1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  14 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+1, it.y()+0, it.z()+1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, 1, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  15 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()+1, it.z()-1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(0, -1, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  16 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+0, it.y()-1, it.z()-1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(-1, 0, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  17 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()-1, it.y()+0, it.z()-1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
                indexVectorAll.push_back( element );
                if( inner.contains( it.x(), it.y(), it.z() ) )
                    indexVectorInner.push_back( element );
                else
                    indexVectorOuter.push_back( element );
            }
            
            if ( isFlagSet( it.neighbor(1, 0, -1 , 0 ), boundaryFlag ) )
            {
                auto element = IndexInfo(it.x(), it.y(),  it.z(),  18 );
                Vector3<double> InitialisatonAdditionalData = elementInitaliser(Cell(it.x()+1, it.y()+0, it.z()-1), blocks, *block);
                element.vel_0 = InitialisatonAdditionalData[0];
                element.vel_1 = InitialisatonAdditionalData[1];
                element.vel_2 = InitialisatonAdditionalData[2];
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
    std::function<Vector3<double>(const Cell &, const shared_ptr<StructuredBlockForest>&, IBlock&)> elementInitaliser; 
public:
    BlockDataID pdfsID;
};



} // namespace lbm
} // namespace walberla
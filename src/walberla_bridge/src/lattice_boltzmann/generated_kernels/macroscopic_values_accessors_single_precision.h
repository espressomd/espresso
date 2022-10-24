// kernel generated with pystencils v1.0+21.g8bd3cef, lbmpy v1.0+8.gac750b5, lbmpy_walberla/pystencils_walberla from commit e1fe2ad1dcbe8f31ea79d95e8a5a5cc0ee3691f3

// file extracted from the walberla lattice model
// (python/lbmpy_walberla/templates/LatticeModel.tmpl.h)

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
//! \\author Martin Bauer <martin.bauer@fau.de>
//
//======================================================================================================================

#pragma once
#include "core/DataTypes.h"
#include "core/logging/Logging.h"
#include "core/math/Matrix3.h"

#include "field/GhostLayerField.h"
#include "field/SwapableCompare.h"
#include "domain_decomposition/BlockDataID.h"
#include "domain_decomposition/IBlock.h"
#include "stencil/D3Q19.h"

#include <type_traits>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif

namespace walberla {
namespace lbm {
namespace accessor {


//======================================================================================================================
//
//  Implementation of macroscopic value backend
//
//======================================================================================================================

namespace EquilibriumDistribution
{
   float get( const stencil::Direction direction,
               const Vector3< float > & u = Vector3< float >( float(0.0) ),
               float rho = float(1.0) )
   {
        
        
    using namespace stencil;
    switch( direction ) {
        case C: return rho*-0.33333333333333331f*(u[0]*u[0]) + rho*-0.33333333333333331f*(u[1]*u[1]) + rho*-0.33333333333333331f*(u[2]*u[2]) + rho*0.33333333333333331f;
        case N: return rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[1] + rho*0.16666666666666666f*(u[1]*u[1]);
        case S: return rho*-0.16666666666666666f*u[1] + rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[1]*u[1]);
        case W: return rho*-0.16666666666666666f*u[0] + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[0]*u[0]);
        case E: return rho*-0.16666666666666666f*(u[1]*u[1]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[0] + rho*0.16666666666666666f*(u[0]*u[0]);
        case T: return rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[2] + rho*0.16666666666666666f*(u[2]*u[2]);
        case B: return rho*-0.16666666666666666f*u[2] + rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[2]*u[2]);
        case NW: return rho*-0.083333333333333329f*u[0] + rho*-0.25f*u[0]*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]);
        case NE: return rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.25f*u[0]*u[1];
        case SW: return rho*-0.083333333333333329f*u[0] + rho*-0.083333333333333329f*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.25f*u[0]*u[1];
        case SE: return rho*-0.083333333333333329f*u[1] + rho*-0.25f*u[0]*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]);
        case TN: return rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[1]*u[2];
        case TS: return rho*-0.083333333333333329f*u[1] + rho*-0.25f*u[1]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]);
        case TW: return rho*-0.083333333333333329f*u[0] + rho*-0.25f*u[0]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]);
        case TE: return rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[0]*u[2];
        case BN: return rho*-0.083333333333333329f*u[2] + rho*-0.25f*u[1]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]);
        case BS: return rho*-0.083333333333333329f*u[1] + rho*-0.083333333333333329f*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[1]*u[2];
        case BW: return rho*-0.083333333333333329f*u[0] + rho*-0.083333333333333329f*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[0]*u[2];
        case BE: return rho*-0.083333333333333329f*u[2] + rho*-0.25f*u[0]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]);
        default:
            WALBERLA_ABORT("Invalid Direction")
    }
    
   }
} // EquilibriumDistribution


namespace Equilibrium
{
   template< typename FieldPtrOrIterator, std::enable_if_t<std::is_same<decltype(*(FieldPtrOrIterator())), float>::value, bool> = true >
   void set( FieldPtrOrIterator & it,
             const Vector3< float > & u = Vector3< float >( float(0.0) ),
             float rho = float(1.0) )
   {
        

       it[0] = rho*-0.33333333333333331f*(u[0]*u[0]) + rho*-0.33333333333333331f*(u[1]*u[1]) + rho*-0.33333333333333331f*(u[2]*u[2]) + rho*0.33333333333333331f;
       it[1] = rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[1] + rho*0.16666666666666666f*(u[1]*u[1]);
       it[2] = rho*-0.16666666666666666f*u[1] + rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[1]*u[1]);
       it[3] = rho*-0.16666666666666666f*u[0] + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[0]*u[0]);
       it[4] = rho*-0.16666666666666666f*(u[1]*u[1]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[0] + rho*0.16666666666666666f*(u[0]*u[0]);
       it[5] = rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[2] + rho*0.16666666666666666f*(u[2]*u[2]);
       it[6] = rho*-0.16666666666666666f*u[2] + rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[2]*u[2]);
       it[7] = rho*-0.083333333333333329f*u[0] + rho*-0.25f*u[0]*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]);
       it[8] = rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.25f*u[0]*u[1];
       it[9] = rho*-0.083333333333333329f*u[0] + rho*-0.083333333333333329f*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.25f*u[0]*u[1];
       it[10] = rho*-0.083333333333333329f*u[1] + rho*-0.25f*u[0]*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]);
       it[11] = rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[1]*u[2];
       it[12] = rho*-0.083333333333333329f*u[1] + rho*-0.25f*u[1]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]);
       it[13] = rho*-0.083333333333333329f*u[0] + rho*-0.25f*u[0]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]);
       it[14] = rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[0]*u[2];
       it[15] = rho*-0.083333333333333329f*u[2] + rho*-0.25f*u[1]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]);
       it[16] = rho*-0.083333333333333329f*u[1] + rho*-0.083333333333333329f*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[1]*u[2];
       it[17] = rho*-0.083333333333333329f*u[0] + rho*-0.083333333333333329f*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[0]*u[2];
       it[18] = rho*-0.083333333333333329f*u[2] + rho*-0.25f*u[0]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]);
       }

   template< typename PdfField_T, std::enable_if_t<std::is_same<typename PdfField_T::value_type, float>::value, bool> = true >
   void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
             const Vector3< float > & u = Vector3< float >( float(0.0) ),
             float rho = float(1.0) )
   {
      

      float & xyz0 = pdf(x,y,z,0);
      pdf.getF( &xyz0, 0)= rho*-0.33333333333333331f*(u[0]*u[0]) + rho*-0.33333333333333331f*(u[1]*u[1]) + rho*-0.33333333333333331f*(u[2]*u[2]) + rho*0.33333333333333331f;
      pdf.getF( &xyz0, 1)= rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[1] + rho*0.16666666666666666f*(u[1]*u[1]);
      pdf.getF( &xyz0, 2)= rho*-0.16666666666666666f*u[1] + rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[1]*u[1]);
      pdf.getF( &xyz0, 3)= rho*-0.16666666666666666f*u[0] + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[0]*u[0]);
      pdf.getF( &xyz0, 4)= rho*-0.16666666666666666f*(u[1]*u[1]) + rho*-0.16666666666666666f*(u[2]*u[2]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[0] + rho*0.16666666666666666f*(u[0]*u[0]);
      pdf.getF( &xyz0, 5)= rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*0.055555555555555552f + rho*0.16666666666666666f*u[2] + rho*0.16666666666666666f*(u[2]*u[2]);
      pdf.getF( &xyz0, 6)= rho*-0.16666666666666666f*u[2] + rho*-0.16666666666666666f*(u[0]*u[0]) + rho*-0.16666666666666666f*(u[1]*u[1]) + rho*0.055555555555555552f + rho*0.16666666666666666f*(u[2]*u[2]);
      pdf.getF( &xyz0, 7)= rho*-0.083333333333333329f*u[0] + rho*-0.25f*u[0]*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]);
      pdf.getF( &xyz0, 8)= rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.25f*u[0]*u[1];
      pdf.getF( &xyz0, 9)= rho*-0.083333333333333329f*u[0] + rho*-0.083333333333333329f*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.25f*u[0]*u[1];
      pdf.getF( &xyz0, 10)= rho*-0.083333333333333329f*u[1] + rho*-0.25f*u[0]*u[1] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[1]*u[1]);
      pdf.getF( &xyz0, 11)= rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[1]*u[2];
      pdf.getF( &xyz0, 12)= rho*-0.083333333333333329f*u[1] + rho*-0.25f*u[1]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]);
      pdf.getF( &xyz0, 13)= rho*-0.083333333333333329f*u[0] + rho*-0.25f*u[0]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]);
      pdf.getF( &xyz0, 14)= rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*u[2] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[0]*u[2];
      pdf.getF( &xyz0, 15)= rho*-0.083333333333333329f*u[2] + rho*-0.25f*u[1]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[1] + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]);
      pdf.getF( &xyz0, 16)= rho*-0.083333333333333329f*u[1] + rho*-0.083333333333333329f*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[1]*u[1]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[1]*u[2];
      pdf.getF( &xyz0, 17)= rho*-0.083333333333333329f*u[0] + rho*-0.083333333333333329f*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]) + rho*0.25f*u[0]*u[2];
      pdf.getF( &xyz0, 18)= rho*-0.083333333333333329f*u[2] + rho*-0.25f*u[0]*u[2] + rho*0.027777777777777776f + rho*0.083333333333333329f*u[0] + rho*0.083333333333333329f*(u[0]*u[0]) + rho*0.083333333333333329f*(u[2]*u[2]);
      }
} // Equilibrium


namespace Density
{
   template< typename FieldPtrOrIterator, std::enable_if_t<std::is_same<decltype(*(FieldPtrOrIterator())), float>::value, bool> = true >
   inline float get( const FieldPtrOrIterator & it )
   {
        const float f_0 = it[0];
        const float f_1 = it[1];
        const float f_2 = it[2];
        const float f_3 = it[3];
        const float f_4 = it[4];
        const float f_5 = it[5];
        const float f_6 = it[6];
        const float f_7 = it[7];
        const float f_8 = it[8];
        const float f_9 = it[9];
        const float f_10 = it[10];
        const float f_11 = it[11];
        const float f_12 = it[12];
        const float f_13 = it[13];
        const float f_14 = it[14];
        const float f_15 = it[15];
        const float f_16 = it[16];
        const float f_17 = it[17];
        const float f_18 = it[18];
        const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const float vel1Term = f_1 + f_11 + f_15 + f_7;
        const float vel2Term = f_12 + f_13 + f_5;
        const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        return rho;
   }

   template< typename PdfField_T, std::enable_if_t<std::is_same<typename PdfField_T::value_type, float>::value, bool> = true >
   inline float get( const PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const float & xyz0 = pdf(x,y,z,0);
        const float f_0 = pdf.getF( &xyz0, 0);
        const float f_1 = pdf.getF( &xyz0, 1);
        const float f_2 = pdf.getF( &xyz0, 2);
        const float f_3 = pdf.getF( &xyz0, 3);
        const float f_4 = pdf.getF( &xyz0, 4);
        const float f_5 = pdf.getF( &xyz0, 5);
        const float f_6 = pdf.getF( &xyz0, 6);
        const float f_7 = pdf.getF( &xyz0, 7);
        const float f_8 = pdf.getF( &xyz0, 8);
        const float f_9 = pdf.getF( &xyz0, 9);
        const float f_10 = pdf.getF( &xyz0, 10);
        const float f_11 = pdf.getF( &xyz0, 11);
        const float f_12 = pdf.getF( &xyz0, 12);
        const float f_13 = pdf.getF( &xyz0, 13);
        const float f_14 = pdf.getF( &xyz0, 14);
        const float f_15 = pdf.getF( &xyz0, 15);
        const float f_16 = pdf.getF( &xyz0, 16);
        const float f_17 = pdf.getF( &xyz0, 17);
        const float f_18 = pdf.getF( &xyz0, 18);
        const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const float vel1Term = f_1 + f_11 + f_15 + f_7;
        const float vel2Term = f_12 + f_13 + f_5;
        const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        return rho;
   }
} // Density


namespace DensityAndVelocity
{
    template< typename FieldPtrOrIterator >
    void set( FieldPtrOrIterator & it,
              const GhostLayerField<float, 3u> & force_field,
              const Vector3< float > & u = Vector3< float >( float(0.0) ),
              const float rho_in = float(1.0) )
    {
        auto x = it.x();
        auto y = it.y();
        auto z = it.z();

        const float rho = rho_in;
        const float delta_rho = rho - 1;
        const float u_0 = -force_field.get(x,y,z,0)*0.50000000000000000f/rho + u[0];
        const float u_1 = -force_field.get(x,y,z,1)*0.50000000000000000f/rho + u[1];
        const float u_2 = -force_field.get(x,y,z,2)*0.50000000000000000f/rho + u[2];

        Equilibrium::set(it, Vector3<float>(u_0, u_1, u_2), rho);
    }

    template< typename PdfField_T >
    void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
              const GhostLayerField<float, 3u> & force_field,
              const Vector3< float > & u = Vector3< float >( float(0.0) ),
              const float rho_in = float(1.0) )
    {
        const float rho = rho_in;
        const float delta_rho = rho - 1;
        const float u_0 = -force_field.get(x,y,z,0)*0.50000000000000000f/rho + u[0];
        const float u_1 = -force_field.get(x,y,z,1)*0.50000000000000000f/rho + u[1];
        const float u_2 = -force_field.get(x,y,z,2)*0.50000000000000000f/rho + u[2];

        Equilibrium::set(pdf, x, y, z, Vector3<float>(u_0, u_1, u_2), rho );
    }
} // DensityAndVelocity


namespace DensityAndVelocityRange
{

   template< typename FieldIteratorXYZ >
   void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end,
             const GhostLayerField<float, 3u> & force_field,
             const Vector3< float > & u = Vector3< float >( float(0.0) ),
             const float rho_in = float(1.0) )
   {
        for( auto cellIt = begin; cellIt != end; ++cellIt )
        {
            const auto x = cellIt.x();
            const auto y = cellIt.y();
            const auto z = cellIt.z();
            const float rho = rho_in;
            const float delta_rho = rho - 1;
            const float u_0 = -force_field.get(x,y,z,0)*0.50000000000000000f/rho + u[0];
            const float u_1 = -force_field.get(x,y,z,1)*0.50000000000000000f/rho + u[1];
            const float u_2 = -force_field.get(x,y,z,2)*0.50000000000000000f/rho + u[2];

            Equilibrium::set(cellIt, Vector3<float>(u_0, u_1, u_2), rho);
        }
   }
} // DensityAndVelocityRange



namespace DensityAndMomentumDensity
{
   template< typename FieldPtrOrIterator >
   float get( Vector3< float > & momentumDensity,
               const GhostLayerField<float, 3u> & force_field,
               const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        const float f_0 = it[0];
        const float f_1 = it[1];
        const float f_2 = it[2];
        const float f_3 = it[3];
        const float f_4 = it[4];
        const float f_5 = it[5];
        const float f_6 = it[6];
        const float f_7 = it[7];
        const float f_8 = it[8];
        const float f_9 = it[9];
        const float f_10 = it[10];
        const float f_11 = it[11];
        const float f_12 = it[12];
        const float f_13 = it[13];
        const float f_14 = it[14];
        const float f_15 = it[15];
        const float f_16 = it[16];
        const float f_17 = it[17];
        const float f_18 = it[18];
        const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
        const float vel1Term = f_1 + f_11 + f_15 + f_7;
        const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
        const float vel2Term = f_12 + f_13 + f_5;
        const float momdensity_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
        const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        const float md_0 = force_field.get(x,y,z,0)*0.50000000000000000f + momdensity_0;
        const float md_1 = force_field.get(x,y,z,1)*0.50000000000000000f + momdensity_1;
        const float md_2 = force_field.get(x,y,z,2)*0.50000000000000000f + momdensity_2;
        momentumDensity[0] = md_0;
        momentumDensity[1] = md_1;
        momentumDensity[2] = md_2;
        
        return rho;
   }

   template< typename PdfField_T >
   float get( Vector3< float > & momentumDensity,
               const GhostLayerField<float, 3u> & force_field,
               const PdfField_T & pdf,
               const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const float & xyz0 = pdf(x,y,z,0);
        const float f_0 = pdf.getF( &xyz0, 0);
        const float f_1 = pdf.getF( &xyz0, 1);
        const float f_2 = pdf.getF( &xyz0, 2);
        const float f_3 = pdf.getF( &xyz0, 3);
        const float f_4 = pdf.getF( &xyz0, 4);
        const float f_5 = pdf.getF( &xyz0, 5);
        const float f_6 = pdf.getF( &xyz0, 6);
        const float f_7 = pdf.getF( &xyz0, 7);
        const float f_8 = pdf.getF( &xyz0, 8);
        const float f_9 = pdf.getF( &xyz0, 9);
        const float f_10 = pdf.getF( &xyz0, 10);
        const float f_11 = pdf.getF( &xyz0, 11);
        const float f_12 = pdf.getF( &xyz0, 12);
        const float f_13 = pdf.getF( &xyz0, 13);
        const float f_14 = pdf.getF( &xyz0, 14);
        const float f_15 = pdf.getF( &xyz0, 15);
        const float f_16 = pdf.getF( &xyz0, 16);
        const float f_17 = pdf.getF( &xyz0, 17);
        const float f_18 = pdf.getF( &xyz0, 18);
        const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
        const float vel1Term = f_1 + f_11 + f_15 + f_7;
        const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
        const float vel2Term = f_12 + f_13 + f_5;
        const float momdensity_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
        const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        const float md_0 = force_field.get(x,y,z,0)*0.50000000000000000f + momdensity_0;
        const float md_1 = force_field.get(x,y,z,1)*0.50000000000000000f + momdensity_1;
        const float md_2 = force_field.get(x,y,z,2)*0.50000000000000000f + momdensity_2;
        momentumDensity[0] = md_0;
        momentumDensity[1] = md_1;
        momentumDensity[2] = md_2;
        
       return rho;
   }
} // DensityAndMomentumDensity


namespace MomentumDensity
{
   template< typename FieldPtrOrIterator >
   void get( Vector3< float > & momentumDensity,
             const GhostLayerField<float, 3u> & force_field,
             const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        const float f_0 = it[0];
        const float f_1 = it[1];
        const float f_2 = it[2];
        const float f_3 = it[3];
        const float f_4 = it[4];
        const float f_5 = it[5];
        const float f_6 = it[6];
        const float f_7 = it[7];
        const float f_8 = it[8];
        const float f_9 = it[9];
        const float f_10 = it[10];
        const float f_11 = it[11];
        const float f_12 = it[12];
        const float f_13 = it[13];
        const float f_14 = it[14];
        const float f_15 = it[15];
        const float f_16 = it[16];
        const float f_17 = it[17];
        const float f_18 = it[18];
        const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
        const float vel1Term = f_1 + f_11 + f_15 + f_7;
        const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
        const float vel2Term = f_12 + f_13 + f_5;
        const float momdensity_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
        const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        const float md_0 = force_field.get(x,y,z,0)*0.50000000000000000f + momdensity_0;
        const float md_1 = force_field.get(x,y,z,1)*0.50000000000000000f + momdensity_1;
        const float md_2 = force_field.get(x,y,z,2)*0.50000000000000000f + momdensity_2;
        momentumDensity[0] = md_0;
        momentumDensity[1] = md_1;
        momentumDensity[2] = md_2;
        
   }

   template< typename PdfField_T >
   void get( Vector3< float > & momentumDensity,
             const GhostLayerField<float, 3u> & force_field,
             const PdfField_T & pdf,
             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const float & xyz0 = pdf(x,y,z,0);
        const float f_0 = pdf.getF( &xyz0, 0);
        const float f_1 = pdf.getF( &xyz0, 1);
        const float f_2 = pdf.getF( &xyz0, 2);
        const float f_3 = pdf.getF( &xyz0, 3);
        const float f_4 = pdf.getF( &xyz0, 4);
        const float f_5 = pdf.getF( &xyz0, 5);
        const float f_6 = pdf.getF( &xyz0, 6);
        const float f_7 = pdf.getF( &xyz0, 7);
        const float f_8 = pdf.getF( &xyz0, 8);
        const float f_9 = pdf.getF( &xyz0, 9);
        const float f_10 = pdf.getF( &xyz0, 10);
        const float f_11 = pdf.getF( &xyz0, 11);
        const float f_12 = pdf.getF( &xyz0, 12);
        const float f_13 = pdf.getF( &xyz0, 13);
        const float f_14 = pdf.getF( &xyz0, 14);
        const float f_15 = pdf.getF( &xyz0, 15);
        const float f_16 = pdf.getF( &xyz0, 16);
        const float f_17 = pdf.getF( &xyz0, 17);
        const float f_18 = pdf.getF( &xyz0, 18);
        const float vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const float momdensity_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + vel0Term;
        const float vel1Term = f_1 + f_11 + f_15 + f_7;
        const float momdensity_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + vel1Term;
        const float vel2Term = f_12 + f_13 + f_5;
        const float momdensity_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + vel2Term;
        const float rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        const float md_0 = force_field.get(x,y,z,0)*0.50000000000000000f + momdensity_0;
        const float md_1 = force_field.get(x,y,z,1)*0.50000000000000000f + momdensity_1;
        const float md_2 = force_field.get(x,y,z,2)*0.50000000000000000f + momdensity_2;
        momentumDensity[0] = md_0;
        momentumDensity[1] = md_1;
        momentumDensity[2] = md_2;
        
   }
} // MomentumDensity


namespace PressureTensor
{
   template< typename FieldPtrOrIterator, std::enable_if_t<std::is_same<decltype(*(FieldPtrOrIterator())), float>::value, bool> = true >
   void get( Matrix3< float > & pressureTensor, const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        const float f_0 = it[0];
        const float f_1 = it[1];
        const float f_2 = it[2];
        const float f_3 = it[3];
        const float f_4 = it[4];
        const float f_5 = it[5];
        const float f_6 = it[6];
        const float f_7 = it[7];
        const float f_8 = it[8];
        const float f_9 = it[9];
        const float f_10 = it[10];
        const float f_11 = it[11];
        const float f_12 = it[12];
        const float f_13 = it[13];
        const float f_14 = it[14];
        const float f_15 = it[15];
        const float f_16 = it[16];
        const float f_17 = it[17];
        const float f_18 = it[18];
        const float p_0 = f_10 + f_13 + f_14 + f_17 + f_18 + f_3 + f_4 + f_7 + f_8 + f_9;
        const float p_1 = -f_10 - f_7 + f_8 + f_9;
        const float p_2 = -f_13 + f_14 + f_17 - f_18;
        const float p_3 = -f_10 - f_7 + f_8 + f_9;
        const float p_4 = f_1 + f_10 + f_11 + f_12 + f_15 + f_16 + f_2 + f_7 + f_8 + f_9;
        const float p_5 = f_11 - f_12 - f_15 + f_16;
        const float p_6 = -f_13 + f_14 + f_17 - f_18;
        const float p_7 = f_11 - f_12 - f_15 + f_16;
        const float p_8 = f_11 + f_12 + f_13 + f_14 + f_15 + f_16 + f_17 + f_18 + f_5 + f_6;
        pressureTensor[0] = p_0;
            pressureTensor[1] = p_1;
            pressureTensor[2] = p_2;
            
        pressureTensor[3] = p_3;
            pressureTensor[4] = p_4;
            pressureTensor[5] = p_5;
            
        pressureTensor[6] = p_6;
            pressureTensor[7] = p_7;
            pressureTensor[8] = p_8;
            
        
   }

   template< typename PdfField_T, std::enable_if_t<std::is_same<typename PdfField_T::value_type, float>::value, bool> = true >
   void get( Matrix3< float > & pressureTensor, const PdfField_T & pdf,
             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const float & xyz0 = pdf(x,y,z,0);
        const float f_0 = pdf.getF( &xyz0, 0);
        const float f_1 = pdf.getF( &xyz0, 1);
        const float f_2 = pdf.getF( &xyz0, 2);
        const float f_3 = pdf.getF( &xyz0, 3);
        const float f_4 = pdf.getF( &xyz0, 4);
        const float f_5 = pdf.getF( &xyz0, 5);
        const float f_6 = pdf.getF( &xyz0, 6);
        const float f_7 = pdf.getF( &xyz0, 7);
        const float f_8 = pdf.getF( &xyz0, 8);
        const float f_9 = pdf.getF( &xyz0, 9);
        const float f_10 = pdf.getF( &xyz0, 10);
        const float f_11 = pdf.getF( &xyz0, 11);
        const float f_12 = pdf.getF( &xyz0, 12);
        const float f_13 = pdf.getF( &xyz0, 13);
        const float f_14 = pdf.getF( &xyz0, 14);
        const float f_15 = pdf.getF( &xyz0, 15);
        const float f_16 = pdf.getF( &xyz0, 16);
        const float f_17 = pdf.getF( &xyz0, 17);
        const float f_18 = pdf.getF( &xyz0, 18);
        const float p_0 = f_10 + f_13 + f_14 + f_17 + f_18 + f_3 + f_4 + f_7 + f_8 + f_9;
        const float p_1 = -f_10 - f_7 + f_8 + f_9;
        const float p_2 = -f_13 + f_14 + f_17 - f_18;
        const float p_3 = -f_10 - f_7 + f_8 + f_9;
        const float p_4 = f_1 + f_10 + f_11 + f_12 + f_15 + f_16 + f_2 + f_7 + f_8 + f_9;
        const float p_5 = f_11 - f_12 - f_15 + f_16;
        const float p_6 = -f_13 + f_14 + f_17 - f_18;
        const float p_7 = f_11 - f_12 - f_15 + f_16;
        const float p_8 = f_11 + f_12 + f_13 + f_14 + f_15 + f_16 + f_17 + f_18 + f_5 + f_6;
        pressureTensor[0] = p_0;
            pressureTensor[1] = p_1;
            pressureTensor[2] = p_2;
            
        pressureTensor[3] = p_3;
            pressureTensor[4] = p_4;
            pressureTensor[5] = p_5;
            
        pressureTensor[6] = p_6;
            pressureTensor[7] = p_7;
            pressureTensor[8] = p_8;
            
        
   }
} // PressureTensor


namespace ShearRate
{
   template< typename FieldPtrOrIterator, std::enable_if_t<std::is_same<decltype(*(FieldPtrOrIterator())), float>::value, bool> = true >
   inline float get( const FieldPtrOrIterator & /* it */,
                          const Vector3< float > & /* velocity */,
                          const float /* rho */)
   {
       WALBERLA_ABORT("Not implemented");
       return float(0.0);
   }

   template< typename PdfField_T, std::enable_if_t<std::is_same<typename PdfField_T::value_type, float>::value, bool> = true >
   inline float get( const PdfField_T & /* pdf */,
                          const cell_idx_t /* x */, const cell_idx_t /* y */, const cell_idx_t /* z */,
                          const Vector3< float > & /* velocity */, const float /* rho */ )
   {
       WALBERLA_ABORT("Not implemented");
       return float(0.0);
   }
} // ShearRate


} // namespace accessor
} // namespace lbm
} // namespace walberla



#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif
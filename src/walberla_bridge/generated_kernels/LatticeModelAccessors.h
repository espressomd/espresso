// kernel generated with pystencils v0.3.3+39.g587a822, lbmpy v0.3.3+33.g036fe13, lbmpy_walberla/pystencils_walberla from commit ref: refs/heads/LeesEdwards

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


#include "core/DataTypes.h"
#include "core/logging/Logging.h"

#include "lbm/lattice_model/EquilibriumDistribution.h"
#include "lbm/field/Density.h"
#include "lbm/field/DensityAndMomentumDensity.h"
#include "lbm/field/DensityAndVelocity.h"
#include "lbm/field/PressureTensor.h"
#include "lbm/field/ShearRate.h"

#include <vector>

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



//======================================================================================================================
//
//  Implementation of macroscopic value backend
//
//======================================================================================================================


template <typename CollisionModel, typename FloatType>
class EquilibriumDistribution< LBWalberlaImpl<CollisionModel, FloatType>, void>
{
public:
   typedef typename LBWalberlaImpl<CollisionModel, FloatType>::Stencil Stencil;

   static real_t get( const stencil::Direction direction,
                      const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ),
                      real_t rho = real_t(1.0) )
   {
        
        
    using namespace stencil;
    switch( direction ) {
        case C: return rho*-0.5*(u[0]*u[0]) + rho*-0.5*(u[1]*u[1]) + rho*-0.5*(u[2]*u[2]) + rho*0.333333333333333;
        case N: return rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[1]*u[1]) + rho*0.166666666666667*u[1];
        case S: return rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*-0.166666666666667*u[1] + rho*0.0555555555555556 + rho*0.166666666666667*(u[1]*u[1]);
        case W: return rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*-0.166666666666667*u[0] + rho*0.0555555555555556 + rho*0.166666666666667*(u[0]*u[0]);
        case E: return rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[0]*u[0]) + rho*0.166666666666667*u[0];
        case T: return rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[1]*u[1]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[2]*u[2]) + rho*0.166666666666667*u[2];
        case B: return rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.166666666666667*u[2] + rho*0.0555555555555556 + rho*0.166666666666667*(u[2]*u[2]);
        case NW: return rho*-0.0833333333333333*u[0] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*-0.25*u[0]*u[1] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0277777777777778;
        case NE: return rho*-0.0416666666666667*(u[2]*u[2]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0277777777777778 + rho*0.25*u[0]*u[1];
        case SW: return rho*-0.0833333333333333*u[0] + rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[1];
        case SE: return rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*-0.25*u[0]*u[1] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778;
        case TN: return rho*-0.0416666666666667*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778 + rho*0.25*u[1]*u[2];
        case TS: return rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*-0.25*u[1]*u[2] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778;
        case TW: return rho*-0.0833333333333333*u[0] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*-0.25*u[0]*u[2] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778;
        case TE: return rho*-0.0416666666666667*(u[1]*u[1]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778 + rho*0.25*u[0]*u[2];
        case BN: return rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*-0.25*u[1]*u[2] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
        case BS: return rho*-0.0833333333333333*u[1] + rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[1]*u[2];
        case BW: return rho*-0.0833333333333333*u[0] + rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[2];
        case BE: return rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*-0.25*u[0]*u[2] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
        default:
            WALBERLA_ABORT("Invalid Direction");
    }
    
   }

   static real_t getSymmetricPart( const stencil::Direction direction,
                                   const Vector3<real_t> & u = Vector3< real_t >(real_t(0.0)),
                                   real_t rho = real_t(1.0) )
   {
        
        
    using namespace stencil;
    switch( direction ) {
        case C: return rho*((double)(-0.5))*(u[0]*u[0]) + rho*((double)(-0.5))*(u[1]*u[1]) + rho*((double)(-0.5))*(u[2]*u[2]) + rho*0.333333333333333;
        case N: return rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[1]*u[1]);
        case S: return rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[1]*u[1]);
        case W: return rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[0]*u[0]);
        case E: return rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[0]*u[0]);
        case T: return rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[1]*u[1]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[2]*u[2]);
        case B: return rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[1]*u[1]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[2]*u[2]);
        case NW: return rho*-0.0416666666666667*(u[2]*u[2]) + rho*-0.25*u[0]*u[1] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778;
        case NE: return rho*-0.0416666666666667*(u[2]*u[2]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[1];
        case SW: return rho*-0.0416666666666667*(u[2]*u[2]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[1];
        case SE: return rho*-0.0416666666666667*(u[2]*u[2]) + rho*-0.25*u[0]*u[1] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778;
        case TN: return rho*-0.0416666666666667*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[1]*u[2];
        case TS: return rho*-0.0416666666666667*(u[0]*u[0]) + rho*-0.25*u[1]*u[2] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
        case TW: return rho*-0.0416666666666667*(u[1]*u[1]) + rho*-0.25*u[0]*u[2] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
        case TE: return rho*-0.0416666666666667*(u[1]*u[1]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[2];
        case BN: return rho*-0.0416666666666667*(u[0]*u[0]) + rho*-0.25*u[1]*u[2] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
        case BS: return rho*-0.0416666666666667*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[1]*u[2];
        case BW: return rho*-0.0416666666666667*(u[1]*u[1]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[2];
        case BE: return rho*-0.0416666666666667*(u[1]*u[1]) + rho*-0.25*u[0]*u[2] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
        default:
            WALBERLA_ABORT("Invalid Direction");
    }
    
   }

   static real_t getAsymmetricPart( const stencil::Direction direction,
                                    const Vector3< real_t > & u = Vector3<real_t>( real_t(0.0) ),
                                    real_t rho = real_t(1.0) )
   {
        
        
    using namespace stencil;
    switch( direction ) {
        case C: return rho*-0.5*(u[0]*u[0]) + rho*-0.5*(u[1]*u[1]) + rho*-0.5*(u[2]*u[2]) - rho*((double)(-0.5))*(u[0]*u[0]) - rho*((double)(-0.5))*(u[1]*u[1]) - rho*((double)(-0.5))*(u[2]*u[2]);
        case N: return rho*0.166666666666667*u[1];
        case S: return rho*-0.166666666666667*u[1];
        case W: return rho*-0.166666666666667*u[0];
        case E: return rho*0.166666666666667*u[0];
        case T: return rho*0.166666666666667*u[2];
        case B: return rho*-0.166666666666667*u[2];
        case NW: return rho*-0.0833333333333333*u[0] + rho*0.0833333333333333*u[1];
        case NE: return rho*0.0833333333333333*u[0] + rho*0.0833333333333333*u[1];
        case SW: return rho*-0.0833333333333333*u[0] + rho*-0.0833333333333333*u[1];
        case SE: return rho*-0.0833333333333333*u[1] + rho*0.0833333333333333*u[0];
        case TN: return rho*0.0833333333333333*u[1] + rho*0.0833333333333333*u[2];
        case TS: return rho*-0.0833333333333333*u[1] + rho*0.0833333333333333*u[2];
        case TW: return rho*-0.0833333333333333*u[0] + rho*0.0833333333333333*u[2];
        case TE: return rho*0.0833333333333333*u[0] + rho*0.0833333333333333*u[2];
        case BN: return rho*-0.0833333333333333*u[2] + rho*0.0833333333333333*u[1];
        case BS: return rho*-0.0833333333333333*u[1] + rho*-0.0833333333333333*u[2];
        case BW: return rho*-0.0833333333333333*u[0] + rho*-0.0833333333333333*u[2];
        case BE: return rho*-0.0833333333333333*u[2] + rho*0.0833333333333333*u[0];
        default:
            WALBERLA_ABORT("Invalid Direction");
    }
    
   }

   static std::vector< real_t > get( const Vector3< real_t > & u = Vector3<real_t>( real_t(0.0) ),
                                     real_t rho = real_t(1.0) )
   {
      

      std::vector< real_t > equilibrium( Stencil::Size );
      for( auto d = Stencil::begin(); d != Stencil::end(); ++d )
      {
         equilibrium[d.toIdx()] = get(*d, u, rho);
      }
      return equilibrium;
   }
};


namespace internal {

template<typename CollisionModel, typename FloatType>
struct AdaptVelocityToForce<LBWalberlaImpl<CollisionModel, FloatType>, void>
{
   template< typename FieldPtrOrIterator >
   static Vector3<real_t> get( FieldPtrOrIterator & it, const LBWalberlaImpl<CollisionModel, FloatType> & lm,
                               const Vector3< real_t > & velocity, const real_t rho )
   {
      auto x = it.x();
      auto y = it.y();
      auto z = it.z();
      
      return velocity - Vector3<real_t>(lm.force_->get(x,y,z,0)*0.5/rho,lm.force_->get(x,y,z,1)*0.5/rho,lm.force_->get(x,y,z,2)*0.5/rho  );
      
   }

   static Vector3<real_t> get( const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const LBWalberlaImpl<CollisionModel, FloatType> & lm,
                               const Vector3< real_t > & velocity, const real_t rho )
   {
      

      return velocity - Vector3<real_t>(lm.force_->get(x,y,z,0)*0.5/rho,lm.force_->get(x,y,z,1)*0.5/rho,lm.force_->get(x,y,z,2)*0.5/rho  );
      
   }
};
} // namespace internal




template <typename CollisionModel, typename FloatType>
struct Equilibrium< LBWalberlaImpl<CollisionModel, FloatType>, void >
{

   template< typename FieldPtrOrIterator >
   static void set( FieldPtrOrIterator & it,
                    const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), real_t rho = real_t(1.0) )
   {
        

       it[0] = rho*-0.5*(u[0]*u[0]) + rho*-0.5*(u[1]*u[1]) + rho*-0.5*(u[2]*u[2]) + rho*0.333333333333333;
       it[1] = rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[1]*u[1]) + rho*0.166666666666667*u[1];
       it[2] = rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*-0.166666666666667*u[1] + rho*0.0555555555555556 + rho*0.166666666666667*(u[1]*u[1]);
       it[3] = rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*-0.166666666666667*u[0] + rho*0.0555555555555556 + rho*0.166666666666667*(u[0]*u[0]);
       it[4] = rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[0]*u[0]) + rho*0.166666666666667*u[0];
       it[5] = rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[1]*u[1]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[2]*u[2]) + rho*0.166666666666667*u[2];
       it[6] = rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.166666666666667*u[2] + rho*0.0555555555555556 + rho*0.166666666666667*(u[2]*u[2]);
       it[7] = rho*-0.0833333333333333*u[0] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*-0.25*u[0]*u[1] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0277777777777778;
       it[8] = rho*-0.0416666666666667*(u[2]*u[2]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0277777777777778 + rho*0.25*u[0]*u[1];
       it[9] = rho*-0.0833333333333333*u[0] + rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[1];
       it[10] = rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*-0.25*u[0]*u[1] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778;
       it[11] = rho*-0.0416666666666667*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778 + rho*0.25*u[1]*u[2];
       it[12] = rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*-0.25*u[1]*u[2] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778;
       it[13] = rho*-0.0833333333333333*u[0] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*-0.25*u[0]*u[2] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778;
       it[14] = rho*-0.0416666666666667*(u[1]*u[1]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778 + rho*0.25*u[0]*u[2];
       it[15] = rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*-0.25*u[1]*u[2] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
       it[16] = rho*-0.0833333333333333*u[1] + rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[1]*u[2];
       it[17] = rho*-0.0833333333333333*u[0] + rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[2];
       it[18] = rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*-0.25*u[0]*u[2] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
       }

   template< typename PdfField_T >
   static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
                    const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), real_t rho = real_t(1.0) )
   {
      

      real_t & xyz0 = pdf(x,y,z,0);
      pdf.getF( &xyz0, 0)= rho*-0.5*(u[0]*u[0]) + rho*-0.5*(u[1]*u[1]) + rho*-0.5*(u[2]*u[2]) + rho*0.333333333333333;
      pdf.getF( &xyz0, 1)= rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[1]*u[1]) + rho*0.166666666666667*u[1];
      pdf.getF( &xyz0, 2)= rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*-0.166666666666667*u[1] + rho*0.0555555555555556 + rho*0.166666666666667*(u[1]*u[1]);
      pdf.getF( &xyz0, 3)= rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*-0.166666666666667*u[0] + rho*0.0555555555555556 + rho*0.166666666666667*(u[0]*u[0]);
      pdf.getF( &xyz0, 4)= rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.0833333333333333*(u[2]*u[2]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[0]*u[0]) + rho*0.166666666666667*u[0];
      pdf.getF( &xyz0, 5)= rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[1]*u[1]) + rho*0.0555555555555556 + rho*0.166666666666667*(u[2]*u[2]) + rho*0.166666666666667*u[2];
      pdf.getF( &xyz0, 6)= rho*-0.0833333333333333*(u[0]*u[0]) + rho*-0.0833333333333333*(u[1]*u[1]) + rho*-0.166666666666667*u[2] + rho*0.0555555555555556 + rho*0.166666666666667*(u[2]*u[2]);
      pdf.getF( &xyz0, 7)= rho*-0.0833333333333333*u[0] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*-0.25*u[0]*u[1] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0277777777777778;
      pdf.getF( &xyz0, 8)= rho*-0.0416666666666667*(u[2]*u[2]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0277777777777778 + rho*0.25*u[0]*u[1];
      pdf.getF( &xyz0, 9)= rho*-0.0833333333333333*u[0] + rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[1];
      pdf.getF( &xyz0, 10)= rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[2]*u[2]) + rho*-0.25*u[0]*u[1] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0277777777777778;
      pdf.getF( &xyz0, 11)= rho*-0.0416666666666667*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778 + rho*0.25*u[1]*u[2];
      pdf.getF( &xyz0, 12)= rho*-0.0833333333333333*u[1] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*-0.25*u[1]*u[2] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778;
      pdf.getF( &xyz0, 13)= rho*-0.0833333333333333*u[0] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*-0.25*u[0]*u[2] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778;
      pdf.getF( &xyz0, 14)= rho*-0.0416666666666667*(u[1]*u[1]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0833333333333333*u[2] + rho*0.0277777777777778 + rho*0.25*u[0]*u[2];
      pdf.getF( &xyz0, 15)= rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*-0.25*u[1]*u[2] + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*u[1] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
      pdf.getF( &xyz0, 16)= rho*-0.0833333333333333*u[1] + rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[0]*u[0]) + rho*0.0833333333333333*(u[1]*u[1]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[1]*u[2];
      pdf.getF( &xyz0, 17)= rho*-0.0833333333333333*u[0] + rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778 + rho*0.25*u[0]*u[2];
      pdf.getF( &xyz0, 18)= rho*-0.0833333333333333*u[2] + rho*-0.0416666666666667*(u[1]*u[1]) + rho*-0.25*u[0]*u[2] + rho*0.0833333333333333*(u[0]*u[0]) + rho*0.0833333333333333*u[0] + rho*0.0833333333333333*(u[2]*u[2]) + rho*0.0277777777777778;
      }
};


template <typename CollisionModel, typename FloatType>
struct Density<LBWalberlaImpl<CollisionModel, FloatType>, void>
{
   template< typename FieldPtrOrIterator >
   static inline real_t get( const LBWalberlaImpl<CollisionModel, FloatType> & , const FieldPtrOrIterator & it )
   {
        const real_t f_0 = it[0];
        const real_t f_1 = it[1];
        const real_t f_2 = it[2];
        const real_t f_3 = it[3];
        const real_t f_4 = it[4];
        const real_t f_5 = it[5];
        const real_t f_6 = it[6];
        const real_t f_7 = it[7];
        const real_t f_8 = it[8];
        const real_t f_9 = it[9];
        const real_t f_10 = it[10];
        const real_t f_11 = it[11];
        const real_t f_12 = it[12];
        const real_t f_13 = it[13];
        const real_t f_14 = it[14];
        const real_t f_15 = it[15];
        const real_t f_16 = it[16];
        const real_t f_17 = it[17];
        const real_t f_18 = it[18];
        const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const double vel1Term = f_1 + f_11 + f_15 + f_7;
        const double vel2Term = f_12 + f_13 + f_5;
        const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        return rho;
   }

   template< typename PdfField_T >
   static inline real_t get( const LBWalberlaImpl<CollisionModel, FloatType> & ,
                             const PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const real_t & xyz0 = pdf(x,y,z,0);
        const real_t f_0 = pdf.getF( &xyz0, 0);
        const real_t f_1 = pdf.getF( &xyz0, 1);
        const real_t f_2 = pdf.getF( &xyz0, 2);
        const real_t f_3 = pdf.getF( &xyz0, 3);
        const real_t f_4 = pdf.getF( &xyz0, 4);
        const real_t f_5 = pdf.getF( &xyz0, 5);
        const real_t f_6 = pdf.getF( &xyz0, 6);
        const real_t f_7 = pdf.getF( &xyz0, 7);
        const real_t f_8 = pdf.getF( &xyz0, 8);
        const real_t f_9 = pdf.getF( &xyz0, 9);
        const real_t f_10 = pdf.getF( &xyz0, 10);
        const real_t f_11 = pdf.getF( &xyz0, 11);
        const real_t f_12 = pdf.getF( &xyz0, 12);
        const real_t f_13 = pdf.getF( &xyz0, 13);
        const real_t f_14 = pdf.getF( &xyz0, 14);
        const real_t f_15 = pdf.getF( &xyz0, 15);
        const real_t f_16 = pdf.getF( &xyz0, 16);
        const real_t f_17 = pdf.getF( &xyz0, 17);
        const real_t f_18 = pdf.getF( &xyz0, 18);
        const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const double vel1Term = f_1 + f_11 + f_15 + f_7;
        const double vel2Term = f_12 + f_13 + f_5;
        const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        return rho;
   }
};


template <typename CollisionModel, typename FloatType>
struct MyDensityAndVelocity
{
  template< typename FieldPtrOrIterator, typename VectorField_T >
    static void set( FieldPtrOrIterator & it, const VectorField_T & force,
                     const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), const real_t rho_in = real_t(1.0) )
    {
        auto x = it.x();
        auto y = it.y();
        auto z = it.z();

        const double rho = rho_in;
        const double u_0 = -0.5*force.get(x,y,z,0)/rho_in + u[0];
        const double u_1 = -0.5*force.get(x,y,z,1)/rho_in + u[1];
        const double u_2 = -0.5*force.get(x,y,z,2)/rho_in + u[2];
        

        Equilibrium<LBWalberlaImpl<CollisionModel, FloatType>>::set(it, Vector3<real_t>(u_0, u_1, u_2), rho);
    }

    template< typename PdfField_T, typename VectorField_T >
    static void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z, const VectorField_T & force,
                     const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), const real_t rho_in = real_t(1.0) )
    {
        const double rho = rho_in;
        const double u_0 = -0.5*force.get(x,y,z,0)/rho_in + u[0];
        const double u_1 = -0.5*force.get(x,y,z,1)/rho_in + u[1];
        const double u_2 = -0.5*force.get(x,y,z,2)/rho_in + u[2];
        

        Equilibrium<LBWalberlaImpl<CollisionModel, FloatType>>::set(pdf, x, y, z, Vector3<real_t>(u_0, u_1, u_2), rho );
    }
};


template <typename CollisionModel, typename FloatType, typename FieldIteratorXYZ>
struct DensityAndVelocityRange<LBWalberlaImpl<CollisionModel, FloatType>, FieldIteratorXYZ>
{

   static void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end, const LBWalberlaImpl<CollisionModel, FloatType> & lm,
                    const Vector3< real_t > & u = Vector3< real_t >( real_t(0.0) ), const real_t rho_in = real_t(1.0) )
   {
        for( auto cellIt = begin; cellIt != end; ++cellIt )
        {
            const auto x = cellIt.x();
            const auto y = cellIt.y();
            const auto z = cellIt.z();
            const double rho = rho_in;
            const double u_0 = -0.5*lm.force_->get(x,y,z,0)/rho_in + u[0];
            const double u_1 = -0.5*lm.force_->get(x,y,z,1)/rho_in + u[1];
            const double u_2 = -0.5*lm.force_->get(x,y,z,2)/rho_in + u[2];
            

            Equilibrium<LBWalberlaImpl<CollisionModel, FloatType>>::set(cellIt, Vector3<real_t>(u_0, u_1, u_2), rho);
        }
   }
};



template <typename CollisionModel, typename FloatType>
struct MyDensityAndMomentumDensity
{
  template< typename FieldPtrOrIterator, typename VectorField_T >
   static real_t get( Vector3< real_t > & momentumDensity, const VectorField_T & force,
                      const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        const real_t f_0 = it[0];
        const real_t f_1 = it[1];
        const real_t f_2 = it[2];
        const real_t f_3 = it[3];
        const real_t f_4 = it[4];
        const real_t f_5 = it[5];
        const real_t f_6 = it[6];
        const real_t f_7 = it[7];
        const real_t f_8 = it[8];
        const real_t f_9 = it[9];
        const real_t f_10 = it[10];
        const real_t f_11 = it[11];
        const real_t f_12 = it[12];
        const real_t f_13 = it[13];
        const real_t f_14 = it[14];
        const real_t f_15 = it[15];
        const real_t f_16 = it[16];
        const real_t f_17 = it[17];
        const real_t f_18 = it[18];
        const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const double vel1Term = f_1 + f_11 + f_15 + f_7;
        const double vel2Term = f_12 + f_13 + f_5;
        const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        const double md_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + (0.5)*force.get(x,y,z,0) + vel0Term;
        const double md_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + (0.5)*force.get(x,y,z,1) + vel1Term;
        const double md_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + (0.5)*force.get(x,y,z,2) + vel2Term;
        momentumDensity[0] = md_0;
        momentumDensity[1] = md_1;
        momentumDensity[2] = md_2;
        
        return rho;
   }

   template< typename PdfField_T, typename VectorField_T >
   static real_t get( Vector3< real_t > & momentumDensity, const VectorField_T & force, const PdfField_T & pdf,
                      const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const real_t & xyz0 = pdf(x,y,z,0);
        const real_t f_0 = pdf.getF( &xyz0, 0);
        const real_t f_1 = pdf.getF( &xyz0, 1);
        const real_t f_2 = pdf.getF( &xyz0, 2);
        const real_t f_3 = pdf.getF( &xyz0, 3);
        const real_t f_4 = pdf.getF( &xyz0, 4);
        const real_t f_5 = pdf.getF( &xyz0, 5);
        const real_t f_6 = pdf.getF( &xyz0, 6);
        const real_t f_7 = pdf.getF( &xyz0, 7);
        const real_t f_8 = pdf.getF( &xyz0, 8);
        const real_t f_9 = pdf.getF( &xyz0, 9);
        const real_t f_10 = pdf.getF( &xyz0, 10);
        const real_t f_11 = pdf.getF( &xyz0, 11);
        const real_t f_12 = pdf.getF( &xyz0, 12);
        const real_t f_13 = pdf.getF( &xyz0, 13);
        const real_t f_14 = pdf.getF( &xyz0, 14);
        const real_t f_15 = pdf.getF( &xyz0, 15);
        const real_t f_16 = pdf.getF( &xyz0, 16);
        const real_t f_17 = pdf.getF( &xyz0, 17);
        const real_t f_18 = pdf.getF( &xyz0, 18);
        const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const double vel1Term = f_1 + f_11 + f_15 + f_7;
        const double vel2Term = f_12 + f_13 + f_5;
        const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        const double md_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + (0.5)*force.get(x,y,z,0) + vel0Term;
        const double md_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + (0.5)*force.get(x,y,z,1) + vel1Term;
        const double md_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + (0.5)*force.get(x,y,z,2) + vel2Term;
        momentumDensity[0] = md_0;
        momentumDensity[1] = md_1;
        momentumDensity[2] = md_2;
        
       return rho;
   }
};


template <typename CollisionModel, typename FloatType>
struct MomentumDensity< LBWalberlaImpl<CollisionModel, FloatType>>
{
   template< typename FieldPtrOrIterator >
   static void get( Vector3< real_t > & momentumDensity, const LBWalberlaImpl<CollisionModel, FloatType> & lm, const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        const real_t f_0 = it[0];
        const real_t f_1 = it[1];
        const real_t f_2 = it[2];
        const real_t f_3 = it[3];
        const real_t f_4 = it[4];
        const real_t f_5 = it[5];
        const real_t f_6 = it[6];
        const real_t f_7 = it[7];
        const real_t f_8 = it[8];
        const real_t f_9 = it[9];
        const real_t f_10 = it[10];
        const real_t f_11 = it[11];
        const real_t f_12 = it[12];
        const real_t f_13 = it[13];
        const real_t f_14 = it[14];
        const real_t f_15 = it[15];
        const real_t f_16 = it[16];
        const real_t f_17 = it[17];
        const real_t f_18 = it[18];
        const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const double vel1Term = f_1 + f_11 + f_15 + f_7;
        const double vel2Term = f_12 + f_13 + f_5;
        const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        const double md_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + (0.5)*lm.force_->get(x,y,z,0) + vel0Term;
        const double md_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + (0.5)*lm.force_->get(x,y,z,1) + vel1Term;
        const double md_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + (0.5)*lm.force_->get(x,y,z,2) + vel2Term;
        momentumDensity[0] = md_0;
        momentumDensity[1] = md_1;
        momentumDensity[2] = md_2;
        
   }

   template< typename PdfField_T >
   static void get( Vector3< real_t > & momentumDensity, const LBWalberlaImpl<CollisionModel, FloatType> & lm, const PdfField_T & pdf,
                    const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const real_t & xyz0 = pdf(x,y,z,0);
        const real_t f_0 = pdf.getF( &xyz0, 0);
        const real_t f_1 = pdf.getF( &xyz0, 1);
        const real_t f_2 = pdf.getF( &xyz0, 2);
        const real_t f_3 = pdf.getF( &xyz0, 3);
        const real_t f_4 = pdf.getF( &xyz0, 4);
        const real_t f_5 = pdf.getF( &xyz0, 5);
        const real_t f_6 = pdf.getF( &xyz0, 6);
        const real_t f_7 = pdf.getF( &xyz0, 7);
        const real_t f_8 = pdf.getF( &xyz0, 8);
        const real_t f_9 = pdf.getF( &xyz0, 9);
        const real_t f_10 = pdf.getF( &xyz0, 10);
        const real_t f_11 = pdf.getF( &xyz0, 11);
        const real_t f_12 = pdf.getF( &xyz0, 12);
        const real_t f_13 = pdf.getF( &xyz0, 13);
        const real_t f_14 = pdf.getF( &xyz0, 14);
        const real_t f_15 = pdf.getF( &xyz0, 15);
        const real_t f_16 = pdf.getF( &xyz0, 16);
        const real_t f_17 = pdf.getF( &xyz0, 17);
        const real_t f_18 = pdf.getF( &xyz0, 18);
        const double vel0Term = f_10 + f_14 + f_18 + f_4 + f_8;
        const double vel1Term = f_1 + f_11 + f_15 + f_7;
        const double vel2Term = f_12 + f_13 + f_5;
        const double rho = f_0 + f_16 + f_17 + f_2 + f_3 + f_6 + f_9 + vel0Term + vel1Term + vel2Term;
        const double md_0 = -f_13 - f_17 - f_3 - f_7 - f_9 + (0.5)*lm.force_->get(x,y,z,0) + vel0Term;
        const double md_1 = -f_10 - f_12 - f_16 - f_2 + f_8 - f_9 + (0.5)*lm.force_->get(x,y,z,1) + vel1Term;
        const double md_2 = f_11 + f_14 - f_15 - f_16 - f_17 - f_18 - f_6 + (0.5)*lm.force_->get(x,y,z,2) + vel2Term;
        momentumDensity[0] = md_0;
        momentumDensity[1] = md_1;
        momentumDensity[2] = md_2;
        
   }
};


template <typename CollisionModel, typename FloatType>
struct PressureTensor<LBWalberlaImpl<CollisionModel, FloatType>>
{
   template< typename FieldPtrOrIterator >
   static void get( Matrix3< real_t > & /* pressureTensor */, const LBWalberlaImpl<CollisionModel, FloatType> & /* latticeModel */, const FieldPtrOrIterator & /* it */ )
   {
       WALBERLA_ABORT("Not implemented");
   }

   template< typename PdfField_T >
   static void get( Matrix3< real_t > & /* pressureTensor */, const LBWalberlaImpl<CollisionModel, FloatType> & /* latticeModel */, const PdfField_T & /* pdf */,
                    const cell_idx_t /* x */, const cell_idx_t /* y */, const cell_idx_t /* z */ )
   {
       WALBERLA_ABORT("Not implemented");
   }
};


template <typename CollisionModel, typename FloatType>
struct ShearRate<LBWalberlaImpl<CollisionModel, FloatType>>
{
   template< typename FieldPtrOrIterator >
   static inline real_t get( const LBWalberlaImpl<CollisionModel, FloatType> & /* latticeModel */, const FieldPtrOrIterator & /* it */,
                             const Vector3< real_t > & /* velocity */, const real_t /* rho */)
   {
       WALBERLA_ABORT("Not implemented");
       return real_t(0.0);
   }

   template< typename PdfField_T >
   static inline real_t get( const LBWalberlaImpl<CollisionModel, FloatType> & latticeModel,
                             const PdfField_T & /* pdf */, const cell_idx_t /* x */, const cell_idx_t /* y */, const cell_idx_t /* z */,
                             const Vector3< real_t > & /* velocity */, const real_t /* rho */ )
   {
       WALBERLA_ABORT("Not implemented");
       return real_t(0.0);
   }

   static inline real_t get( const std::vector< real_t > & /* nonEquilibrium */, const real_t /* relaxationParam */,
                             const real_t /* rho */ = real_t(1) )
   {
       WALBERLA_ABORT("Not implemented");
       return real_t(0.0);
   }
};


} // namespace lbm
} // namespace walberla



#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif

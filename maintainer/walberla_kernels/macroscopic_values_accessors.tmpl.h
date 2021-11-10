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
#include "stencil/{{stencil_name}}.h"

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
namespace {{namespace}} {
namespace accessor {


//======================================================================================================================
//
//  Implementation of macroscopic value backend
//
//======================================================================================================================

namespace EquilibriumDistribution
{
   {{real_t}} get( const stencil::Direction direction,
               const Vector3< {{real_t}} > & u = Vector3< {{real_t}} >( {{real_t}}(0.0) ),
               {{real_t}} rho = {{real_t}}(1.0) )
   {
        {% if not compressible %}
        rho -= {{real_t}}(1.0);
        {% endif %}
        {{equilibrium_from_direction}}
   }
} // EquilibriumDistribution


namespace Equilibrium
{
   template< typename FieldPtrOrIterator, std::enable_if_t<std::is_same<decltype(*(FieldPtrOrIterator())), {{real_t}}>::value, bool> = true >
   void set( FieldPtrOrIterator & it,
             const Vector3< {{real_t}} > & u = Vector3< {{real_t}} >( {{real_t}}(0.0) ),
             {{real_t}} rho = {{real_t}}(1.0) )
   {
        {%if not compressible %}
        rho -= {{real_t}}(1.0);
        {%endif %}

       {% for eqTerm in equilibrium -%}
       it[{{loop.index0 }}] = {{eqTerm}};
       {% endfor -%}
   }

   template< typename PdfField_T, std::enable_if_t<std::is_same<typename PdfField_T::value_type, {{real_t}}>::value, bool> = true >
   void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
             const Vector3< {{real_t}} > & u = Vector3< {{real_t}} >( {{real_t}}(0.0) ),
             {{real_t}} rho = {{real_t}}(1.0) )
   {
      {%if not compressible %}
      rho -= {{real_t}}(1.0);
      {%endif %}

      {{real_t}} & xyz0 = pdf(x,y,z,0);
      {% for eqTerm in equilibrium -%}
      pdf.getF( &xyz0, {{loop.index0 }})= {{eqTerm}};
      {% endfor -%}
   }
} // Equilibrium


namespace Density
{
   template< typename FieldPtrOrIterator, std::enable_if_t<std::is_same<decltype(*(FieldPtrOrIterator())), {{real_t}}>::value, bool> = true >
   inline {{real_t}} get( const FieldPtrOrIterator & it )
   {
        {% for i in range(Q) -%}
            const {{real_t}} f_{{i}} = it[{{i}}];
        {% endfor -%}
        {{density_getters | indent(8)}}
        return rho;
   }

   template< typename PdfField_T, std::enable_if_t<std::is_same<typename PdfField_T::value_type, {{real_t}}>::value, bool> = true >
   inline {{real_t}} get( const PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const {{real_t}} & xyz0 = pdf(x,y,z,0);
        {% for i in range(Q) -%}
            const {{real_t}} f_{{i}} = pdf.getF( &xyz0, {{i}});
        {% endfor -%}
        {{density_getters | indent(8)}}
        return rho;
   }
} // Density


namespace DensityAndVelocity
{
    template< typename FieldPtrOrIterator >
    void set( FieldPtrOrIterator & it,
              const GhostLayerField<{{real_t}}, 3u> & force_field,
              const Vector3< {{real_t}} > & u = Vector3< {{real_t}} >( {{real_t}}(0.0) ),
              const {{real_t}} rho_in = {{real_t}}(1.0) )
    {
        auto x = it.x();
        auto y = it.y();
        auto z = it.z();

        {{density_velocity_setter_macroscopic_values | indent(8)}}

        Equilibrium::set(it, Vector3<{{real_t}}>(u_0, u_1, u_2), rho{%if not compressible %} + {{real_t}}(1) {%endif%});
    }

    template< typename PdfField_T >
    void set( PdfField_T & pdf, const cell_idx_t x, const cell_idx_t y, const cell_idx_t z,
              const GhostLayerField<{{real_t}}, 3u> & force_field,
              const Vector3< {{real_t}} > & u = Vector3< {{real_t}} >( {{real_t}}(0.0) ),
              const {{real_t}} rho_in = {{real_t}}(1.0) )
    {
        {{density_velocity_setter_macroscopic_values | indent(8)}}

        Equilibrium::set(pdf, x, y, z, Vector3<{{real_t}}>(u_0, u_1, u_2), rho {%if not compressible %} + {{real_t}}(1) {%endif%});
    }
} // DensityAndVelocity


namespace DensityAndVelocityRange
{

   template< typename FieldIteratorXYZ >
   void set( FieldIteratorXYZ & begin, const FieldIteratorXYZ & end,
             const GhostLayerField<{{real_t}}, 3u> & force_field,
             const Vector3< {{real_t}} > & u = Vector3< {{real_t}} >( {{real_t}}(0.0) ),
             const {{real_t}} rho_in = {{real_t}}(1.0) )
   {
        for( auto cellIt = begin; cellIt != end; ++cellIt )
        {
            const auto x = cellIt.x();
            const auto y = cellIt.y();
            const auto z = cellIt.z();
            {{density_velocity_setter_macroscopic_values | indent(12)}}

            Equilibrium::set(cellIt, Vector3<{{real_t}}>(u_0, u_1, u_2), rho{%if not compressible %} + {{real_t}}(1) {%endif%});
        }
   }
} // DensityAndVelocityRange



namespace DensityAndMomentumDensity
{
   template< typename FieldPtrOrIterator >
   {{real_t}} get( Vector3< {{real_t}} > & momentumDensity,
               const GhostLayerField<{{real_t}}, 3u> & force_field,
               const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        {% for i in range(Q) -%}
            const {{real_t}} f_{{i}} = it[{{i}}];
        {% endfor -%}

        {{momentum_density_getter | indent(8) }}
        {% for i in range(D) -%}
            momentumDensity[{{i}}] = md_{{i}};
        {% endfor %}
        return rho;
   }

   template< typename PdfField_T >
   {{real_t}} get( Vector3< {{real_t}} > & momentumDensity,
               const GhostLayerField<{{real_t}}, 3u> & force_field,
               const PdfField_T & pdf,
               const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const {{real_t}} & xyz0 = pdf(x,y,z,0);
        {% for i in range(Q) -%}
            const {{real_t}} f_{{i}} = pdf.getF( &xyz0, {{i}});
        {% endfor -%}

        {{momentum_density_getter | indent(8) }}
        {% for i in range(D) -%}
            momentumDensity[{{i}}] = md_{{i}};
        {% endfor %}
       return rho;
   }
} // DensityAndMomentumDensity


namespace MomentumDensity
{
   template< typename FieldPtrOrIterator >
   void get( Vector3< {{real_t}} > & momentumDensity,
             const GhostLayerField<{{real_t}}, 3u> & force_field,
             const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        {% for i in range(Q) -%}
            const {{real_t}} f_{{i}} = it[{{i}}];
        {% endfor -%}

        {{momentum_density_getter | indent(8) }}
        {% for i in range(D) -%}
            momentumDensity[{{i}}] = md_{{i}};
        {% endfor %}
   }

   template< typename PdfField_T >
   void get( Vector3< {{real_t}} > & momentumDensity,
             const GhostLayerField<{{real_t}}, 3u> & force_field,
             const PdfField_T & pdf,
             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const {{real_t}} & xyz0 = pdf(x,y,z,0);
        {% for i in range(Q) -%}
            const {{real_t}} f_{{i}} = pdf.getF( &xyz0, {{i}});
        {% endfor -%}

        {{momentum_density_getter | indent(8) }}
        {% for i in range(D) -%}
            momentumDensity[{{i}}] = md_{{i}};
        {% endfor %}
   }
} // MomentumDensity


namespace PressureTensor
{
   template< typename FieldPtrOrIterator, std::enable_if_t<std::is_same<decltype(*(FieldPtrOrIterator())), {{real_t}}>::value, bool> = true >
   void get( Matrix3< {{real_t}} > & pressureTensor, const FieldPtrOrIterator & it )
   {
        const auto x = it.x();
        const auto y = it.y();
        const auto z = it.z();

        {% for i in range(Q) -%}
            const {{real_t}} f_{{i}} = it[{{i}}];
        {% endfor -%}

        {{second_momentum_getter | indent(8) }}
        {% for i in range(D) -%}
            {% for j in range(D) -%}
                pressureTensor[{{i*D+j}}] = p_{{i*D+j}};
            {% endfor %}
        {% endfor %}
   }

   template< typename PdfField_T, std::enable_if_t<std::is_same<typename PdfField_T::value_type, {{real_t}}>::value, bool> = true >
   void get( Matrix3< {{real_t}} > & pressureTensor, const PdfField_T & pdf,
             const cell_idx_t x, const cell_idx_t y, const cell_idx_t z )
   {
        const {{real_t}} & xyz0 = pdf(x,y,z,0);
        {% for i in range(Q) -%}
            const {{real_t}} f_{{i}} = pdf.getF( &xyz0, {{i}});
        {% endfor -%}

        {{second_momentum_getter | indent(8) }}
        {% for i in range(D) -%}
            {% for j in range(D) -%}
                pressureTensor[{{i*D+j}}] = p_{{i*D+j}};
            {% endfor %}
        {% endfor %}
   }
} // PressureTensor


namespace ShearRate
{
   template< typename FieldPtrOrIterator, std::enable_if_t<std::is_same<decltype(*(FieldPtrOrIterator())), {{real_t}}>::value, bool> = true >
   inline {{real_t}} get( const FieldPtrOrIterator & /* it */,
                          const Vector3< {{real_t}} > & /* velocity */,
                          const {{real_t}} /* rho */)
   {
       WALBERLA_ABORT("Not implemented");
       return {{real_t}}(0.0);
   }

   template< typename PdfField_T, std::enable_if_t<std::is_same<typename PdfField_T::value_type, {{real_t}}>::value, bool> = true >
   inline {{real_t}} get( const PdfField_T & /* pdf */,
                          const cell_idx_t /* x */, const cell_idx_t /* y */, const cell_idx_t /* z */,
                          const Vector3< {{real_t}} > & /* velocity */, const {{real_t}} /* rho */ )
   {
       WALBERLA_ABORT("Not implemented");
       return {{real_t}}(0.0);
   }
} // ShearRate


} // namespace accessor
} // namespace {{namespace}}
} // namespace walberla



#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif

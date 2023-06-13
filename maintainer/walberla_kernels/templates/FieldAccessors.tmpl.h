/*
 * Copyright (C) 2021-2023 The ESPResSo project
 * Copyright (C) 2020 The waLBerla project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/*
 * Lattice field accessors.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/a16141524c58ab88386e2a0f8fdd7c63c5edd704/python/lbmpy_walberla/templates/LatticeModel.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Matrix{{D}}.h>
#include <core/math/Vector{{D}}.h>

#include <field/GhostLayerField.h>
#include <stencil/{{stencil_name}}.h>

#include <array>
#include <cassert>
#include <tuple>
#include <vector>

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

namespace Population
{
    inline std::array<{{dtype}}, {{Q}}u>
    get( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > const * pdf_field,
         Cell const & cell )
    {
        {{dtype}} const & xyz0 = pdf_field->get(cell, uint_t{ 0u });
        std::array<{{dtype}}, {{Q}}u> pop;
        {% for i in range(Q) -%}
            pop[{{i}}u] = pdf_field->getF( &xyz0, uint_t{ {{i}}u });
        {% endfor -%}
        return pop;
    }

    inline void
    set( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > * pdf_field,
         std::array<{{dtype}}, {{Q}}u> const & pop,
         Cell const & cell )
    {
        {{dtype}} & xyz0 = pdf_field->get(cell, uint_t{ 0u });
        {% for i in range(Q) -%}
            pdf_field->getF( &xyz0, uint_t{ {{i}}u }) = pop[{{i}}u];
        {% endfor -%}
    }

    inline void
    broadcast( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > * pdf_field,
               std::array<{{dtype}}, {{Q}}u> const & pop)
     {
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(pdf_field, {
             {{dtype}} & xyz0 = pdf_field->get(x, y, z, uint_t{ 0u });
             {% for i in range(Q) -%}
                 pdf_field->getF( &xyz0, uint_t{ {{i}}u }) = pop[{{i}}u];
             {% endfor -%}
         });
     }

    inline std::vector< {{dtype}} >
    get( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > const * pdf_field,
         CellInterval const & ci )
    {
        std::vector< {{dtype}} > out;
        out.reserve(ci.numCells() * uint_t({{Q}}u));
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    {{dtype}} const & xyz0 = pdf_field->get(x, y, z, uint_t{ 0u });
                    {% for i in range(Q) -%}
                        out.emplace_back(pdf_field->getF( &xyz0, uint_t{ {{i}}u }));
                    {% endfor -%}
                }
            }
        }
        return out;
    }

    inline void
    set( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > * pdf_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci )
    {
        assert(uint_c(values.size()) == ci.numCells() * uint_t({{Q}}u));
        auto values_ptr = values.data();
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    {{dtype}} & xyz0 = pdf_field->get(x, y, z, uint_t{ 0u });
                    {% for i in range(Q) -%}
                        pdf_field->getF( &xyz0, uint_t{ {{i}}u }) = values_ptr[{{i}}u];
                    {% endfor -%}
                    values_ptr += {{Q}}u;
                }
            }
        }
    }
} // namespace Population

namespace Vector
{
    inline Vector{{D}}< {{dtype}} >
    get( GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > const * vec_field,
         Cell const & cell )
    {
        const {{dtype}} & xyz0 = vec_field->get(cell, uint_t{ 0u });
        Vector{{D}}< {{dtype}} > vec;
        {% for i in range(D) -%}
            vec[{{i}}] = vec_field->getF( &xyz0, uint_t{ {{i}}u });
        {% endfor -%}
        return vec;
    }

    inline void
    set( GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > * vec_field,
         Vector{{D}}< {{dtype}} > const & vec,
         Cell const & cell )
    {
        {{dtype}} & xyz0 = vec_field->get(cell, uint_t{ 0u });
        {% for i in range(D) -%}
            vec_field->getF( &xyz0, uint_t{ {{i}}u }) = vec[{{i}}u];
        {% endfor -%}
    }

    inline void
    add( GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > * vec_field,
         Vector{{D}}< {{dtype}} > const & vec,
         Cell const & cell )
    {
        {{dtype}} & xyz0 = vec_field->get(cell, uint_t{ 0u });
        {% for i in range(D) -%}
            vec_field->getF( &xyz0, uint_t{ {{i}}u }) += vec[{{i}}u];
        {% endfor -%}
    }

    inline void
    broadcast( GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > * vec_field,
               Vector{{D}}< {{dtype}} > const & vec)
     {
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
             {{dtype}} & xyz0 = vec_field->get(x, y, z, uint_t{ 0u });
             {% for i in range(D) -%}
                 vec_field->getF( &xyz0, uint_t{ {{i}}u }) = vec[{{i}}u];
             {% endfor -%}
         });
     }

    inline void
    add_to_all( GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > * vec_field,
                Vector{{D}}< {{dtype}} > const & vec)
     {
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(vec_field, {
             {{dtype}} & xyz0 = vec_field->get(x, y, z, uint_t{ 0u });
             {% for i in range(D) -%}
                 vec_field->getF( &xyz0, uint_t{ {{i}}u }) += vec[{{i}}u];
             {% endfor -%}
         });
     }

    inline std::vector< {{dtype}} >
    get( GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > const * vec_field,
         CellInterval const & ci )
    {
        std::vector< {{dtype}} > out;
        out.reserve(ci.numCells() * uint_t({{D}}u));
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    const {{dtype}} & xyz0 = vec_field->get(x, y, z, uint_t{ 0u });
                    {% for i in range(D) -%}
                      out.emplace_back(vec_field->getF( &xyz0, uint_t{ {{i}}u }));
                    {% endfor -%}
                }
            }
        }
        return out;
    }

    inline void
    set( GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > * vec_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci )
    {
        assert(uint_c(values.size()) == ci.numCells() * uint_t({{D}}u));
        auto values_ptr = values.data();
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    {{dtype}} & xyz0 = vec_field->get(x, y, z, uint_t{ 0u });
                    {% for i in range(D) -%}
                        vec_field->getF( &xyz0, uint_t{ {{i}}u }) = values_ptr[{{i}}u];
                    {% endfor -%}
                    values_ptr += {{D}}u;
                }
            }
        }
    }
} // namespace Vector

namespace EquilibriumDistribution
{
    inline {{dtype}}
    get( stencil::Direction const direction,
         Vector{{D}}< {{dtype}} > const & u = Vector{{D}}< {{dtype}} >( {{dtype}}(0.0) ),
         {{dtype}} rho = {{dtype}}(1.0) )
    {
        {% if not compressible %}
        rho -= {{dtype}}(1.0);
        {% endif %}
        {{equilibrium_from_direction}}
    }
} // namespace EquilibriumDistribution

namespace Equilibrium
{
    inline void
    set( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > * pdf_field,
         Vector{{D}}< {{dtype}} > const & u,
         {{dtype}} const rho,
         Cell const & cell )
    {
        {%if not compressible %}
        rho -= {{dtype}}(1.0);
        {%endif %}

        {{dtype}} & xyz0 = pdf_field->get(cell, uint_t{ 0u });
        {% for eqTerm in equilibrium -%}
            pdf_field->getF( &xyz0, uint_t{ {{ loop.index0 }}u }) = {{eqTerm}};
        {% endfor -%}
    }
} // namespace Equilibrium

namespace Density
{
    inline {{dtype}}
    get( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > const * pdf_field,
         Cell const & cell )
    {
        const {{dtype}} & xyz0 = pdf_field->get(cell, uint_t{ 0u });
        {% for i in range(Q) -%}
            const {{dtype}} f_{{i}} = pdf_field->getF( &xyz0, uint_t{ {{i}}u });
        {% endfor -%}
        {{density_getters | indent(8)}}
        return rho;
    }

    inline void
    set( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > * pdf_field,
         {{dtype}} const rho_in,
         Cell const & cell )
    {
        const {{dtype}} & xyz0 = pdf_field->get(cell, uint_t{ 0u });
        {% for i in range(Q) -%}
            const {{dtype}} f_{{i}} = pdf_field->getF( &xyz0, uint_t{ {{i}}u });
        {% endfor -%}

        {{unshifted_momentum_density_getter | indent(8)}}

        // calculate current velocity (before density change)
        const {{dtype}} conversion = {{dtype}}(1) / rho;
        Vector{{D}}< {{dtype}} > velocity;
        {% for i in range(D) -%}
            velocity[{{i}}u] = momdensity_{{i}} * conversion;
        {% endfor %}

        Equilibrium::set(pdf_field, velocity, rho_in {%if not compressible %} + {{dtype}}(1) {%endif%}, cell);
    }

    inline std::vector< {{dtype}} >
    get( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > const * pdf_field,
         CellInterval const & ci )
    {
        std::vector< {{dtype}} > out;
        out.reserve(ci.numCells());
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    const {{dtype}} & xyz0 = pdf_field->get(x, y, z, uint_t{ 0u });
                    {% for i in range(Q) -%}
                        const {{dtype}} f_{{i}} = pdf_field->getF( &xyz0, uint_t{ {{i}}u });
                    {% endfor -%}
                    {{density_getters | indent(12)}}
                    out.emplace_back(rho);
                }
            }
        }
        return out;
    }

    inline void
    set( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > * pdf_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci )
    {
        assert(uint_c(values.size()) == ci.numCells());
        auto values_it = values.begin();
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    const {{dtype}} & xyz0 = pdf_field->get(x, y, z, uint_t{ 0u });
                    {% for i in range(Q) -%}
                        const {{dtype}} f_{{i}} = pdf_field->getF( &xyz0, uint_t{ {{i}}u });
                    {% endfor -%}

                    {{unshifted_momentum_density_getter | indent(12)}}

                    // calculate current velocity (before density change)
                    const {{dtype}} conversion = {{dtype}}(1) / rho;
                    Vector{{D}}< {{dtype}} > velocity;
                    {% for i in range(D) -%}
                        velocity[{{i}}u] = momdensity_{{i}} * conversion;
                    {% endfor %}

                    Equilibrium::set(pdf_field, velocity, *values_it {%if not compressible %} + {{dtype}}(1) {%endif%}, Cell{x, y, z});
                    ++values_it;
                }
            }
        }
    }
} // namespace Density

namespace Velocity
{
    inline void
    set( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > * pdf_field,
         GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > const * force_field,
         Vector{{D}}< {{dtype}} > const & u,
         Cell const & cell )
    {
        const {{dtype}} & xyz0 = pdf_field->get(cell, uint_t{ 0u });
        {% for i in range(Q) -%}
            const {{dtype}} f_{{i}} = pdf_field->getF( &xyz0, uint_t{ {{i}}u });
        {% endfor -%}
        {{density_getters | indent(8)}}

        {% for c in "xyz" -%}
            const auto {{c}} = cell.{{c}}();
        {% endfor -%}
        {{density_velocity_setter_macroscopic_values | substitute_force_getter_cpp | indent(8)}}

        Equilibrium::set(pdf_field, Vector{{D}}<{{dtype}}>({% for i in range(D) %}u_{{i}}{% if not loop.last %}, {% endif %}{% endfor %}), rho {%if not compressible %} + {{dtype}}(1) {%endif%}, cell);
    }
} // namespace Velocity

namespace MomentumDensity
{
    inline Vector{{D}}< {{dtype}} >
    reduce( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > const * pdf_field,
            GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > const * force_field )
    {
        Vector{{D}}< {{dtype}} > momentumDensity({{dtype}} {0});
        WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
            const {{dtype}} & xyz0 = pdf_field->get(x, y, z, uint_t{ 0u });
            {% for i in range(Q) -%}
                const {{dtype}} f_{{i}} = pdf_field->getF( &xyz0, uint_t{ {{i}}u });
            {% endfor -%}

            {{momentum_density_getter | substitute_force_getter_cpp | indent(8) }}

            {% for i in range(D) -%}
                momentumDensity[{{i}}u] += md_{{i}};
            {% endfor %}
        });
        return momentumDensity;
    }
} // namespace MomentumDensity

namespace PressureTensor
{
    inline Matrix{{D}}< {{dtype}} >
    get( GhostLayerField< {{dtype}}, uint_t{ {{Q}}u } > const * pdf_field,
         Cell const & cell )
   {
        const {{dtype}} & xyz0 = pdf_field->get(cell, uint_t{ 0u });
        {% for i in range(Q) -%}
            const {{dtype}} f_{{i}} = pdf_field->getF( &xyz0, uint_t{ {{i}}u });
        {% endfor -%}

        {{second_momentum_getter | indent(8) }}

        Matrix{{D}}< {{dtype}} > pressureTensor;
        {% for i in range(D) -%}
            {% for j in range(D) -%}
                pressureTensor[{{i*D+j}}u] = p_{{i*D+j}};
            {% endfor %}
        {% endfor %}
        return pressureTensor;
   }
} // namespace PressureTensor

} // namespace accessor
} // namespace {{namespace}}
} // namespace walberla

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif

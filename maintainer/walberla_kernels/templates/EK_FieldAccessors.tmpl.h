/*
 * Copyright (C) 2021-2026 The ESPResSo project
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
#include <core/math/Vector{{D}}.h>

#include <field/GhostLayerField.h>

#include <array>
#include <cassert>
#include <iterator>
#include <tuple>
#include <vector>

namespace walberla {
namespace {{namespace}} {
namespace accessor {

namespace Scalar
{
    inline auto
    get( GhostLayerField< {{dtype}}, 1u > const * scalar_field,
         Cell const & cell )
    {
        return scalar_field->get(cell);
    }

    inline void
    set( GhostLayerField< {{dtype}}, 1u > * scalar_field,
        {{dtype}} const & value,
         Cell const & cell )
    {
        scalar_field->get(cell) = value;
    }

    inline void
    add( GhostLayerField< {{dtype}},1u > * scalar_field,
         {{dtype}} const & value,
         Cell const & cell )
    {
        scalar_field->get(cell) += value;
    }

    inline void
    initialize( GhostLayerField< {{dtype}}, 1u > * scalar_field,
                {{dtype}} const & value)
     {
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(scalar_field, {
             scalar_field->get( x, y, z) = value;
         });
     }

    inline void
    add_to_all( GhostLayerField< {{dtype}}, 1u > * scalar_field,
                {{dtype}} const & value)
     {
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(scalar_field, {
             scalar_field->get( x, y, z) += value;
         });
     }

    inline auto
    get( GhostLayerField< {{dtype}}, 1u > const * scalar_field,
         CellInterval const & ci )
    {
        std::vector< {{dtype}} > out;
        out.reserve(ci.numCells());
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    out.emplace_back(scalar_field->get(x, y, z));
                }
            }
        }
        return out;
    }

    inline void
    set( GhostLayerField< {{dtype}}, 1u > * scalar_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci )
    {
        assert(uint_c(values.size()) == ci.numCells());
        auto values_ptr = values.data();
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    scalar_field->get(x, y, z) = values_ptr[0u];
                    std::advance(values_ptr, 1);
                }
            }
        }
    }
} // namespace Scalar

namespace Vector
{
    inline auto
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
    initialize( GhostLayerField< {{dtype}}, uint_t{ {{D}}u } > * vec_field,
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

    inline auto
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
                    std::advance(values_ptr, {{D}});
                }
            }
        }
    }
} // namespace Vector

namespace Flux
{
    inline auto
    get( GhostLayerField< {{dtype}}, uint_t{ {{FluxCount}}u } > const * flux_field,
         Cell const & cell )
    {
        const {{dtype}} & xyz0 = flux_field->get(cell, uint_t{ 0u });
        std::array< {{dtype}}, {{FluxCount}}u > value;
        {% for i in range(FluxCount) -%}
            value[{{i}}] = flux_field->getF( &xyz0, uint_t{ {{i}}u });
        {% endfor -%}
        return value;
    }

    inline void
    initialize( GhostLayerField< {{dtype}}, uint_t{ {{FluxCount}}u } > * flux_field,
                std::array< {{dtype}}, {{FluxCount}} > const & values)
     {
         WALBERLA_FOR_ALL_CELLS_INCLUDING_GHOST_LAYER_XYZ(flux_field, {
             {{dtype}} & xyz0 = flux_field->get(x, y, z, uint_t{ 0u });
             {% for i in range(FluxCount) -%}
                 flux_field->getF( &xyz0, uint_t{ {{i}}u }) = values[{{i}}u];
             {% endfor -%}
         });
     }

    inline auto
    get( GhostLayerField< {{dtype}}, uint_t{ {{FluxCount}}u } > const * flux_field,
         CellInterval const & ci )
    {
        std::vector< {{dtype}} > out;
        out.reserve(ci.numCells() * uint_t({{FluxCount}}u));
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    const {{dtype}} & xyz0 = flux_field->get(x, y, z, uint_t{ 0u });
                    {% for i in range(D) -%}
                      out.emplace_back(flux_field->getF( &xyz0, uint_t{ {{i}}u }));
                    {% endfor -%}
                }
            }
        }
        return out;
    }

    inline auto
    get_vector( GhostLayerField< {{dtype}}, uint_t{ {{FluxCount}}u } > const * flux_field,
         Cell const & cell )
    {
        Vector{{D}}< {{dtype}} > result = Vector{{D}}< {{dtype}} >(0,0,0);
        const {{dtype}} & xyz0 = flux_field->get(cell, uint_t{ 0u });
        std::array< {{dtype}}, {{2*FluxCount+1}}u > local_value;
        // get fluxes in all directions
        {% for i in range(FluxCount*2+1) -%}
            {% if i == 0 -%}
                local_value[{{i}}] = {{dtype}}(0.0);
            {% elif Stencils[i] in StaggeredStencils -%}
                local_value[{{i}}] = {{dtype}}(0.5) * flux_field->getF( &xyz0, uint_t{ {{StaggeredStencils[Stencils[i]]}}u });
            {% else -%}
                local_value[{{i}}] = {{dtype}}(-0.5) * flux_field->getNeighbor(cell.x(), cell.y(), cell.z(), uint_t{ {{InverseStencils[Stencils[i]]}}u }, stencil::Direction(uint_t{ {{i}}u }));
            {% endif -%}
        {% endfor %}

        {% for i in range(FluxCount*2+1) -%}
            {% for j in range(3) -%}
                {% if Stencils[i][j] != 0-%}
                    result[{{j}}] += local_value[{{i}}] * {{Stencils[i][j]}};
                {% endif -%}
            {% endfor %}  
        {% endfor %}
        return result;
    }

    inline auto
    get_vector( GhostLayerField< {{dtype}}, uint_t{ {{FluxCount}}u } > const * flux_field,
         CellInterval const & ci )
    {
        std::vector< {{dtype}} > out;
        out.reserve(ci.numCells() * uint_t({{D}}u));
        for (auto x = ci.xMin(); x <= ci.xMax(); ++x) {
            for (auto y = ci.yMin(); y <= ci.yMax(); ++y) {
                for (auto z = ci.zMin(); z <= ci.zMax(); ++z) {
                    Vector{{D}}< {{dtype}} > result = Vector{{D}}< {{dtype}} >(0,0,0);
                    const {{dtype}} & xyz0 = flux_field->get(x, y, z, uint_t{ 0u });
                    std::array< {{dtype}}, {{2*FluxCount+1}}u > local_value;
                    // get fluxes in all directions
                    {% for i in range(FluxCount*2+1) -%}
                        {% if i == 0 -%}
                            local_value[{{i}}] = {{dtype}}(0.0);
                        {% elif Stencils[i] in StaggeredStencils -%}
                            local_value[{{i}}] = {{dtype}}(0.5) * flux_field->getF( &xyz0, uint_t{ {{StaggeredStencils[Stencils[i]]}}u });
                        {% else -%}
                            local_value[{{i}}] = {{dtype}}(-0.5) * flux_field->getNeighbor(x, y, z, uint_t{ {{InverseStencils[Stencils[i]]}}u }, stencil::Direction(uint_t{ {{i}}u }));
                        {% endif -%}
                    {% endfor %}

                    {% for i in range(FluxCount*2+1) -%}
                        {% for j in range(3) -%}
                            {% if Stencils[i][j] != 0-%}
                                result[{{j}}] += local_value[{{i}}] * {{Stencils[i][j]}};
                            {% endif -%}
                        {% endfor %}  
                    {% endfor %}
                    {% for i in range(D) -%}
                        out.emplace_back(result[{{i}}u]);
                    {% endfor -%}
                }
            }
        }
        return out;
    }
} // namespace Flux

} // namespace accessor
} // namespace {{namespace}}
} // namespace walberla

/*
 * Copyright (C) 2023-2025 The ESPResSo project
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

#include <gpu/GPUField.h>

#include <array>
#include <tuple>
#include <vector>

namespace walberla {
namespace {{namespace}} {
namespace accessor {

namespace Scalar {
    /** @brief Get value from a single cell. */
    {{dtype}}
    get( gpu::GPUField< {{dtype}} > const * field,
         Cell const & cell );
    /** @brief Set value on a single cell. */
    void set( gpu::GPUField< {{dtype}} > * field,
              {{dtype}} value,
              Cell const & cell );
    /** @brief Add value to a single cell. */
    void add( gpu::GPUField< {{dtype}} > * field,
              {{dtype}} value,
              Cell const & cell );
    /** @brief Initialize all cells with the same value. */
    void initialize( gpu::GPUField< {{dtype}} > * field,
                    {{dtype}} value);
    /** @brief Add value to all cells. */
    void add_to_all( gpu::GPUField< {{dtype}} > * field,
                     {{dtype}} value);
    /** @brief Get values from a cell interval. */
    std::vector< {{dtype}} >
    get( gpu::GPUField< {{dtype}} > const * ield,
         CellInterval const & ci);
    /** @brief Set values on a cell interval. */
    void
    set( gpu::GPUField< {{dtype}} > * field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci );

} // namespace Scalar

namespace Vector {
    /** @brief Get value from a single cell. */
    Vector{{D}}< {{dtype}} >
    get( gpu::GPUField< {{dtype}} > const * field,
         Cell const & cell );
    /** @brief Set value on a single cell. */
    void set( gpu::GPUField< {{dtype}} > * field,
              Vector{{D}}< {{dtype}} > const & vec,
              Cell const & cell );
    /** @brief Add value to a single cell. */
    void add( gpu::GPUField< {{dtype}} > * field,
              Vector{{D}}< {{dtype}} > const & vec,
              Cell const & cell );
    /** @brief Initialize all cells with the same value. */
    void initialize( gpu::GPUField< {{dtype}} > * field,
                    Vector{{D}}< {{dtype}} > const & vec);
    /** @brief Add value to all cells. */
    void add_to_all( gpu::GPUField< {{dtype}} > * field,
                     Vector{{D}}< {{dtype}} > const & vec);
    /** @brief Get values from a cell interval. */
    std::vector< {{dtype}} >
    get( gpu::GPUField< {{dtype}} > const * vec_field,
         CellInterval const & ci);
    /** @brief Set values on a cell interval. */
    void
    set( gpu::GPUField< {{dtype}} > * vec_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci );

} // namespace Vector

namespace Flux {
    /** @brief Get value from a single cell. */
    std::array< {{dtype}}, {{FluxCount}} >
    get( gpu::GPUField< {{dtype}} > const * flux_field,
         Cell const & cell );

    /** @brief Initialize all cells with the same value. */
    void initialize( gpu::GPUField< {{dtype}} > * flux_field,
                     std::array< {{dtype}}, {{FluxCount}} > const & flux);
    /** @brief Get values from a cell interval. */
    std::vector< {{dtype}} >
    get( gpu::GPUField< {{dtype}} > const * flux_field,
         CellInterval const & ci);

    /** @brief Get flux vector from a single cell. */
    Vector{{D}}< {{dtype}} >
    get_vector( gpu::GPUField< {{dtype}} > const * flux_field,
         Cell const & cell );

     /** @brief Get flux vector from a cell interval. */
    std::vector< {{dtype}} >
    get_vector( gpu::GPUField< {{dtype}} > const * flux_field,
         CellInterval const & ci);

} // namespace Flux

} // namespace accessor
} // namespace {{namespace}}
} // namespace walberla

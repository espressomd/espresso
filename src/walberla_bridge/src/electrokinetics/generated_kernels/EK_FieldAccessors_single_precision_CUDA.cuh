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

// kernel generated with pystencils v1.3.7, lbmpy v1.3.7, sympy v1.12.1,
// lbmpy_walberla/pystencils_walberla from waLBerla commit
// c69cb11d6a95d32b2280544d3d9abde1fe5fdbb5

/*
 * Lattice field accessors.
 * Adapted from the waLBerla source file
 * https://i10git.cs.fau.de/walberla/walberla/-/blob/a16141524c58ab88386e2a0f8fdd7c63c5edd704/python/lbmpy_walberla/templates/LatticeModel.tmpl.h
 */

#pragma once

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Matrix3.h>
#include <core/math/Vector3.h>

#include <gpu/GPUField.h>

#include <array>
#include <tuple>
#include <vector>

namespace walberla {
namespace ek {
namespace accessor {

namespace Scalar {
/** @brief Get value from a single cell. */
float get(gpu::GPUField<float> const *field, Cell const &cell);
/** @brief Set value on a single cell. */
void set(gpu::GPUField<float> *field, float value, Cell const &cell);
/** @brief Add value to a single cell. */
void add(gpu::GPUField<float> *field, float value, Cell const &cell);
/** @brief Initialize all cells with the same value. */
void initialize(gpu::GPUField<float> *field, float value);
/** @brief Add value to all cells. */
void add_to_all(gpu::GPUField<float> *field, float value);
/** @brief Get values from a cell interval. */
std::vector<float> get(gpu::GPUField<float> const *ield,
                       CellInterval const &ci);
/** @brief Set values on a cell interval. */
void set(gpu::GPUField<float> *field, std::vector<float> const &values,
         CellInterval const &ci);

} // namespace Scalar

namespace Vector {
/** @brief Get value from a single cell. */
Vector3<float> get(gpu::GPUField<float> const *field, Cell const &cell);
/** @brief Set value on a single cell. */
void set(gpu::GPUField<float> *field, Vector3<float> const &vec,
         Cell const &cell);
/** @brief Add value to a single cell. */
void add(gpu::GPUField<float> *field, Vector3<float> const &vec,
         Cell const &cell);
/** @brief Initialize all cells with the same value. */
void initialize(gpu::GPUField<float> *field, Vector3<float> const &vec);
/** @brief Add value to all cells. */
void add_to_all(gpu::GPUField<float> *field, Vector3<float> const &vec);
/** @brief Get values from a cell interval. */
std::vector<float> get(gpu::GPUField<float> const *vec_field,
                       CellInterval const &ci);
/** @brief Set values on a cell interval. */
void set(gpu::GPUField<float> *vec_field, std::vector<float> const &values,
         CellInterval const &ci);

} // namespace Vector

namespace Flux {
/** @brief Get value from a single cell. */
std::array<float, 13> get(gpu::GPUField<float> const *flux_field,
                          Cell const &cell);

/** @brief Initialize all cells with the same value. */
void initialize(gpu::GPUField<float> *flux_field,
                std::array<float, 13> const &flux);
/** @brief Get values from a cell interval. */
std::vector<float> get(gpu::GPUField<float> const *flux_field,
                       CellInterval const &ci);

/** @brief Get flux vector from a single cell. */
Vector3<float> get_vector(gpu::GPUField<float> const *flux_field,
                          Cell const &cell);

/** @brief Get flux vector from a cell interval. */
std::vector<float> get_vector(gpu::GPUField<float> const *flux_field,
                              CellInterval const &ci);

} // namespace Flux

} // namespace accessor
} // namespace ek
} // namespace walberla

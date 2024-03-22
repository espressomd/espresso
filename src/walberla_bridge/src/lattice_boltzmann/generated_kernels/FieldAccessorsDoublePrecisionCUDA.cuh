/*
 * Copyright (C) 2023 The ESPResSo project
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

// kernel generated with pystencils v1.2, lbmpy v1.2,
// lbmpy_walberla/pystencils_walberla from waLBerla commit
// 065ce5f311850371a97ac4766f47dbb5ca8424ba

/**
 * @file
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

#include <cuda/GPUField.h>

#include <array>
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
namespace lbm {
namespace accessor {

namespace Population {
std::array<double, 19u> get(cuda::GPUField<double> const *pdf_field,
                            Cell const &cell);
void set(cuda::GPUField<double> *pdf_field, std::array<double, 19u> const &pop,
         Cell const &cell);
void broadcast(cuda::GPUField<double> *pdf_field,
               std::array<double, 19u> const &pop);
std::vector<double> get(cuda::GPUField<double> const *pdf_field,
                        CellInterval const &ci);
void set(cuda::GPUField<double> *pdf_field, std::vector<double> const &values,
         CellInterval const &ci);
} // namespace Population

namespace Vector {
Vector3<double> get(cuda::GPUField<double> const *field, Cell const &cell);
void set(cuda::GPUField<double> *field, Vector3<double> const &vec,
         Cell const &cell);
void add(cuda::GPUField<double> *field, Vector3<double> const &vec,
         Cell const &cell);
void add_at(cuda::GPUField<double> *vec_field, std::vector<double> const &vecs,
            std::vector<cell_idx_t> const &cells);
std::vector<double> get_at(cuda::GPUField<double> *vec_field,
                           std::vector<cell_idx_t> const &cells);
void broadcast(cuda::GPUField<double> *field, Vector3<double> const &vec);
void add_to_all(cuda::GPUField<double> *field, Vector3<double> const &vec);
std::vector<double> get(cuda::GPUField<double> const *vec_field,
                        CellInterval const &ci);
void set(cuda::GPUField<double> *vec_field, std::vector<double> const &values,
         CellInterval const &ci);
} // namespace Vector

namespace Density {
double get(cuda::GPUField<double> const *pdf_field, Cell const &cell);
void set(cuda::GPUField<double> *pdf_field, double const rho, Cell const &cell);
std::vector<double> get(cuda::GPUField<double> const *pdf_field,
                        CellInterval const &ci);
void set(cuda::GPUField<double> *pdf_field, std::vector<double> const &values,
         CellInterval const &ci);
} // namespace Density

namespace Velocity {
void set(cuda::GPUField<double> *pdf_field, cuda::GPUField<double> *force_field,
         Vector3<double> const &u, Cell const &cell);
} // namespace Velocity

namespace DensityAndVelocity {
std::tuple<double, Vector3<double>>
get(cuda::GPUField<double> const *pdf_field,
    cuda::GPUField<double> const *force_field, Cell const &cell);
void set(cuda::GPUField<double> *pdf_field, cuda::GPUField<double> *force_field,
         Vector3<double> const &u, double const rho, Cell const &cell);
} // namespace DensityAndVelocity

namespace DensityAndMomentumDensity {
std::tuple<double, Vector3<double>>
get(cuda::GPUField<double> const *pdf_field,
    cuda::GPUField<double> const *force_field, Cell const &cell);
} // namespace DensityAndMomentumDensity

namespace MomentumDensity {
Vector3<double> reduce(cuda::GPUField<double> const *pdf_field,
                       cuda::GPUField<double> const *force_field);
} // namespace MomentumDensity

namespace PressureTensor {
Matrix3<double> get(cuda::GPUField<double> const *pdf_field, Cell const &cell);
} // namespace PressureTensor

} // namespace accessor
} // namespace lbm
} // namespace walberla

#ifdef WALBERLA_CXX_COMPILER_IS_GNU
#pragma GCC diagnostic pop
#endif

#ifdef WALBERLA_CXX_COMPILER_IS_CLANG
#pragma clang diagnostic pop
#endif
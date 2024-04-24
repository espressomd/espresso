/*
 * Copyright (C) 2023-2024 The ESPResSo project
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
// 0c8b4b926c6979288fd8a6846d02ec0870e1fe41

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

#include <gpu/GPUField.h>

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
std::array<float, 19u> get(gpu::GPUField<float> const *pdf_field,
                           Cell const &cell);
void set(gpu::GPUField<float> *pdf_field, std::array<float, 19u> const &pop,
         Cell const &cell);
void broadcast(gpu::GPUField<float> *pdf_field,
               std::array<float, 19u> const &pop);
std::vector<float> get(gpu::GPUField<float> const *pdf_field,
                       CellInterval const &ci);
void set(gpu::GPUField<float> *pdf_field, std::vector<float> const &values,
         CellInterval const &ci);
} // namespace Population

namespace Vector {
Vector3<float> get(gpu::GPUField<float> const *field, Cell const &cell);
void set(gpu::GPUField<float> *field, Vector3<float> const &vec,
         Cell const &cell);
void add(gpu::GPUField<float> *field, Vector3<float> const &vec,
         Cell const &cell);
void broadcast(gpu::GPUField<float> *field, Vector3<float> const &vec);
void add_to_all(gpu::GPUField<float> *field, Vector3<float> const &vec);
std::vector<float> get(gpu::GPUField<float> const *vec_field,
                       CellInterval const &ci);
void set(gpu::GPUField<float> *vec_field, std::vector<float> const &values,
         CellInterval const &ci);

} // namespace Vector

namespace Coupling {
std::vector<float> get_interpolated(gpu::GPUField<float> const *vec_field,
                                    std::vector<float> const &pos, uint gl);
void set_interpolated(gpu::GPUField<float> const *vec_field,
                      std::vector<float> const &pos,
                      std::vector<float> const &forces, uint gl);
} // namespace Coupling

namespace Density {
float get(gpu::GPUField<float> const *pdf_field, Cell const &cell);
void set(gpu::GPUField<float> *pdf_field, float const rho, Cell const &cell);
std::vector<float> get(gpu::GPUField<float> const *pdf_field,
                       CellInterval const &ci);
void set(gpu::GPUField<float> *pdf_field, std::vector<float> const &values,
         CellInterval const &ci);
} // namespace Density

namespace Velocity {
void set(gpu::GPUField<float> *pdf_field, gpu::GPUField<float> *force_field,
         Vector3<float> const &u, Cell const &cell);
} // namespace Velocity

namespace DensityAndVelocity {
std::tuple<float, Vector3<float>> get(gpu::GPUField<float> const *pdf_field,
                                      gpu::GPUField<float> const *force_field,
                                      Cell const &cell);
void set(gpu::GPUField<float> *pdf_field, gpu::GPUField<float> *force_field,
         Vector3<float> const &u, float const rho, Cell const &cell);
} // namespace DensityAndVelocity

namespace DensityAndMomentumDensity {
std::tuple<float, Vector3<float>> get(gpu::GPUField<float> const *pdf_field,
                                      gpu::GPUField<float> const *force_field,
                                      Cell const &cell);
} // namespace DensityAndMomentumDensity

namespace MomentumDensity {
Vector3<float> reduce(gpu::GPUField<float> const *pdf_field,
                      gpu::GPUField<float> const *force_field);
} // namespace MomentumDensity

namespace PressureTensor {
Matrix3<float> get(gpu::GPUField<float> const *pdf_field, Cell const &cell);
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

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
#include <core/math/Matrix{{D}}.h>
#include <core/math/Vector{{D}}.h>

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
namespace {{namespace}} {
namespace accessor {

namespace Population {
    std::array<{{dtype}}, {{Q}}u>
    get( cuda::GPUField< {{dtype}} > const * pdf_field,
         Cell const & cell );
    void
    set( cuda::GPUField< {{dtype}} > * pdf_field,
         std::array< {{dtype}}, {{Q}}u > const & pop,
         Cell const & cell );
    void broadcast(
         cuda::GPUField< {{dtype}} > * pdf_field,
         std::array< {{dtype}}, {{Q}}u > const & pop );
    std::vector< {{dtype}} >
    get( cuda::GPUField< {{dtype}} > const * pdf_field,
         CellInterval const & ci );
    void
    set( cuda::GPUField< {{dtype}} > * pdf_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci );
} // namespace Population

namespace Vector {
    Vector{{D}}< {{dtype}} >
    get( cuda::GPUField< {{dtype}} > const * field,
         Cell const & cell );
    void set( cuda::GPUField< {{dtype}} > * field,
              Vector{{D}}< {{dtype}} > const & vec,
              Cell const & cell );
    void add( cuda::GPUField< {{dtype}} > * field,
              Vector{{D}}< {{dtype}} > const & vec,
              Cell const & cell );
    void broadcast( cuda::GPUField< {{dtype}} > * field,
                    Vector{{D}}< {{dtype}} > const & vec);
    void add_to_all( cuda::GPUField< {{dtype}} > * field,
                     Vector{{D}}< {{dtype}} > const & vec);
    std::vector< {{dtype}} >
    get( cuda::GPUField< {{dtype}} > const * vec_field,
         CellInterval const & ci);
    void
    set( cuda::GPUField< {{dtype}} > * vec_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci );
} // namespace Vector

namespace Density {
    {{dtype}}
    get( cuda::GPUField< {{dtype}} > const * pdf_field,
         Cell const & cell );
    void
    set( cuda::GPUField< {{dtype}} > * pdf_field,
         {{dtype}} const rho,
         Cell const & cell );
    std::vector< {{dtype}} >
    get( cuda::GPUField< {{dtype}} > const * pdf_field,
         CellInterval const & ci );
    void
    set( cuda::GPUField< {{dtype}} > * pdf_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci );
} // namespace Density

namespace Velocity {
    void
    set( cuda::GPUField< {{dtype}} > * pdf_field,
         cuda::GPUField< {{dtype}} > * force_field,
         Vector{{D}}< {{dtype}} > const & u,
         Cell const & cell );
} // namespace Velocity

namespace DensityAndVelocity {
    std::tuple< {{dtype}} , Vector{{D}}< {{dtype}} > >
    get( cuda::GPUField< {{dtype}} > const * pdf_field,
         cuda::GPUField< {{dtype}} > const * force_field,
         Cell const & cell );
    void
    set( cuda::GPUField< {{dtype}} > * pdf_field,
         cuda::GPUField< {{dtype}} > * force_field,
         Vector{{D}}< {{dtype}} > const & u,
         {{dtype}} const rho,
         Cell const & cell );
} // namespace DensityAndVelocity

namespace DensityAndMomentumDensity {
    std::tuple< {{dtype}} , Vector{{D}}< {{dtype}} > >
    get( cuda::GPUField< {{dtype}} > const * pdf_field,
         cuda::GPUField< {{dtype}} > const * force_field,
         Cell const & cell );
} // namespace DensityAndMomentumDensity

namespace MomentumDensity {
    Vector{{D}}< {{dtype}} >
    reduce( cuda::GPUField< {{dtype}} > const * pdf_field,
            cuda::GPUField< {{dtype}} > const * force_field );
} // namespace MomentumDensity

namespace PressureTensor {
    Matrix{{D}}< {{dtype}} >
    get( cuda::GPUField< {{dtype}} > const * pdf_field,
         Cell const & cell );
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

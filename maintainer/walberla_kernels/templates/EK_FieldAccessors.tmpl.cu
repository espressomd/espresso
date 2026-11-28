/*
 * Copyright (C) 2023-2026 The ESPResSo project
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

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Vector{{D}}.h>

#include <field/iterators/IteratorMacros.h>

#include <gpu/FieldAccessor.h>
#include <gpu/FieldIndexing.h>
#include <gpu/GPUField.h>
#include <gpu/Kernel.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

#include <array>
#include <vector>

#if defined(__NVCC__)
#define RESTRICT __restrict__
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
#define RESTRICT __restrict__
#else
// clang compiling CUDA code in host mode
#define RESTRICT __restrict__
#endif
#endif
#elif defined(__GNUC__) or defined(__GNUG__)
#define RESTRICT __restrict__
#elif defined(_MSC_VER)
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

/** @brief Get linear index of flattened data with original layout @c fzyx. */
static __forceinline__ __device__ uint getLinearIndex( uint3 blockIdx, uint3 threadIdx, uint3 gridDim, uint3 blockDim, uint fOffset ) {
  auto const x = threadIdx.x;
  auto const y = blockIdx.x;
  auto const z = blockIdx.y;
  auto const f = blockIdx.z;
  auto const ySize = gridDim.x;
  auto const zSize = gridDim.y;
  auto const fSize = fOffset;
  return f                         +
         z * fSize                 +
         y * fSize * zSize         +
         x * fSize * zSize * ySize ;
}

namespace walberla {
namespace {{namespace}} {
namespace accessor {

namespace Scalar
{
// LCOV_EXCL_START
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > scalar_field,
        {{dtype}} * out )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 1u);
        scalar_field.set( blockIdx, threadIdx );
        out += offset;
        if (scalar_field.isValidPosition()) {
            out[0u] = scalar_field.get(0u);
        }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > scalar_field,
        {{dtype}} const * RESTRICT in )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 1u);
        scalar_field.set( blockIdx, threadIdx );
        in += offset;
        if (scalar_field.isValidPosition()) {
            scalar_field.get(0u) = in[0u];
        }
    }

    __global__ void kernel_broadcast(
        gpu::FieldAccessor< {{dtype}} > scalar_field,
        {{dtype}} const in )
    {
        scalar_field.set( blockIdx, threadIdx );
        if (scalar_field.isValidPosition()) {
            scalar_field.get(0u) = in;
        }
    }

    __global__ void kernel_add(
        gpu::FieldAccessor< {{dtype}} > scalar_field,
        {{dtype}} const * RESTRICT in )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 1u);
        scalar_field.set( blockIdx, threadIdx );
        in += offset;
        if (scalar_field.isValidPosition()) {
            scalar_field.get(0u) += in[0u];
        }
    }

    __global__ void kernel_broadcast_add(
        gpu::FieldAccessor< {{dtype}} > scalar_field,
        {{dtype}} const in )
    {
        scalar_field.set( blockIdx, threadIdx );
        if (scalar_field.isValidPosition()) {
            scalar_field.get(0u) += in;
        }
    }
// LCOV_EXCL_STOP

    {{dtype}} get(
        gpu::GPUField< {{dtype}} > const * scalar_field,
        Cell const & cell)
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(1u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *scalar_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        {{dtype}} result{};
        thrust::copy(dev_data.begin(), dev_data.end(), &result);
        return result;
    }

    void set(
        gpu::GPUField< {{dtype}} > * scalar_field,
        {{dtype}} const value,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        auto kernel = gpu::make_kernel( kernel_broadcast );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *scalar_field, ci ) );
        kernel.addParam( value );
        kernel();
    }

    void add(
        gpu::GPUField< {{dtype}} > * scalar_field,
        {{dtype}} const value,
        Cell const &cell )
    {
        CellInterval ci ( cell, cell );
        auto kernel = gpu::make_kernel( kernel_add );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *scalar_field, ci ) );
        kernel.addParam( value );
        kernel();
    }

    void initialize(
        gpu::GPUField< {{dtype}} > * scalar_field,
        {{dtype}} const value )
    {
        CellInterval ci = scalar_field->xyzSizeWithGhostLayer();
        auto kernel = gpu::make_kernel( kernel_broadcast );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *scalar_field, ci ) );
        kernel.addParam( value );
        kernel();
   }

    void add_to_all(
        gpu::GPUField< {{dtype}} > * scalar_field,
        Vector{{D}}< {{dtype}} > const value )
    {
        CellInterval ci = scalar_field->xyzSizeWithGhostLayer();
        auto kernel = gpu::make_kernel( kernel_broadcast_add );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *scalar_field, ci ) );
        kernel.addParam( value );
        kernel();
    }

    std::vector< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * scalar_field,
        CellInterval const & ci)
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells() * 1u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *scalar_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(ci.numCells());
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
    }

    void set(
        gpu::GPUField< {{dtype}} > * scalar_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *scalar_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Scalar

namespace Vector
{
// LCOV_EXCL_START
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} * u_out )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
        vec.set( blockIdx, threadIdx );
        u_out += offset;
        if (vec.isValidPosition()) {
            {% for i in range(D) -%}
                u_out[{{i}}u] = vec.get({{i}}u);
            {% endfor %}
        }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const * RESTRICT u_in )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
        vec.set( blockIdx, threadIdx );
        u_in += offset;
        if (vec.isValidPosition()) {
            {% for i in range(D) -%}
                vec.get({{i}}u) = u_in[{{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_broadcast(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const * RESTRICT u_in )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            {% for i in range(D) -%}
                vec.get({{i}}u) = u_in[{{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_add(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const * RESTRICT u_in )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
        vec.set( blockIdx, threadIdx );
        u_in += offset;
        if (vec.isValidPosition()) {
            {% for i in range(D) -%}
                vec.get({{i}}u) += u_in[{{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_broadcast_add(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const * RESTRICT u_in )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            {% for i in range(D) -%}
                vec.get({{i}}u) += u_in[{{i}}u];
            {% endfor %}
        }
    }
// LCOV_EXCL_STOP

    Vector{{D}}< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * vec_field,
        Cell const & cell)
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        Vector{{D}}< {{dtype}} > vec;
        thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
        return vec;
    }

    void set(
        gpu::GPUField< {{dtype}} > * vec_field,
        Vector{{D}}< {{dtype}} > const & vec,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(vec.data(), vec.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void add(
        gpu::GPUField< {{dtype}} > * vec_field,
        Vector{{D}}< {{dtype}} > const & vec,
        Cell const &cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(vec.data(), vec.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_add );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void initialize(
        gpu::GPUField< {{dtype}} > * vec_field,
        Vector{{D}}< {{dtype}} > const & vec )
    {
        CellInterval ci = vec_field->xyzSizeWithGhostLayer();
        thrust::device_vector< {{dtype}} > dev_data(vec.data(), vec.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_broadcast );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
   }

    void add_to_all(
        gpu::GPUField< {{dtype}} > * vec_field,
        Vector{{D}}< {{dtype}} > const & vec )
    {
        CellInterval ci = vec_field->xyzSizeWithGhostLayer();
        thrust::device_vector< {{dtype}} > dev_data(vec.data(), vec.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_broadcast_add );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    std::vector< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * vec_field,
        CellInterval const & ci)
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells() * {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(ci.numCells() * {{D}}u);
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
    }

    void set(
        gpu::GPUField< {{dtype}} > * vec_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Vector

namespace Flux
{
// LCOV_EXCL_START
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > flux_field,
        {{dtype}} * j_out )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{FluxCount}}u);
        flux_field.set( blockIdx, threadIdx );
        j_out += offset;
        if (flux_field.isValidPosition()) {
            {% for i in range(FluxCount) -%}
                j_out[{{i}}u] = flux_field.get({{i}}u);
            {% endfor %}
        }
    }

    __global__ void kernel_broadcast(
        gpu::FieldAccessor< {{dtype}} > flux_field,
        {{dtype}} const * RESTRICT j_in )
    {
        flux_field.set( blockIdx, threadIdx );
        if (flux_field.isValidPosition()) {
            {% for i in range(FluxCount) -%}
                flux_field.get({{i}}u) = j_in[{{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_get_vector(
        gpu::FieldAccessor< {{dtype}} > flux_field,
        {{dtype}} * j_out )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
        flux_field.set( blockIdx, threadIdx );
        j_out += offset;
        if (flux_field.isValidPosition()) {
            {% for i in range(D) -%}
                j_out[{{i}}u] = {{dtype}}(0.0);
            {% endfor %}
            int cx = 0;
            int cy = 0;
            int cz = 0;
            {{dtype}} add_flux;

            {% for i in range(1,2*FluxCount+1) -%}
                {% if Stencils[i] in StaggeredStencils -%}
                    add_flux = {{dtype}}(0.5) * flux_field.get({{StaggeredStencils[Stencils[i]]}}u);
                {% else -%}
                    cx = {{Stencils[i][0]}};
                    cy = {{Stencils[i][1]}};
                    cz = {{Stencils[i][2]}};
                    add_flux = {{dtype}}(-0.5) * flux_field.getNeighbor(cx, cy, cz, {{InverseStencils[Stencils[i]]}}u);
                {% endif -%}
                {% for j in range(3) -%}
                    {% if Stencils[i][j] != 0-%}
                        j_out[{{j}}u] += add_flux * {{Stencils[i][j]}};
                    {% endif -%}
                {% endfor -%}
            {% endfor %}
        }
    }
// LCOV_EXCL_STOP

    std::array< {{dtype}}, {{FluxCount}} > get(
        gpu::GPUField< {{dtype}} > const * flux_field,
        Cell const & cell)
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{FluxCount}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *flux_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::array< {{dtype}}, {{FluxCount}} > vec;
        thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
        return vec;
    }

    void initialize(
        gpu::GPUField< {{dtype}} > * flux_field,
        std::array< {{dtype}}, {{FluxCount}}> const & flux )
    {
        CellInterval ci = flux_field->xyzSizeWithGhostLayer();
        thrust::device_vector< {{dtype}} > dev_data(flux.data(), flux.data() + {{FluxCount}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_broadcast );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *flux_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
   }

    std::vector< {{dtype}}> get(
        gpu::GPUField< {{dtype}} > const * flux_field,
        CellInterval const & ci)
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells() * {{FluxCount}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *flux_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(ci.numCells() * {{FluxCount}}u);
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
    }

    Vector{{D}}< {{dtype}} > get_vector(
        gpu::GPUField< {{dtype}} > const * flux_field,
        Cell const & cell)
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get_vector );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *flux_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        Vector{{D}}< {{dtype}} > vec;
        thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
        return vec;
    }

    std::vector< {{dtype}}> get_vector(
        gpu::GPUField< {{dtype}} > const * flux_field,
        CellInterval const & ci)
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells() * {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get_vector );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *flux_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(ci.numCells() * {{D}}u);
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
    }
} // namespace Flux

} // namespace accessor
} // namespace {{namespace}}
} // namespace walberla

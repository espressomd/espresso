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

#include <core/DataTypes.h>
#include <core/cell/Cell.h>
#include <core/cell/CellInterval.h>
#include <core/math/Matrix{{D}}.h>
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
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // unused variable
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#else
// clang compiling CUDA code in host mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#endif
#endif
#elif defined(__GNUC__) or defined(__GNUG__)
#define RESTRICT __restrict__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
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

namespace Population
{
// LCOV_EXCL_START
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT pop )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
        pdf.set( blockIdx, threadIdx );
        pop += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                pop[{{i}}u] = pdf.get({{i}}u);
            {% endfor -%}
        }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} const * RESTRICT pop )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
        pdf.set( blockIdx, threadIdx );
        pop += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                pdf.get({{i}}u) = pop[{{i}}u];
            {% endfor -%}
        }
    }

    __global__ void kernel_broadcast(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} const * RESTRICT pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                pdf.get({{i}}u) = pop[{{i}}u];
            {% endfor -%}
        }
    }

    __global__ void kernel_set_vel(
        gpu::FieldAccessor< {{dtype}} > pdf,
        gpu::FieldAccessor< {{dtype}} > velocity,
        gpu::FieldAccessor< {{dtype}} > force,
        {{dtype}} const * RESTRICT pop )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
        pdf.set( blockIdx, threadIdx );
        velocity.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        pop += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                const {{dtype}} f_{{i}} = pdf.get({{i}}u) = pop[{{i}}u];
            {% endfor -%}
            {{momentum_density_getter | substitute_force_getter_cu | indent(8) }}
            const {{dtype}} rho_inv = {{dtype}} {1} / rho;
            {% for i in range(D) -%}
                velocity.get({{i}}u) = md_{{i}} * rho_inv;
            {% endfor %}
        }
    }
// LCOV_EXCL_STOP

    std::array<{{dtype}}, {{Q}}u> get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{Q}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::array<{{dtype}}, {{Q}}u> pop;
        thrust::copy(dev_data.begin(), dev_data.end(), pop.data());
        return pop;
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        std::array< {{dtype}}, {{Q}}u > const & pop,
        Cell const & cell )
    {
        thrust::device_vector< {{dtype}} > dev_data(pop.begin(), pop.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        CellInterval ci ( cell, cell );
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        gpu::GPUField< {{dtype}} > * velocity_field,
        gpu::GPUField< {{dtype}} > const * force_field,
        std::array< {{dtype}}, {{Q}}u > const & pop,
        Cell const & cell )
    {
        thrust::device_vector< {{dtype}} > dev_data(pop.begin(), pop.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        CellInterval ci ( cell, cell );
        auto kernel = gpu::make_kernel( kernel_set_vel );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *velocity_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void initialize(
        gpu::GPUField< {{dtype}} > * pdf_field,
        std::array< {{dtype}}, {{Q}}u > const & pop )
    {
        CellInterval ci = pdf_field->xyzSizeWithGhostLayer();
        thrust::device_vector< {{dtype}} > dev_data(pop.begin(), pop.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_broadcast );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
   }

    std::vector< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells() * {{Q}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(ci.numCells() * {{Q}}u);
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        gpu::GPUField< {{dtype}} > * velocity_field,
        gpu::GPUField< {{dtype}} > const * force_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set_vel );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *velocity_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Population

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

namespace Interpolation
{
// LCOV_EXCL_START
    /** @brief Calculate interpolation weights. */
    static __forceinline__ __device__ void calculate_weights(
        {{dtype}} const *RESTRICT const pos,
        int *RESTRICT const corner,
        {{dtype}} *RESTRICT const weights,
        uint gl)
    {
      #pragma unroll
      for (int dim = 0; dim < {{D}}; ++dim) {
        auto const fractional_index = pos[dim] - {{dtype}}{0.5};
        auto const nmp = floorf(fractional_index);
        auto const distance = fractional_index - nmp - {{dtype}}{0.5};
        corner[dim] = __{{dtype}}2int_rn(nmp) + static_cast<int>(gl);
        weights[dim * 2 + 0] = {{dtype}}{0.5} - distance;
        weights[dim * 2 + 1] = {{dtype}}{0.5} + distance;
      }
    }

    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const *RESTRICT const pos,
        {{dtype}} *RESTRICT const vel,
        uint n_pos,
        uint gl)
    {

      uint pos_index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

      vec.set({0u, 0u, 0u}, {0u, 0u, 0u});
      if (vec.isValidPosition() and pos_index < n_pos) {
        auto const array_offset = pos_index * uint({{D}}u);
        int corner[{{D}}];
        {{dtype}} weights[{{D}}][2];
        calculate_weights(pos + array_offset, corner, &weights[0][0], gl);
        #pragma unroll
        for (int i = 0; i < 2; i++) {
          auto const cx = corner[0] + i;
          auto const wx = weights[0][i];
          #pragma unroll
          for (int j = 0; j < 2; j++) {
            auto const cy = corner[1] + j;
            auto const wxy = wx * weights[1][j];
            #pragma unroll
            for (int k = 0; k < 2; k++) {
              auto const cz = corner[2] + k;
              auto const weight = wxy * weights[2][k];
              {% for cf in range(D) -%}
                vel[array_offset + {{cf}}u] += weight * vec.getNeighbor(cx, cy, cz, {{cf}}u);
              {% endfor %}
            }
          }
        }
      }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const *RESTRICT const pos,
        {{dtype}} const *RESTRICT const forces,
        uint n_pos,
        uint gl )
    {

      uint pos_index = blockIdx.y * gridDim.x * blockDim.x +
                       blockDim.x * blockIdx.x + threadIdx.x;

      vec.set({0u, 0u, 0u}, {0u, 0u, 0u});
      if (vec.isValidPosition() and pos_index < n_pos) {
        auto const array_offset = pos_index * uint({{D}}u);
        int corner[{{D}}];
        {{dtype}} weights[{{D}}][2];
        calculate_weights(pos + array_offset, corner, &weights[0][0], gl);
        #pragma unroll
        for (int i = 0; i < 2; i++) {
          auto const cx = corner[0] + i;
          auto const wx = weights[0][i];
          #pragma unroll
          for (int j = 0; j < 2; j++) {
            auto const cy = corner[1] + j;
            auto const wxy = wx * weights[1][j];
            #pragma unroll
            for (int k = 0; k < 2; k++) {
              auto const cz = corner[2] + k;
              auto const weight = wxy * weights[2][k];
              {% for cf in range(D) -%}
                atomicAdd(&vec.getNeighbor(cx, cy, cz, {{cf}}u),
                          weight * forces[array_offset + {{cf}}u]);
              {% endfor %}
            }
          }
        }
      }
    }
// LCOV_EXCL_STOP

    static dim3 calculate_dim_grid(uint const threads_x,
                                   uint const blocks_per_grid_y,
                                   uint const threads_per_block) {
      assert(threads_x >= 1u);
      assert(blocks_per_grid_y >= 1u);
      assert(threads_per_block >= 1u);
      auto const threads_y = threads_per_block * blocks_per_grid_y;
      auto const blocks_per_grid_x = (threads_x + threads_y - 1) / threads_y;
      return make_uint3(blocks_per_grid_x, blocks_per_grid_y, 1);
    }

    std::vector< {{dtype}} >
    get(
        gpu::GPUField< {{dtype}} > const *vec_field,
        std::vector< {{dtype}} > const &pos,
        uint gl )
    {
      thrust::device_vector< {{dtype}} > dev_pos(pos.begin(), pos.end());
      thrust::device_vector< {{dtype}} > dev_vel(pos.size());
      auto const dev_pos_ptr = thrust::raw_pointer_cast(dev_pos.data());
      auto const dev_vel_ptr = thrust::raw_pointer_cast(dev_vel.data());

      auto const threads_per_block = uint(64u);
      auto const n_pos = static_cast<uint>(pos.size() / {{D}}ul);
      auto const dim_grid = calculate_dim_grid(n_pos, 4u, threads_per_block);
      kernel_get<<<dim_grid, threads_per_block, 0u, nullptr>>>(
          gpu::FieldIndexing< {{dtype}} >::withGhostLayerXYZ(*vec_field, gl).gpuAccess(),
          dev_pos_ptr, dev_vel_ptr, n_pos, gl);

      std::vector< {{dtype}} > out(pos.size());
      thrust::copy(dev_vel.begin(), dev_vel.end(), out.data());
      return out;
    }

    void set(
        gpu::GPUField< {{dtype}} > const *vec_field,
        std::vector< {{dtype}} > const &pos,
        std::vector< {{dtype}} > const &forces,
        uint gl )
    {
      thrust::device_vector< {{dtype}} > dev_pos(pos.begin(), pos.end());
      thrust::device_vector< {{dtype}} > dev_for(forces.begin(), forces.end());
      auto const dev_pos_ptr = thrust::raw_pointer_cast(dev_pos.data());
      auto const dev_for_ptr = thrust::raw_pointer_cast(dev_for.data());

      auto const threads_per_block = uint(64u);
      auto const n_pos = static_cast<uint>(pos.size() / {{D}}ul);
      auto const dim_grid = calculate_dim_grid(n_pos, 4u, threads_per_block);
      kernel_set<<<dim_grid, threads_per_block, 0u, nullptr>>>(
          gpu::FieldIndexing< {{dtype}} >::withGhostLayerXYZ(*vec_field, gl).gpuAccess(),
          dev_pos_ptr, dev_for_ptr, n_pos, gl);
    }
} // namespace Interpolation

namespace Equilibrium
{
// LCOV_EXCL_START
    __device__ void kernel_set_device(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} const * RESTRICT const u,
        {{dtype}} rho )
    {
        {%if not compressible %}
        rho -= {{dtype}}(1.0);
        {%endif %}

        {% for eqTerm in equilibrium -%}
            pdf.get({{loop.index0 }}u) = {{eqTerm}};
        {% endfor -%}
    }
// LCOV_EXCL_STOP
} // namespace Equilibrium

namespace Density
{
// LCOV_EXCL_START
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT rho_out )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 1u);
        pdf.set( blockIdx, threadIdx );
        rho_out += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}}u);
            {% endfor -%}
            {{density_getters | indent(12)}}
            rho_out[0u] = rho;
        }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} const * RESTRICT rho_in )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, 1u);
        pdf.set( blockIdx, threadIdx );
        rho_in += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}}u);
            {% endfor -%}
            {{unshifted_momentum_density_getter | indent(12)}}

            // calculate current velocity (before density change)
            {{dtype}} const rho_inv = {{dtype}} {1} / rho;
            {{dtype}} const u_old[{{D}}] = { {% for i in range(D) %}momdensity_{{i}} * rho_inv{% if not loop.last %}, {% endif %}{% endfor %} };

            Equilibrium::kernel_set_device(pdf, u_old, rho_in[0u] {%if not compressible %} + {{dtype}} {1} {%endif%});
        }
    }
// LCOV_EXCL_STOP

    {{dtype}} get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(1u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        {{dtype}} rho = dev_data[0u];
        return rho;
    }

    std::vector< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(dev_data.size());
        thrust::copy(dev_data.begin(), dev_data.end(), out.begin());
        return out;
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        const {{dtype}} rho,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(1u, rho);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Density

namespace Velocity
{
// LCOV_EXCL_START
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > pdf,
        gpu::FieldAccessor< {{dtype}} > force,
        {{dtype}} * RESTRICT u_out )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
        pdf.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        u_out += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}}u);
            {% endfor -%}
            {{momentum_density_getter | substitute_force_getter_cu | indent(8) }}
            auto const rho_inv = {{dtype}} {1} / rho;
            {% for i in range(D) -%}
                u_out[{{i}}u] = md_{{i}} * rho_inv;
            {% endfor %}
        }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > pdf,
        gpu::FieldAccessor< {{dtype}} > velocity,
        gpu::FieldAccessor< {{dtype}} > force,
        {{dtype}} const * RESTRICT u_in )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
        pdf.set( blockIdx, threadIdx );
        velocity.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        u_in += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}}u);
            {% endfor -%}
            {{dtype}} const * RESTRICT const u = u_in;
            {{density_getters | indent(8)}}
            {{density_velocity_setter_macroscopic_values | substitute_force_getter_cu | indent(8)}}
            {% for i in range(D) -%}
                velocity.get({{i}}u) = u_in[{{i}}u];
            {% endfor %}
            {{dtype}} u_new[{{D}}] = { {% for i in range(D) %}u_{{i}}{% if not loop.last %}, {% endif %}{% endfor %} };

            Equilibrium::kernel_set_device(pdf, u_new, rho {%if not compressible %} + {{dtype}}(1) {%endif%});
        }
    }
// LCOV_EXCL_STOP

    Vector{{D}}< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        gpu::GPUField< {{dtype}} > const * force_field,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        Vector{{D}}< {{dtype}} > vec;
        thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
        return vec;
    }

    std::vector< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        gpu::GPUField< {{dtype}} > const * force_field,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data({{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(dev_data.size());
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        gpu::GPUField< {{dtype}} > * velocity_field,
        gpu::GPUField< {{dtype}} > const * force_field,
        Vector{{D}}< {{dtype}} > const & u,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(u.data(), u.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *velocity_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        gpu::GPUField< {{dtype}} > * velocity_field,
        gpu::GPUField< {{dtype}} > const * force_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *velocity_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Velocity

namespace Force {
// LCOV_EXCL_START
    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > pdf,
        gpu::FieldAccessor< {{dtype}} > velocity,
        gpu::FieldAccessor< {{dtype}} > force,
        {{dtype}} const * RESTRICT f_in )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
        pdf.set( blockIdx, threadIdx );
        velocity.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        f_in += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}}u);
            {% endfor -%}

            {{momentum_density_getter | substitute_force_getter_pattern("force->get\(x, ?y, ?z, ?([0-9])u?\)", "f_in[\g<1>u]") | indent(8) }}
            auto const rho_inv = {{dtype}} {1} / rho;

            {% for i in range(D) -%}
                force.get({{i}}u) = f_in[{{i}}u];
            {% endfor %}

            {% for i in range(D) -%}
                velocity.get({{i}}u) = md_{{i}} * rho_inv;
            {% endfor %}
        }
    }
// LCOV_EXCL_STOP

    void
    set( gpu::GPUField< {{dtype}} > const * pdf_field,
         gpu::GPUField< {{dtype}} > * velocity_field,
         gpu::GPUField< {{dtype}} > * force_field,
         Vector{{D}}< {{dtype}} > const & u,
         Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(u.data(), u.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *velocity_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void
    set( gpu::GPUField< {{dtype}} > const * pdf_field,
         gpu::GPUField< {{dtype}} > * velocity_field,
         gpu::GPUField< {{dtype}} > * force_field,
         std::vector< {{dtype}} > const & values,
         CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *velocity_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Force

namespace MomentumDensity
{
// LCOV_EXCL_START
    __global__ void kernel_sum(
        gpu::FieldAccessor< {{dtype}} > pdf,
        gpu::FieldAccessor< {{dtype}} > force,
        {{dtype}} * RESTRICT out )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
        pdf.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        out += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}}u);
            {% endfor -%}
            {{momentum_density_getter | substitute_force_getter_cu | indent(8) }}
            {% for i in range(D) -%}
                out[{{i}}u] += md_{{i}};
            {% endfor %}
        }
    }
// LCOV_EXCL_STOP

    Vector{{D}}< {{dtype}} > reduce(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        gpu::GPUField< {{dtype}} > const * force_field )
    {
        thrust::device_vector< {{dtype}} > dev_data({{D}}u, {{dtype}} {0});
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
            Cell cell(x, y, z);
            CellInterval ci ( cell, cell );
            auto kernel = gpu::make_kernel( kernel_sum );
            kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
            kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
            kernel.addParam( dev_data_ptr );
            kernel();
        });
        Vector{{D}}< {{dtype}} > mom({{dtype}} {0});
        thrust::copy(dev_data.begin(), dev_data.begin() + {{D}}u, mom.data());
        return mom;
    }
} // namespace MomentumDensity

namespace PressureTensor
{
// LCOV_EXCL_START
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT p_out )
    {
        auto const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D**2}}u);
        pdf.set( blockIdx, threadIdx );
        p_out += offset;
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}}u);
            {% endfor -%}
            {{second_momentum_getter | indent(12) }}
            {% for i in range(D**2) -%}
                p_out[{{i}}u] = p_{{i}};
            {% endfor %}
        }
    }
// LCOV_EXCL_STOP

    Matrix{{D}}< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{D**2}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        Matrix{{D}}< {{dtype}} > out;
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
   }

    std::vector< {{dtype}} > get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data({{D**2}}u * ci.numCells());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(dev_data.size());
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
   }
} // namespace PressureTensor


} // namespace accessor
} // namespace {{namespace}}
} // namespace walberla

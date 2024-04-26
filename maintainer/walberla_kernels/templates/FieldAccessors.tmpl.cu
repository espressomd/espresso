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
#if defined(__NVCC_DIAG_PRAGMA_SUPPORT__)
#pragma nv_diagnostic push
#pragma nv_diag_suppress 177 // unused variable
#else
#pragma push
#pragma diag_suppress 177 // unused variable
#endif
#elif defined(__clang__)
#if defined(__CUDA__)
#if defined(__CUDA_ARCH__)
// clang compiling CUDA code in device mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#else
// clang compiling CUDA code in host mode
#define RESTRICT __restrict__
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Wunused-variable"
#pragma clang diagnostic ignored "-Wunused-parameter"
#endif
#endif
#elif defined(__GNUC__) or defined(__GNUG__)
#define RESTRICT __restrict__
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-variable"
#pragma GCC diagnostic ignored "-Wunused-parameter"
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
    __global__ void kernel_get_interval(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT const pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
            {% for i in range(Q) -%}
                pop[offset + {{i}}u] = pdf.get({{i}});
            {% endfor -%}
        }
    }

    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT const pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                pop[{{i}}u] = pdf.get({{i}});
            {% endfor -%}
        }
    }

    __global__ void kernel_set_interval(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} const * RESTRICT const pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
            {% for i in range(Q) -%}
                pdf.get({{i}}) = pop[offset + {{i}}u];
            {% endfor -%}
        }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} const * RESTRICT const pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            {% for i in range(Q) -%}
                pdf.get({{i}}) = pop[{{i}}u];
            {% endfor -%}
        }
    }

    std::array<{{dtype}}, {{Q}}u> get(
        gpu::GPUField< {{dtype}} > const * pdf_field,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{Q}}u, {{dtype}} {0});
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
        thrust::device_vector< {{dtype}} > dev_data(pop.data(), pop.data() + {{Q}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        CellInterval ci ( cell, cell );
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void initialize(
        gpu::GPUField< {{dtype}} > * pdf_field,
        std::array< {{dtype}}, {{Q}}u > const & pop )
    {
        CellInterval ci = pdf_field->xyzSizeWithGhostLayer();
        thrust::device_vector< {{dtype}} > dev_data(pop.data(), pop.data() + {{Q}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
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
        auto kernel = gpu::make_kernel( kernel_get_interval );
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
        auto kernel = gpu::make_kernel( kernel_set_interval );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Population

namespace Vector
{
    __global__ void kernel_get_interval(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} * const out )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                out[offset + {{i}}u] = vec.get({{i}});
            {% endfor %}
        }
    }

    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} * const out )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            {% for i in range(D) -%}
                out[{{i}}u] = vec.get({{i}});
            {% endfor %}
        }
    }

    __global__ void kernel_set_interval(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const * RESTRICT const u )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                vec.get({{i}}) = u[offset + {{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > vec,
        const {{dtype}} * RESTRICT const u )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            {% for i in range(D) -%}
                vec.get({{i}}) = u[{{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_add_interval(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const * RESTRICT const u )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                vec.get({{i}}) += u[offset + {{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_add(
        gpu::FieldAccessor< {{dtype}} > vec,
        {{dtype}} const * RESTRICT const u )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            {% for i in range(D) -%}
                vec.get({{i}}) += u[{{i}}u];
            {% endfor %}
        }
    }

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
        auto kernel = gpu::make_kernel( kernel_set );
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
        auto kernel = gpu::make_kernel( kernel_add );
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
        auto kernel = gpu::make_kernel( kernel_get_interval );
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
        auto kernel = gpu::make_kernel( kernel_set_interval );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Vector

namespace Interpolation
{
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
                vel[array_offset + {{cf}}u] += weight * vec.getNeighbor(cx, cy, cz, {{cf}});
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
                atomicAdd(&vec.getNeighbor(cx, cy, cz, {{cf}}),
                          weight * forces[array_offset + {{cf}}u]);
              {% endfor %}
            }
          }
        }
      }
    }

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
    __device__ void kernel_set_device(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} const * RESTRICT const u,
        {{dtype}} rho )
    {
        {%if not compressible %}
        rho -= {{dtype}}(1.0);
        {%endif %}

        {% for eqTerm in equilibrium -%}
            pdf.get({{loop.index0 }}) = {{eqTerm}};
        {% endfor -%}
    }
} // namespace Equilibrium

namespace Density
{
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT const out )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, uint(1u));
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{density_getters | indent(12)}}
            out[offset] = rho;
        }
    }

    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} const * RESTRICT const rho_in )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, uint(1u));
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{unshifted_momentum_density_getter | indent(12)}}

            // calculate current velocity (before density change)
            {{dtype}} const conversion = {{dtype}}(1) / rho;
            {{dtype}} const u_old[{{D}}] = { {% for i in range(D) %}momdensity_{{i}} * conversion{% if not loop.last %}, {% endif %}{% endfor %} };

            Equilibrium::kernel_set_device(pdf, u_old, rho_in[offset] {%if not compressible %} + {{dtype}}(1) {%endif%});
        }
    }

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
        std::vector< {{dtype}} > out(ci.numCells());
        thrust::copy(dev_data.begin(), dev_data.end(), out.begin());
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
} // namespace Density

namespace Velocity
{
    __global__ void kernel_set(
        gpu::FieldAccessor< {{dtype}} > pdf,
        gpu::FieldAccessor< {{dtype}} > force,
        {{dtype}} const * RESTRICT const u_in )
    {
        pdf.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, uint({{D}}u));
            uint const bufsize = {{D}}u;
            {{dtype}} const * RESTRICT const u = u_in + bufsize * offset;
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{density_getters | indent(8)}}
            {{density_velocity_setter_macroscopic_values | substitute_force_getter_cu | indent(8)}}
            {{dtype}} u_new[{{D}}] = { {% for i in range(D) %}u_{{i}}{% if not loop.last %}, {% endif %}{% endfor %} };

            Equilibrium::kernel_set_device(pdf, u_new, rho {%if not compressible %} + {{dtype}}(1) {%endif%});
        }
    }

    void set(
        gpu::GPUField< {{dtype}} > * pdf_field,
        gpu::GPUField< {{dtype}} > * force_field,
        Vector{{D}}< {{dtype}} > const & u,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(u.data(), u.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = gpu::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( gpu::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Velocity

namespace MomentumDensity
{
    __global__ void kernel_sum(
        gpu::FieldAccessor< {{dtype}} > pdf,
        gpu::FieldAccessor< {{dtype}} > force,
        {{dtype}} * RESTRICT const out )
    {
        pdf.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            uint const bufsize = {{D}}u;
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, bufsize);
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{momentum_density_getter | substitute_force_getter_cu | indent(8) }}
            {% for i in range(D) -%}
                out[bufsize * offset + {{i}}u] += md_{{i}};
            {% endfor %}
        }
    }

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
    __global__ void kernel_get(
        gpu::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT const out )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            uint const bufsize = {{D**2}}u;
            uint const offset = getLinearIndex(blockIdx, threadIdx, gridDim, blockDim, bufsize);
            {% for i in range(Q) -%}
                {{dtype}} const f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{second_momentum_getter | indent(12) }}
            {% for i in range(D) -%}
                {% for j in range(D) -%}
                    out[bufsize * offset + {{i*D+j}}u] = p_{{i*D+j}};
                {% endfor %}
            {% endfor %}
        }
    }

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
        thrust::copy(dev_data.begin(), dev_data.begin() + {{D**2}}u, out.data());
        return out;
   }
} // namespace PressureTensor


} // namespace accessor
} // namespace {{namespace}}
} // namespace walberla

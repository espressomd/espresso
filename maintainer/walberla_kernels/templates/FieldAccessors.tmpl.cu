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

#include <cuda/FieldAccessor.h>
#include <cuda/FieldIndexing.h>
#include <cuda/GPUField.h>
#include <cuda/Kernel.h>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

#include <array>
#include <tuple>
#include <vector>

#ifdef __GNUC__
#define RESTRICT __restrict__
#elif _MSC_VER
#define RESTRICT __restrict
#else
#define RESTRICT
#endif

__device__ inline uint get_num_threads( uint3 gridDim, uint3 blockDim ) {
  return gridDim.x * gridDim.y * gridDim.z * blockDim.x * blockDim.y * blockDim.z;
}

__device__ inline uint getLinearIndexXYZF( uint3 blockIdx, uint3 threadIdx, uint3 gridDim, uint3 blockDim ) {
  auto const x = threadIdx.x;
  auto const y = blockIdx.x;
  auto const z = blockIdx.y;
  auto const f = blockIdx.z;
  auto const xSize = blockDim.x;
  auto const ySize = gridDim.x;
  auto const zSize = gridDim.y;
  return x                         +
         y * xSize                 +
         z * xSize * ySize         +
         f * xSize * ySize * zSize ;
}

__device__ inline uint getLinearIndexFZYX( uint3 blockIdx, uint3 threadIdx, uint3 gridDim, uint3 blockDim, uint fOffset ) {
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
        cuda::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT const pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
            {% for i in range(Q) -%}
                pop[offset + {{i}}u] = pdf.get({{i}});
            {% endfor -%}
        }
    }

    __global__ void kernel_get(
        cuda::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT const pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
            {% for i in range(Q) -%}
                pop[{{i}}u] = pdf.get({{i}});
            {% endfor -%}
        }
    }

    __global__ void kernel_set_interval(
        cuda::FieldAccessor< {{dtype}} > pdf,
        const {{dtype}} * RESTRICT const pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
            {% for i in range(Q) -%}
                pdf.get({{i}}) = pop[offset + {{i}}u];
            {% endfor -%}
        }
    }

    __global__ void kernel_set(
        cuda::FieldAccessor< {{dtype}} > pdf,
        const {{dtype}} * RESTRICT const pop )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{Q}}u);
            {% for i in range(Q) -%}
                pdf.get({{i}}) = pop[{{i}}u];
            {% endfor -%}
        }
    }

    std::array<{{dtype}}, {{Q}}u> get(
        cuda::GPUField< {{dtype}} > const * pdf_field,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{Q}}u, {{dtype}} {0});
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::array<{{dtype}}, {{Q}}u> pop;
        thrust::copy(dev_data.begin(), dev_data.end(), pop.data());
        return pop;
    }

    void set(
        cuda::GPUField< {{dtype}} > * pdf_field,
        std::array< {{dtype}}, {{Q}}u > const & pop,
        Cell const & cell )
    {
        thrust::device_vector< {{dtype}} > dev_data(pop.data(), pop.data() + {{Q}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        CellInterval ci ( cell, cell );
        auto kernel = cuda::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void broadcast(
        cuda::GPUField< {{dtype}} > * pdf_field,
        std::array< {{dtype}}, {{Q}}u > const & pop )
    {
        CellInterval ci = pdf_field->xyzSizeWithGhostLayer();
        thrust::device_vector< {{dtype}} > dev_data(pop.data(), pop.data() + {{Q}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
   }

    std::vector< {{dtype}} > get(
        cuda::GPUField< {{dtype}} > const * pdf_field,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells() * {{Q}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_get_interval );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(ci.numCells() * {{Q}}u);
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
    }

    void set(
        cuda::GPUField< {{dtype}} > * pdf_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_set_interval );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Population

namespace Vector
{
    __global__ void kernel_get_interval(
        cuda::FieldAccessor< {{dtype}} > vec,
        {{dtype}} * const out )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                out[offset + {{i}}u] = vec.get({{i}});
            {% endfor %}
        }
    }

    __global__ void kernel_get(
        cuda::FieldAccessor< {{dtype}} > vec,
        {{dtype}} * const out )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                out[{{i}}u] = vec.get({{i}});
            {% endfor %}
        }
    }

    __global__ void kernel_set_interval(
        cuda::FieldAccessor< {{dtype}} > vec,
        const {{dtype}} * RESTRICT const u )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                vec.get({{i}}) = u[offset + {{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_set(
        cuda::FieldAccessor< {{dtype}} > vec,
        const {{dtype}} * RESTRICT const u )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                vec.get({{i}}) = u[{{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_add_interval(
        cuda::FieldAccessor< {{dtype}} > vec,
        const {{dtype}} * RESTRICT const u )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                vec.get({{i}}) += u[offset + {{i}}u];
            {% endfor %}
        }
    }

    __global__ void kernel_add(
        cuda::FieldAccessor< {{dtype}} > vec,
        const {{dtype}} * RESTRICT const u )
    {
        vec.set( blockIdx, threadIdx );
        if (vec.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, {{D}}u);
            {% for i in range(D) -%}
                vec.get({{i}}) += u[{{i}}u];
            {% endfor %}
        }
    }

    Vector{{D}}< {{dtype}} > get(
        cuda::GPUField< {{dtype}} > const * vec_field,
        Cell const & cell)
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        Vector{{D}}< {{dtype}} > vec;
        thrust::copy(dev_data.begin(), dev_data.end(), vec.data());
        return vec;
    }

    void set(
        cuda::GPUField< {{dtype}} > * vec_field,
        Vector{{D}}< {{dtype}} > const & vec,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(vec.data(), vec.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void add(
        cuda::GPUField< {{dtype}} > * vec_field,
        Vector{{D}}< {{dtype}} > const & vec,
        Cell const &cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(vec.data(), vec.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_add );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    void broadcast(
        cuda::GPUField< {{dtype}} > * vec_field,
        Vector{{D}}< {{dtype}} > const & vec )
    {
        CellInterval ci = vec_field->xyzSizeWithGhostLayer();
        thrust::device_vector< {{dtype}} > dev_data(vec.data(), vec.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
   }

    void add_to_all(
        cuda::GPUField< {{dtype}} > * vec_field,
        Vector{{D}}< {{dtype}} > const & vec )
    {
        CellInterval ci = vec_field->xyzSizeWithGhostLayer();
        thrust::device_vector< {{dtype}} > dev_data(vec.data(), vec.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_add );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    std::vector< {{dtype}} > get(
        cuda::GPUField< {{dtype}} > const * vec_field,
        CellInterval const & ci)
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells() * {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_get_interval );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(ci.numCells() * {{D}}u);
        thrust::copy(dev_data.begin(), dev_data.end(), out.data());
        return out;
    }

    void set(
        cuda::GPUField< {{dtype}} > * vec_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_set_interval );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *vec_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Vector

namespace Equilibrium
{
    __device__ void kernel_set_device(
        cuda::FieldAccessor< {{dtype}} > pdf,
        const {{dtype}} * RESTRICT const u,
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
        cuda::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT const out )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, uint(1u));
            {% for i in range(Q) -%}
                const {{dtype}} f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{density_getters | indent(12)}}
            out[offset] = rho;
        }
    }

    __global__ void kernel_set(
        cuda::FieldAccessor< {{dtype}} > pdf,
        const {{dtype}} * RESTRICT const rho_in )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, uint(1u));
            {% for i in range(Q) -%}
                const {{dtype}} f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{unshifted_momentum_density_getter | indent(12)}}

            // calculate current velocity (before density change)
            const {{dtype}} conversion = {{dtype}}(1) / rho;
            const {{dtype}} u_old[{{D}}] = { {% for i in range(D) %}momdensity_{{i}} * conversion{% if not loop.last %}, {% endif %}{% endfor %} };

            Equilibrium::kernel_set_device(pdf, u_old, rho_in[offset] {%if not compressible %} + {{dtype}}(1) {%endif%});
        }
    }

    {{dtype}} get(
        cuda::GPUField< {{dtype}} > const * pdf_field,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(1u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        {{dtype}} rho = dev_data[0u];
        return rho;
    }

    void set(
        cuda::GPUField< {{dtype}} > * pdf_field,
        const {{dtype}} rho,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(1u, rho);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }

    std::vector< {{dtype}} > get(
        cuda::GPUField< {{dtype}} > const * pdf_field,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(ci.numCells());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( dev_data_ptr );
        kernel();
        std::vector< {{dtype}} > out(ci.numCells());
        thrust::copy(dev_data.begin(), dev_data.end(), out.begin());
        return out;
    }

    void set(
        cuda::GPUField< {{dtype}} > * pdf_field,
        std::vector< {{dtype}} > const & values,
        CellInterval const & ci )
    {
        thrust::device_vector< {{dtype}} > dev_data(values.begin(), values.end());
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Density

namespace Velocity
{
    __global__ void kernel_set(
        cuda::FieldAccessor< {{dtype}} > pdf,
        cuda::FieldAccessor< {{dtype}} > force,
        const {{dtype}} * RESTRICT const u_in )
    {
        pdf.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, uint({{D}}u));
            const uint_t bufsize = {{D}}u;
            const {{dtype}} * RESTRICT const u = u_in + bufsize * offset;
            {% for i in range(Q) -%}
                const {{dtype}} f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{density_getters | indent(8)}}
            {{density_velocity_setter_macroscopic_values | substitute_force_getter_cu | indent(8)}}
            {{dtype}} u_new[{{D}}] = { {% for i in range(D) %}u_{{i}}{% if not loop.last %}, {% endif %}{% endfor %} };

            Equilibrium::kernel_set_device(pdf, u_new, rho {%if not compressible %} + {{dtype}}(1) {%endif%});
        }
    }

    void set(
        cuda::GPUField< {{dtype}} > * pdf_field,
        cuda::GPUField< {{dtype}} > * force_field,
        Vector{{D}}< {{dtype}} > const & u,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data(u.data(), u.data() + {{D}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_set );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
        kernel.addParam( const_cast<const {{dtype}} *>(dev_data_ptr) );
        kernel();
    }
} // namespace Velocity

namespace MomentumDensity
{
    __global__ void kernel_sum(
        cuda::FieldAccessor< {{dtype}} > pdf,
        cuda::FieldAccessor< {{dtype}} > force,
        {{dtype}} * RESTRICT const out )
    {
        pdf.set( blockIdx, threadIdx );
        force.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint bufsize = {{D}}u;
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, bufsize);
            {% for i in range(Q) -%}
                const {{dtype}} f_{{i}} = pdf.get({{i}});
            {% endfor -%}
            {{momentum_density_getter | substitute_force_getter_cu | indent(8) }}
            {% for i in range(D) -%}
                out[bufsize * offset + {{i}}u] += md_{{i}};
            {% endfor %}
        }
    }

    Vector{{D}}< {{dtype}} > reduce(
        cuda::GPUField< {{dtype}} > const * pdf_field,
        cuda::GPUField< {{dtype}} > const * force_field )
    {
        thrust::device_vector< {{dtype}} > dev_data({{D}}u, {{dtype}} {0});
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        WALBERLA_FOR_ALL_CELLS_XYZ(pdf_field, {
            Cell cell(x, y, z);
            CellInterval ci ( cell, cell );
            auto kernel = cuda::make_kernel( kernel_sum );
            kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
            kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *force_field, ci ) );
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
        cuda::FieldAccessor< {{dtype}} > pdf,
        {{dtype}} * RESTRICT const out )
    {
        pdf.set( blockIdx, threadIdx );
        if (pdf.isValidPosition()) {
            const uint bufsize = {{D**2}}u;
            const uint offset = getLinearIndexFZYX(blockIdx, threadIdx, gridDim, blockDim, bufsize);
            {% for i in range(Q) -%}
                const {{dtype}} f_{{i}} = pdf.get({{i}});
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
        cuda::GPUField< {{dtype}} > const * pdf_field,
        Cell const & cell )
    {
        CellInterval ci ( cell, cell );
        thrust::device_vector< {{dtype}} > dev_data({{D**2}}u);
        auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
        auto kernel = cuda::make_kernel( kernel_get );
        kernel.addFieldIndexingParam( cuda::FieldIndexing< {{dtype}} >::interval( *pdf_field, ci ) );
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

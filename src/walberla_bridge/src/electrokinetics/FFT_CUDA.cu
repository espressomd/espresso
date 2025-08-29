/*
 * Copyright (C) 2025 The ESPResSo project
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

#include "FFT_CUDA.cuh"

#include <gpu/FieldAccessor.h>
#include <gpu/FieldIndexing.h>

#include <gpu/Kernel.h>

#include <heffte.h>
#include <heffte_backends.h>
#include <heffte_geometry.h>

#include <utils/Vector.hpp>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

namespace walberla {

template <typename FloatType> struct FFT_CUDA<FloatType>::heffte_container {
  std::shared_ptr<heffte::box3d<>> m_box_in;
  std::shared_ptr<heffte::box3d<>> m_box_out;
  std::shared_ptr<heffte::fft3d<heffte::backend::cufft>> m_fft;
  std::shared_ptr<
      heffte::fft3d<heffte::backend::cufft>::buffer_container<ComplexType>>
      m_buffer;
};

// LCOV_EXCL_START
template <typename FloatType>
__global__ void
create_greens_function(gpu::FieldAccessor<FloatType> greens_function, int x,
                       int y, int z) {
  using RealType = std::conditional<std::is_same<FloatType, float>::value,
                                    cufftReal, cufftDoubleReal>::type;
  greens_function.set(blockIdx, threadIdx);
  unsigned int index =
      greens_function.getLinearIndex(blockIdx, threadIdx, gridDim, blockDim);
  unsigned int tmp;
  unsigned int coord[3];

  coord[0] = index % x;
  tmp = index / x;
  coord[1] = tmp % y;
  coord[2] = tmp / y;
  if (index < z * y * x) {
    if (index == 0) {
      // setting 0th Fourier mode to 0 enforces charge neutrality
      greens_function.get(0u) = 0.0f;
    } else {
      constexpr RealType two_pi = 2.0f * M_PI;
      greens_function.get(0u) = -0.5f /
                                (cos(two_pi * static_cast<RealType>(coord[0]) /
                                     static_cast<RealType>(x)) +
                                 cos(two_pi * static_cast<RealType>(coord[1]) /
                                     static_cast<RealType>(y)) +
                                 cos(two_pi * static_cast<RealType>(coord[2]) /
                                     static_cast<RealType>(z)) -
                                 3.0f) /
                                static_cast<RealType>(x * y * z);
    }
  }
}

template <typename FloatType, typename ComplexType>
__global__ void
multiply_by_greens_function(gpu::FieldAccessor<ComplexType> potential,
                            gpu::FieldAccessor<FloatType> greens_function) {
  potential.set(blockIdx, threadIdx);
  greens_function.set(blockIdx, threadIdx);
  if (potential.isValidPosition() && greens_function.isValidPosition()) {
    potential.get(0u) =
        ComplexType(potential.get(0u).x * greens_function.get(0u),
                    potential.get(0u).y * greens_function.get(0u));
  }
}

template <typename FloatType>
__global__ void add_fields_with_factor(gpu::FieldAccessor<FloatType> field_out,
                                       gpu::FieldAccessor<FloatType> field_add,
                                       const FloatType factor) {
  field_out.set(blockIdx, threadIdx);
  field_add.set(blockIdx, threadIdx);
  if (field_out.isValidPosition() && field_add.isValidPosition()) {
    field_out.get(0u) += field_add.get(0u) * factor;
  }
}

template <typename FloatType>
__global__ void move_field(gpu::FieldAccessor<FloatType> dest_field,
                           gpu::FieldAccessor<FloatType> src_field) {
  dest_field.set(blockIdx, threadIdx);
  src_field.set(blockIdx, threadIdx);
  if (dest_field.isValidPosition() && src_field.isValidPosition()) {
    dest_field.get(0u) = src_field.get(0u);
    src_field.get(0u) = 0.0;
  }
}
// LCOV_EXCL_STOP

template <typename T, std::size_t N>
auto to_array(Utils::Vector<T, N> const &vec) {
  std::array<T, N> res{};
  std::copy(vec.begin(), vec.end(), res.begin());
  return res;
}

template <typename FloatType>
FFT_CUDA<FloatType>::FFT_CUDA(std::shared_ptr<LatticeWalberla> lattice,
                              double permittivity)
    : m_lattice(std::move(lattice)), m_permittivity(permittivity),
      kernel_greens(gpu::make_kernel(
          multiply_by_greens_function<FloatType, ComplexType>)),
      kernel_move_fields(gpu::make_kernel(move_field<FloatType>)) {
  m_blocks = get_lattice().get_blocks();

  m_potential_field_id = gpu::addGPUFieldToStorage<PotentialField>(
      get_lattice().get_blocks(), "potential field", 1, field::fzyx, 0, false);
  m_potential_field_with_ghosts_id = gpu::addGPUFieldToStorage<PotentialField>(
      get_lattice().get_blocks(), "potential field with ghosts", 1, field::fzyx,
      get_lattice().get_ghost_layers());
  m_greens_function_field_id = gpu::addGPUFieldToStorage<GreenFunctionField>(
      get_lattice().get_blocks(), "greens function", 1, field::fzyx, 0, false);
  m_potential_fourier_id = gpu::addGPUFieldToStorage<PotentialFourier>(
      get_lattice().get_blocks(), "fourier field", 1, field::fzyx, 0, false);
  reset_charge_field();

  heffte = std::make_shared<heffte_container>();
  auto dim = get_lattice().get_grid_dimensions();
  auto offset_vec = Utils::Vector3i({1, 1, 1});
  auto order = Utils::Vector3i({0, 1, 2});
  heffte->m_box_in = std::make_shared<heffte::box3d<>>(
      to_array(Utils::Vector3i({0, 0, 0})), to_array(dim - offset_vec),
      to_array(order));
  heffte->m_box_out = std::make_shared<heffte::box3d<>>(
      to_array(Utils::Vector3i({0, 0, 0})), to_array(dim - offset_vec),
      to_array(order));
  heffte->m_fft = std::make_shared<heffte::fft3d<heffte::backend::cufft>>(
      *(heffte->m_box_in), *(heffte->m_box_out), MPI_COMM_WORLD);
  heffte->m_buffer = std::make_shared<
      heffte::fft3d<heffte::backend::cufft>::buffer_container<ComplexType>>(
      heffte->m_fft->size_workspace());

  auto block = get_lattice().get_blocks()->getBlock(0, 0, 0);
  auto green_field =
      block->template getData<GreenFunctionField>(m_greens_function_field_id);
  auto kernel = gpu::make_kernel(create_greens_function<FloatType>);
  kernel.addFieldIndexingParam(
      gpu::FieldIndexing<FloatType>::xyz(*green_field));
  kernel.addParam(dim[0]);
  kernel.addParam(dim[1]);
  kernel.addParam(dim[2]);
  kernel();

  for (auto &block : *get_lattice().get_blocks()) {
    auto potential =
        block.template getData<PotentialField>(m_potential_field_id);
    auto potential_ghosts = block.template getData<PotentialField>(
        m_potential_field_with_ghosts_id);
    auto green =
        block.template getData<GreenFunctionField>(m_greens_function_field_id);
    auto fourier =
        block.template getData<PotentialFourier>(m_potential_fourier_id);

    kernel_greens =
        gpu::make_kernel(multiply_by_greens_function<FloatType, ComplexType>);
    kernel_greens.addFieldIndexingParam(
        gpu::FieldIndexing<ComplexType>::allInner(*fourier));
    kernel_greens.addFieldIndexingParam(
        gpu::FieldIndexing<FloatType>::allInner(*green));

    kernel_move_fields = gpu::make_kernel(move_field<FloatType>);
    kernel_move_fields.addFieldIndexingParam(
        gpu::FieldIndexing<FloatType>::xyz(*potential_ghosts));
    kernel_move_fields.addFieldIndexingParam(
        gpu::FieldIndexing<FloatType>::xyz(*potential));
  }

  m_full_communication =
      std::make_shared<FullCommunicator>(get_lattice().get_blocks());
  m_full_communication->addPackInfo(
      std::make_shared<gpu::communication::MemcpyPackInfo<PotentialField>>(
          m_potential_field_with_ghosts_id));
}

template <typename FloatType> void FFT_CUDA<FloatType>::reset_charge_field() {
  // the FFT-solver re-uses the potential field for the charge
  auto const potential_id = walberla::BlockDataID(m_potential_field_id);

  for (auto &block : *get_lattice().get_blocks()) {
    auto field = block.template getData<PotentialField>(potential_id);
    ek::accessor::Scalar::initialize(field, FloatType_c(0.0));
  }
}

template <typename FloatType>
void FFT_CUDA<FloatType>::add_fields(PotentialField *field_out,
                                     gpu::GPUField<FloatType> *field_add,
                                     FloatType factor) {
  auto kernel = gpu::make_kernel(add_fields_with_factor<FloatType>);
  kernel.addFieldIndexingParam(gpu::FieldIndexing<FloatType>::xyz(*field_out));
  kernel.addFieldIndexingParam(gpu::FieldIndexing<FloatType>::xyz(*field_add));
  kernel.addParam(factor);
  kernel();
}

template <typename FloatType>
void FFT_CUDA<FloatType>::add_charge_to_field(std::size_t id, double valency,
                                              bool is_double_precision) {
  auto const factor = FloatType_c(valency) / FloatType_c(get_permittivity());
  const auto density_id = walberla::BlockDataID(id);
  for (auto &block : *get_lattice().get_blocks()) {
    auto field = block.template getData<PotentialField>(m_potential_field_id);
    // TODO do we enforce, that the FloatType of the Charge and the
    // species density is the same?
    auto density_field =
        block.template getData<gpu::GPUField<FloatType>>(density_id);
    add_fields(field, density_field, FloatType_c(factor));
  }
}

template <typename FloatType> void FFT_CUDA<FloatType>::solve() {
  for (auto &block : *get_lattice().get_blocks()) {
    auto potential =
        block.template getData<PotentialField>(m_potential_field_id);
    auto fourier =
        block.template getData<PotentialFourier>(m_potential_fourier_id);
    FloatType *_data_potential = potential->dataAt(0, 0, 0, 0);
    ComplexType *_data_fourier = fourier->dataAt(0, 0, 0, 0);
    heffte->m_fft->forward(_data_potential, _data_fourier,
                           heffte->m_buffer->data());
    kernel_greens();
    heffte->m_fft->backward(_data_fourier, _data_potential,
                            heffte->m_buffer->data());
    kernel_move_fields();
    ghost_communication();
  }
}

template class FFT_CUDA<float>;
template class FFT_CUDA<double>;

} // namespace walberla

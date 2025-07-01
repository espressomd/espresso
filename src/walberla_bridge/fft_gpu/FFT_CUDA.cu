#include "FFT_CUDA.cuh"

#include <gpu/FieldAccessor.h>
#include <gpu/FieldIndexing.h>
// #include <gpu/GPUField.h>
#include <gpu/Kernel.h>

#include <heffte_geometry.h>

#include <utils/Vector.hpp>

#include <thrust/device_ptr.h>
#include <thrust/device_vector.h>

namespace walberla {

__device__ unsigned int getThreadIndex() {

  return blockIdx.y * gridDim.x * blockDim.x + blockDim.x * blockIdx.x +
         threadIdx.x;
}

__global__ void
create_greens_function(gpu::FieldAccessor<double> greens_function, int x, int y,
                       int z) {
  unsigned int index = getThreadIndex();
  unsigned int tmp;
  unsigned int coord[3];

  coord[0] = index % (x / 2 + 1);
  tmp = index / (y / 2 + 1);
  coord[1] = tmp % y;
  coord[2] = tmp / y;

  if (index < z * y * (x / 2 + 1)) {
    if (index == 0) {
      // setting 0th Fourier mode to 0 enforces charge neutrality
      greens_function.get(0u) = 0.0f;
    } else {
      constexpr cufftReal two_pi = 2.0f * 3.141592654f; // TODO PI
      greens_function.get(0u) = -2.0f * two_pi * 0.5f /
                                (cos(two_pi * static_cast<cufftReal>(coord[0]) /
                                     static_cast<cufftReal>(x)) +
                                 cos(two_pi * static_cast<cufftReal>(coord[1]) /
                                     static_cast<cufftReal>(y)) +
                                 cos(two_pi * static_cast<cufftReal>(coord[2]) /
                                     static_cast<cufftReal>(z)) -
                                 3.0f) /
                                static_cast<cufftReal>(x * y * z);
    }
  }
}

__global__ void
multiply_by_greens_function(gpu::FieldAccessor<cufftDoubleComplex> potential,
                            gpu::FieldAccessor<double> greens_function) {
  potential.set(blockIdx, threadIdx);
  if (potential.isValidPosition()) {
    potential.get(0u) =
        cufftDoubleComplex(potential.get(0u).x * greens_function.get(0u),
                           potential.get(0u).y * greens_function.get(0u));
  }
  // unsigned int index = fde_getThreadIndex();
  // potential[index].x *= greensfunction[index];
  // potential[index].y *= greensfunction[index];
}

template <typename T, std::size_t N>
auto to_array(Utils::Vector<T, N> const &vec) {
  std::array<T, N> res{};
  std::copy(vec.begin(), vec.end(), res.begin());
  return res;
};

template <typename FloatType>
FFT_CUDA<FloatType>::FFT_CUDA(std::shared_ptr<LatticeWalberla> lattice,
                              double permittivity)
    : m_lattice(std::move(lattice)), m_permittivity(permittivity) {
  m_blocks = get_lattice().get_blocks();

  // Vector3<uint_t> dim(m_blocks->getNumberOfXCells(),
  //                     m_blocks->getNumberOfYCells(),
  //                     m_blocks->getNumberOfZCells());
  // auto const greens = [dim](uint_t x, uint_t y, uint_t z) -> real_t {
  //   if (x == 0u && y == 0u && z == 0u)
  //     return 0.;
  //   return -0.5 /
  //           (std::cos(2. * std::numbers::pi * real_c(x) / real_c(dim[0])) +
  //           std::cos(2. * std::numbers::pi * real_c(y) / real_c(dim[1])) +
  //           std::cos(2. * std::numbers::pi * real_c(z) / real_c(dim[2])) -
  //           3.) /
  //           real_c(dim[0] * dim[1] * dim[2]);
  // };

  m_potential_field_id = gpu::addGPUFieldToStorage<PotentialField>(
      get_lattice().get_blocks(), "potential field", 1, field::fzyx,
      get_lattice().get_ghost_layers());
  m_greens_function_field_id = gpu::addGPUFieldToStorage<GreenFunctionField>(
      get_lattice().get_blocks(), "greens function", 1, field::fzyx,
      get_lattice().get_ghost_layers());
  m_potential_furier_id = gpu::addGPUFieldToStorage<PotentialFurier>(
      get_lattice().get_blocks(), "furier field", 1, field::fzyx,
      get_lattice().get_ghost_layers());
  reset_charge_field();

  auto dim = get_lattice().get_grid_dimensions();
  m_box_in = std::make_shared<heffte::box3d<>>(
      to_array(Utils::Vector3i({0, 0, 0})), to_array(dim));
  m_box_out = std::make_shared<heffte::box3d<>>(
      to_array(Utils::Vector3i({0, 0, 0})), to_array(dim));
  m_fft = std::make_shared<heffte::fft3d<heffte::backend::cufft>>(
      *m_box_in, *m_box_out, MPI_COMM_WORLD);
  m_buffer = std::make_shared<
      heffte::fft3d<heffte::backend::cufft>::buffer_container<ComplexType>>(
      m_fft->size_workspace());
  // m_fft_out =
  // std::make_shared<heffte::gpu::vector<std::complex<FloatType>>>(m_fft->size_outbox());

  auto block = get_lattice().get_blocks()->getBlock(0, 0, 0);
  auto green_field =
      block->template getData<GreenFunctionField>(m_greens_function_field_id);
  auto kernel = gpu::make_kernel(create_greens_function);
  kernel.addFieldIndexingParam(
      gpu::FieldIndexing<FloatType>::xyz(*green_field));
  kernel.addParam<int>(dim[0]);
  kernel.addParam<int>(dim[1]);
  kernel.addParam<int>(dim[2]);
  kernel();

  // m_full_communication =
  //     std::make_shared<FullCommunicator>(get_lattice().get_blocks());
  // m_full_communication->addPackInfo(
  //     std::make_shared<field::communication::PackInfo<PotentialField>>(
  //         m_potential_field_id));
}

template <typename FloatType> void FFT_CUDA<FloatType>::reset_charge_field() {
  // the FFT-solver re-uses the potential field for the charge
  auto const potential_id = walberla::BlockDataID(get_potential_field_id());

  for (auto &block : *get_lattice().get_blocks()) {
    auto field = block.template getData<PotentialField>(potential_id);
    ek::accessor::Scalar::initialize(field, FloatType_c(0.0));
  }
}

template <typename FloatType>
void FFT_CUDA<FloatType>::add_charge_to_field(std::size_t id, double valency,
                                              bool is_double_precision) {
  auto const factor = FloatType_c(valency) / FloatType_c(get_permittivity());
  // the FFT-solver re-uses the potential field for the charge
  const auto charge_id = walberla::BlockDataID(get_potential_field_id());
  // const auto density_id = walberla::BlockDataID(id);
  for (auto &block : *get_lattice().get_blocks()) {
    auto field = block.template getData<PotentialField>(charge_id);
    ek::accessor::Scalar::initialize(field, FloatType_c(0.0));
  }
}

template <typename FloatType> void FFT_CUDA<FloatType>::solve() {
  for (auto &block : *get_lattice().get_blocks()) {
    auto potential =
        block.template getData<PotentialField>(m_potential_field_id);
    auto green =
        block.template getData<GreenFunctionField>(m_greens_function_field_id);
    auto furier =
        block.template getData<PotentialFurier>(m_potential_furier_id);
    FloatType *_data_potential = potential->dataAt(-1, -1, -1, 0);
    ComplexType *_data_furier = furier->dataAt(-1, -1, -1, 0);
    const int64_t _stride_0 = int64_t(potential->xStride());
    const int64_t _stride_1 = int64_t(potential->yStride());
    const int64_t _stride_2 = int64_t(potential->zStride());
    const int64_t _stride_3 = int64_t(1 * int64_t(potential->fStride()));
    // thrust::device_vector<double> dev_data(1u);
    // auto const dev_data_ptr = thrust::raw_pointer_cast(dev_data.data());
    auto kernel = gpu::make_kernel(multiply_by_greens_function);
    kernel.addFieldIndexingParam(gpu::FieldIndexing<ComplexType>::xyz(*furier));
    kernel.addFieldIndexingParam(gpu::FieldIndexing<FloatType>::xyz(*green));
    // kernel.addParam(dev_data_ptr);
    kernel();
    // dim3 dim_grid = calculate_dim_grid(
    //   static_cast<unsigned>(parameters.dim_z * parameters.dim_y *
    //                         (parameters.dim_x / 2 + 1)),
    //   4, threads_per_block);

    m_fft->forward(_data_potential, _data_furier); //, m_buffer->data());
    kernel();
    // KERNELCALL(multiply_by_greens_function, dim_grid, threads_per_block,
    // _data_potential);
    m_fft->backward(_data_furier, _data_potential); //, m_buffer->data());

    ghost_communication();
  }
}

template class FFT_CUDA<float>;
template class FFT_CUDA<double>;

} // namespace walberla

#pragma once

// #include <blockforest/communication/UniformBufferedScheme.h>
#include <domain_decomposition/BlockDataID.h>
// #include <fft/Fft.h>
// #include <field/AddToStorage.h>
// #include <field/GhostLayerField.h>
// #include <field/communication/PackInfo.h>
#include <gpu/AddGPUFieldToStorage.h>
#include <gpu/GPUField.h>
#include <gpu/communication/UniformGPUScheme.h>
// #include <gpu/Kernel.h>
// #include <gpu/FieldAccessor.h>
// #include <gpu/FieldIndexing.h>
#include <heffte.h>
#include <heffte_backends.h>
#include <stencil/D3Q27.h>
#include <walberla_bridge/LatticeWalberla.hpp>
// TODO dirty
#include "../src/electrokinetics/generated_kernels/EK_FieldAccessors_double_precision_CUDA.cuh"
#include "../src/electrokinetics/generated_kernels/EK_FieldAccessors_single_precision_CUDA.cuh"

#include <cmath>
#include <cstddef>
#include <memory>
#include <numbers>
#include <utility>

namespace walberla {

template <typename FloatType> class FFT_CUDA {
private:
  template <typename T> FloatType FloatType_c(T t) {
    return numeric_cast<FloatType>(t);
  }

  std::shared_ptr<LatticeWalberla> m_lattice;
  double m_permittivity;

  domain_decomposition::BlockDataID m_potential_field_id;
  domain_decomposition::BlockDataID m_greens_function_field_id;
  domain_decomposition::BlockDataID m_potential_furier_id;

  using ComplexType = std::conditional<std::is_same<FloatType, float>::value,
                                       cufftComplex, cufftDoubleComplex>::type;
  using PotentialField = gpu::GPUField<FloatType>;
  using GreenFunctionField = gpu::GPUField<FloatType>;
  using PotentialFurier = gpu::GPUField<ComplexType>;

  std::shared_ptr<heffte::box3d<>> m_box_in;
  std::shared_ptr<heffte::box3d<>> m_box_out;
  std::shared_ptr<heffte::fft3d<heffte::backend::cufft>> m_fft;
  std::shared_ptr<
      heffte::fft3d<heffte::backend::cufft>::buffer_container<ComplexType>>
      m_buffer;
  // std::shared_ptr<heffte::gpu::vector<std::complex<FloatType>>> m_fft_out;
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

  // using FullCommunicator = gpu::communication::UniformGPUScheme<
  //     typename stencil::D3Q27>;
  // std::shared_ptr<FullCommunicator> m_full_communication;

public:
  // FFT_CUDA() = default;
  FFT_CUDA(std::shared_ptr<LatticeWalberla> lattice, double permittivity);
  ~FFT_CUDA() = default;

  void reset_charge_field();

  void add_charge_to_field(std::size_t id, double valency,
                           bool is_double_precision);

  [[nodiscard]] std::size_t get_potential_field_id() const noexcept {
    return static_cast<std::size_t>(m_potential_field_id);
  }

  void solve();

  void set_permittivity(double permittivity) noexcept {
    m_permittivity = permittivity;
  }

  [[nodiscard]] double get_permittivity() const noexcept {
    return m_permittivity;
  }

  [[nodiscard]] auto const &get_lattice() const noexcept { return *m_lattice; }

private:
  void add_fields(PotentialField *field_out,
                  gpu::GPUField<FloatType> *field_add, FloatType factor);
  void ghost_communication() {} //(*m_full_communication)(); }
};

} // namespace walberla

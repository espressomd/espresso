#include "LBWalberlaImpl.hpp"
#ifdef __AVX2__
#include "generated_kernels/MRTLatticeModelAvx.h"
#define LatticeModelName lbm::MRTLatticeModelAvx
#else
#include "generated_kernels/MRTLatticeModel.h"
#define LatticeModelName lbm::MRTLatticeModel
#endif

namespace walberla {
class LBWalberlaD3Q19MRT : public LBWalberlaImpl<LatticeModelName> {
  using LatticeModel = LatticeModelName;

public:
  void construct_lattice_model(double viscosity) {
    const real_t omega = 2 / (6 * real_c(viscosity) + 1);
    const real_t magic_number = real_c(3.) / real_c(16.);
    const real_t omega_2 =
        (4 - 2 * omega) / (4 * magic_number * omega + 2 - omega);
    m_lattice_model = std::make_shared<LatticeModel>(
        LatticeModel(m_last_applied_force_field_id,
                     omega,   // bulk
                     omega,   // even
                     omega_2, // odd
                     omega)); // shear
  };
  void set_viscosity(double viscosity) override {
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    const real_t omega = 2 / (6 * real_c(viscosity) + 1);
    const real_t magic_number = real_c(3.) / real_c(16.);
    const real_t omega_2 =
        (4 - 2 * omega) / (4 * magic_number * omega + 2 - omega);
    lm->omega_shear_ = omega;
    lm->omega_odd_ = omega_2;
    lm->omega_even_ = omega;
    lm->omega_bulk_ = omega;
    on_lattice_model_change();
  };
  double get_viscosity() const override {
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return (2 - lm->omega_shear_) / (6 * lm->omega_shear_);
  };
  LBWalberlaD3Q19MRT(double viscosity, double density,
                     const Utils::Vector3i &grid_dimensions,
                     const Utils::Vector3i &node_grid, int n_ghost_layers)
      : LBWalberlaImpl(viscosity, grid_dimensions, node_grid, n_ghost_layers) {
    construct_lattice_model(viscosity);
    setup_with_valid_lattice_model(density);
  };
};

} // namespace walberla

#undef LatticeModelName

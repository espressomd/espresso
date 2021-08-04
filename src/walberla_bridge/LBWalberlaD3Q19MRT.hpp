#include "LBWalberlaImpl.hpp"
#include "relaxation_rates.hpp"
#ifdef __AVX2__
#include "generated_kernels/CollideSweepAVX.h"
#define CollisionModelName walberla::pystencils::CollideSweepAVX
#include "generated_kernels/MRTLatticeModelAvx.h"
#define LatticeModelName lbm::MRTLatticeModelAvx
#else
#include "generated_kernels/CollideSweep.h"
#define CollisionModelName walberla::pystencils::CollideSweep
#include "generated_kernels/MRTLatticeModel.h"
#define LatticeModelName lbm::MRTLatticeModel
#endif
#define CollisionModelFactoryName walberla::pystencils::CollideSweepFactory

namespace walberla {
class LBWalberlaD3Q19MRT
    : public LBWalberlaImpl<LatticeModelName, CollisionModelName, CollisionModelFactoryName> {
  using LatticeModel = LatticeModelName;

public:
  void construct_lattice_model(double viscosity) {
    const real_t omega = shear_mode_relaxation_rate(viscosity);
    const real_t omega_odd = odd_mode_relaxation_rate(omega);
    m_lattice_model = std::make_shared<LatticeModel>(
        LatticeModel(m_last_applied_force_field_id,
                     omega,     // bulk
                     omega,     // even
                     omega_odd, // odd
                     omega));   // shear
  };
  void set_viscosity(double viscosity) override {
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    const real_t omega = shear_mode_relaxation_rate(viscosity);
    const real_t omega_odd = odd_mode_relaxation_rate(omega);
    lm->omega_shear_ = omega;
    lm->omega_odd_ = omega_odd;
    lm->omega_even_ = omega;
    lm->omega_bulk_ = omega;
    on_lattice_model_change();
  };
  double get_viscosity() const override {
    auto *lm = dynamic_cast<LatticeModel *>(m_lattice_model.get());
    return viscosity_from_shear_relaxation_rate(lm->omega_shear_);
  };
  LBWalberlaD3Q19MRT(double viscosity, double density,
                     const Utils::Vector3i &grid_dimensions,
                     const Utils::Vector3i &node_grid, int n_ghost_layers)
      : LBWalberlaImpl(viscosity, grid_dimensions, node_grid, n_ghost_layers) {
    construct_lattice_model(viscosity);
    setup_with_valid_lattice_model(density, 0u, 0u);
  };
};

} // namespace walberla

#undef LatticeModelName
#undef CollisionModelName
#undef CollisionModelFactoryName

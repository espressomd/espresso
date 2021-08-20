#include "LBWalberlaImpl.hpp"
#include "relaxation_rates.hpp"

#ifdef __AVX__
#include "generated_kernels/CollideSweepThermalizedAVX.h"
#define CollisionModelName walberla::pystencils::CollideSweepThermalizedAVX
#include "generated_kernels/FluctuatingMRTLatticeModelAvx.h"
#define LatticeModelName lbm::FluctuatingMRTLatticeModelAvx
#else
#include "generated_kernels/CollideSweepThermalized.h"
#define CollisionModelName walberla::pystencils::CollideSweepThermalized
#include "generated_kernels/FluctuatingMRTLatticeModel.h"
#define LatticeModelName lbm::FluctuatingMRTLatticeModel
#endif

namespace walberla {
class LBWalberlaD3Q19FluctuatingMRT
    : public LBWalberlaImpl<LatticeModelName, CollisionModelName> {

  using LatticeModel = LatticeModelName;

public:
  void construct_lattice_model(double viscosity, double kT, unsigned int seed) {
    const real_t omega = shear_mode_relaxation_rate(viscosity);
    const real_t omega_odd = odd_mode_relaxation_rate(omega);
    m_lattice_model = std::make_shared<LatticeModel>(
        LatticeModel(m_last_applied_force_field_id, real_c(kT),
                     omega,     // bulk
                     omega,     // even
                     omega_odd, // odd
                     omega,     // shear
                     seed,      // RNG seed
                     0          // time_step
                     ));
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
  LBWalberlaD3Q19FluctuatingMRT(double viscosity, double density,
                                const Utils::Vector3i &grid_dimensions,
                                const Utils::Vector3i &node_grid,
                                int n_ghost_layers, double kT,
                                unsigned int seed)
      : LBWalberlaImpl(viscosity, grid_dimensions, node_grid, n_ghost_layers,
                       kT, seed) {
    construct_lattice_model(viscosity, kT, seed);
    setup_with_valid_lattice_model(density, seed, 0u);
  };
};

} // namespace walberla
#undef LatticeModelName
#undef CollisionModelName

#include "LBWalberlaImpl.hpp"

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
  LBWalberlaD3Q19FluctuatingMRT(double viscosity, double density,
                                const Utils::Vector3i &grid_dimensions,
                                const Utils::Vector3i &node_grid,
                                int n_ghost_layers, double kT,
                                unsigned int seed)
      : LBWalberlaImpl(viscosity, grid_dimensions, node_grid, n_ghost_layers,
                       kT, seed) {
    m_lattice_model = std::make_shared<LatticeModel>(
        m_last_applied_force_field_id, real_c(kT), -1., -1., -1., -1., seed, 0);
    setup_with_valid_lattice_model(density, seed, 0u);
  };
};

} // namespace walberla
#undef LatticeModelName
#undef CollisionModelName

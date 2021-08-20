#include "LBWalberlaImpl.hpp"
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

namespace walberla {
class LBWalberlaD3Q19MRT
    : public LBWalberlaImpl<LatticeModelName, CollisionModelName> {
  using LatticeModel = LatticeModelName;

public:
  LBWalberlaD3Q19MRT(double viscosity, double density,
                     const Utils::Vector3i &grid_dimensions,
                     const Utils::Vector3i &node_grid, int n_ghost_layers,
                     double kT, unsigned int seed)
      : LBWalberlaImpl(viscosity, grid_dimensions, node_grid, n_ghost_layers,
                       kT, seed) {
    m_lattice_model = std::make_shared<LatticeModel>(
        m_last_applied_force_field_id, -1., -1., -1., -1.);
    setup_with_valid_lattice_model(density, 0u, 0u);
  };
};

} // namespace walberla

#undef LatticeModelName
#undef CollisionModelName

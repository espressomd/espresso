#include "LbWalberla_impl.hpp"

namespace walberla {
using LatticeModelD3Q19TRT = walberla::lbm::D3Q19<
    walberla::lbm::collision_model::TRT, false,
    walberla::lbm::force_model::GuoField<
        walberla::GhostLayerField<walberla::Vector3<walberla::real_t>, 1>>>;

class LbWalberlaD3Q19TRT : public LbWalberla<LatticeModelD3Q19TRT> {
public:
  using LbWalberla::LbWalberla;
  void construct_lattice_model(double viscosity) {
    m_lattice_model =
        std::make_shared<LatticeModelD3Q19TRT>(LatticeModelD3Q19TRT(
            lbm::collision_model::TRT::constructWithMagicNumber(
                lbm::collision_model::omegaFromViscosity((real_t)viscosity)),
            lbm::force_model::GuoField<VectorField>(
                m_last_applied_force_field_id)));
  };
  void set_viscosity(double viscosity) override {
    get_lattice_model()->collisionModel().resetWithMagicNumber(
        lbm::collision_model::omegaFromViscosity((real_t)viscosity));
  };
};

} // namespace walberla

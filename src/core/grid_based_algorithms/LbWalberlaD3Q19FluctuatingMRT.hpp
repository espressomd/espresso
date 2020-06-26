#include "FluctuatingMRT_LatticeModel.h"
#include "LbWalberla_impl.hpp"

namespace walberla {
class LbWalberlaD3Q19FluctuatingMRT
    : public LbWalberla<lbm::FluctuatingMRT_LatticeModel> {
public:
  void construct_lattice_model(double viscosity, double kT, unsigned int seed) {
    const double omega = 2 / (6* viscosity + 1);
    const double magic_number = 3. / 16.;
    const double omega_2 =
        (4 - 2 * omega) / (4 * magic_number * omega + 2 - omega);
    m_lattice_model = std::make_shared<lbm::FluctuatingMRT_LatticeModel>(
        lbm::FluctuatingMRT_LatticeModel(m_last_applied_force_field_id, 
                                         omega_2,omega_2,omega_2,omega, 1,kT,seed));

  };
  void set_viscosity(double viscosity) override {
    lbm::FluctuatingMRT_LatticeModel *lm =
        dynamic_cast<lbm::FluctuatingMRT_LatticeModel *>(m_lattice_model.get());
    const double omega = 2 / (6 * viscosity + 1);
    const double magic_number = 3. / 16.;
    const double omega_2 =
        (4 - 2 * omega) / (4 * magic_number * omega + 2 - omega);
    lm->omega_shear_ = omega;
    lm->omega_odd_ = omega_2;
    lm->omega_even_ =omega_2;
    lm->omega_bulk_ = omega_2;
  };
  double get_viscosity() const override {
    lbm::FluctuatingMRT_LatticeModel *lm =
        dynamic_cast<lbm::FluctuatingMRT_LatticeModel *>(m_lattice_model.get());
    return (2 - lm->omega_shear_) / (6 * lm->omega_shear_);
  };
  LbWalberlaD3Q19FluctuatingMRT(double viscosity, double density, double agrid, double tau,
                     const Utils::Vector3d &box_dimensions,
                     const Utils::Vector3i &node_grid, int n_ghost_layers,
                     double kT, unsigned int seed)
      : LbWalberla(viscosity, density, agrid, tau, box_dimensions, node_grid,
                   n_ghost_layers) {
    construct_lattice_model(viscosity, kT, seed);
    setup_with_valid_lattice_model();
  };
  void integrate() override {
    m_time_loop->singleStep();
    lbm::FluctuatingMRT_LatticeModel *lm =
        dynamic_cast<lbm::FluctuatingMRT_LatticeModel *>(m_lattice_model.get());
    lm->time_step_ +=1;
    printf("Time step is %d\n",lm->time_step_);
  };
};

} // namespace walberla

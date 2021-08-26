#include "ekin_walberla_init.hpp"

#include "EKinWalberlaBase.hpp"
#include "EKinWalberlaImpl.hpp"

#include <memory>

EKinWalberlaBase<double> *
new_ekin_walberla(const walberla::WalberlaBlockForest *blockforest,
                  double diffusion, double kT, double valency, double density,
                  bool advection, bool friction_coupling) {
  return new walberla::EKinWalberlaImpl(
      walberla::EKinWalberlaImpl(blockforest, diffusion, kT, valency, density,
                                 advection, friction_coupling));
}

namespace walberla {
std::shared_ptr<EKinWalberlaBase<double>>
new_ek_walberla(const walberla::WalberlaBlockForest *blockforest,
                double diffusion, double kT, double valency, double density,
                bool advection, bool friction_coupling) {
  return std::make_shared<EKinWalberlaImpl<13, double>>(
      blockforest, diffusion, kT, valency, density, advection,
      friction_coupling);
}
} // namespace walberla
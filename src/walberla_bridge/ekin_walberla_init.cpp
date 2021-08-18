#include "ekin_walberla_init.hpp"

#include "EKinWalberlaBase.hpp"
#include "EKinWalberlaImpl.hpp"
#include "PoissonSolver/PoissonSolver.hpp"

#include <memory>

EKinWalberlaBase<double> *
new_ekin_walberla(const walberla::WalberlaBlockForest *blockforest,
                  double diffusion, double kT, double valency, double density) {
  return new walberla::EKinWalberlaImpl(
      walberla::EKinWalberlaImpl(blockforest, diffusion, kT, valency, density));
}

namespace walberla {
std::shared_ptr<EKinWalberlaBase<double>>
new_ek_walberla(const walberla::WalberlaBlockForest *blockforest,
                double diffusion, double kT, double valency, double density) {
  return std::make_shared<EKinWalberlaImpl<13, double>>(blockforest, diffusion,
                                                        kT, valency, density);
}

std::unique_ptr<EKWalberlaCharge<double>>
new_ek_charge(const walberla::WalberlaBlockForest *blockforest,
              const walberla::PoissonSolver<double> *poissonsolver) {
  return std::make_unique<EKWalberlaCharge<double>>(blockforest, poissonsolver);
}
} // namespace walberla
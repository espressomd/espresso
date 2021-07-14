#include "ekin_walberla_init.hpp"

#include "EKinWalberlaBase.hpp"
#include "EKinWalberlaImpl.hpp"

EKinWalberlaBase<double> *
new_ekin_walberla(const walberla::WalberlaBlockForest *blockforest,
                  double diffusion, double kT, double density) {
  return new walberla::EKinWalberlaImpl(
      walberla::EKinWalberlaImpl(blockforest, diffusion, kT, density));
}
#ifndef ESPRESSO_EKIN_WALBERLA_INIT_HPP
#define ESPRESSO_EKIN_WALBERLA_INIT_HPP

#include "EKinWalberlaBase.hpp"

EKinWalberlaBase<double> *
new_ekin_walberla(const walberla::WalberlaBlockForest *blockforest,
                  double diffusion, double kT, double density);

namespace walberla {
std::shared_ptr<EKinWalberlaBase<double>>
new_ek_walberla(const walberla::WalberlaBlockForest *blockforest,
                double diffusion, double kT, double density);
}

#endif // ESPRESSO_EKIN_WALBERLA_INIT_HPP

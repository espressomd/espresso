#ifndef ESPRESSO_EKIN_WALBERLA_INIT_HPP
#define ESPRESSO_EKIN_WALBERLA_INIT_HPP

#include "EKinWalberlaBase.hpp"

#include <memory>

EKinWalberlaBase<double> *
new_ekin_walberla(const walberla::WalberlaBlockForest *blockforest,
                  double diffusion, double kT, double valency, double density);

namespace walberla {
std::shared_ptr<EKinWalberlaBase<double>>
new_ek_walberla(const walberla::WalberlaBlockForest *blockforest,
                double diffusion, double kT, double valency, double density);
} // namespace walberla

#endif // ESPRESSO_EKIN_WALBERLA_INIT_HPP

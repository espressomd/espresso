#ifndef ESPRESSO_WALBERLA_BLOCKFOREST_HPP
#define ESPRESSO_WALBERLA_BLOCKFOREST_HPP

#include <memory>

#include "WalberlaBlockForest.hpp"

namespace walberla {
struct WalberlaBlockForestParams {
  explicit WalberlaBlockForestParams(double agrid) : m_agrid(agrid) {}
  [[nodiscard]] double get_agrid() const { return m_agrid; };

private:
  double m_agrid;
};
} // namespace walberla

std::shared_ptr<walberla::WalberlaBlockForest> get_walberla_blockforest();
const walberla::WalberlaBlockForestParams *get_walberla_blockforest_params();

void mpi_init_walberla_blockforest(const Utils::Vector3d &box_size,
                                   double agrid, int n_ghost_layers);

#endif // ESPRESSO_WALBERLA_BLOCKFOREST_HPP

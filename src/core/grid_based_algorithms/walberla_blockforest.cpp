#include "walberla_blockforest.hpp"
#include "MpiCallbacks.hpp"
#include "communication.hpp"
#include "grid.hpp"

namespace {
/** @brief Access the per-MPI-node Walberla instance */
std::shared_ptr<walberla::WalberlaBlockForest> walberla_blockforest;
} // namespace

std::shared_ptr<walberla::WalberlaBlockForest> get_walberla_blockforest() {
  if (!walberla_blockforest) {
    throw std::runtime_error(
        "Attempted access to uninitialized BlockForest instance.");
  }
  return walberla_blockforest;
}

void init_walberla_blockforest_local(const Utils::Vector3i &grid_dimensions,
                                     const Utils::Vector3i &node_grid,
                                     int n_ghost_layers) {
  walberla_blockforest = std::make_shared<walberla::WalberlaBlockForest>(
      grid_dimensions, node_grid, n_ghost_layers);
}

REGISTER_CALLBACK(init_walberla_blockforest_local)

void mpi_init_walberla_blockforest(const Utils::Vector3d &box_size,
                                   double agrid, int n_ghost_layers) {
  const Utils::Vector3i grid_dimensions{
      static_cast<int>(std::round(box_size[0] / agrid)),
      static_cast<int>(std::round(box_size[1] / agrid)),
      static_cast<int>(std::round(box_size[2] / agrid))};
  for (int i : {0, 1, 2}) {
    if (fabs(grid_dimensions[i] * agrid - box_size[i]) / box_size[i] >
        std::numeric_limits<double>::epsilon()) {
      throw std::runtime_error(
          "Box length not commensurate with agrid in direction " +
          std::to_string(i) + " length " + std::to_string(box_size[i]) +
          " agrid " + std::to_string(agrid));
    }
  }
  mpi_call_all(init_walberla_blockforest_local, grid_dimensions, node_grid,
               n_ghost_layers);
}
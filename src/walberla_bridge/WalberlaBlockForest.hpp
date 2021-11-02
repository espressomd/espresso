#ifndef ESPRESSO_WALBERLABLOCKFOREST_HPP
#define ESPRESSO_WALBERLABLOCKFOREST_HPP

#include "blockforest/Initialization.h"
#include "blockforest/StructuredBlockForest.h"

#include "walberla_utils.hpp"

namespace walberla {

/** Class that runs and controls the BlockForest in walberla
 */
class WalberlaBlockForest {
private:
  /** Member variables */
  Utils::Vector3i m_grid_dimensions;
  int m_n_ghost_layers;

  /** Block forest */
  std::shared_ptr<blockforest::StructuredBlockForest> m_blocks;

public:
  WalberlaBlockForest(const Utils::Vector3i &grid_dimensions,
                      const Utils::Vector3i &node_grid, int n_ghost_layers)
      : m_grid_dimensions{grid_dimensions}, m_n_ghost_layers{n_ghost_layers} {

    if (m_n_ghost_layers <= 0)
      throw std::runtime_error("At least one ghost layer must be used");
    for (int i : {0, 1, 2}) {
      if (m_grid_dimensions[i] % node_grid[i] != 0) {
        throw std::runtime_error(
            "LB grid dimensions and mpi node grid are not compatible.");
      }
    }

    m_blocks = blockforest::createUniformBlockGrid(
        uint_c(node_grid[0]), // blocks in x direction
        uint_c(node_grid[1]), // blocks in y direction
        uint_c(node_grid[2]), // blocks in z direction
        uint_c(m_grid_dimensions[0] /
               node_grid[0]), // number of cells per block in x direction
        uint_c(m_grid_dimensions[1] /
               node_grid[1]), // number of cells per block in y direction
        uint_c(m_grid_dimensions[2] /
               node_grid[2]), // number of cells per block in z direction
        1,                    // Lattice constant
        uint_c(node_grid[0]), uint_c(node_grid[1]),
        uint_c(node_grid[2]), // cpus per direction
        true, true, true);
  };

  // Grid, domain, halo
  [[nodiscard]] auto get_ghost_layers() const { return m_n_ghost_layers; };
  [[nodiscard]] auto get_grid_dimensions() const { return m_grid_dimensions; }
  [[nodiscard]] auto get_block() const { return m_blocks; }
  [[nodiscard]] std::pair<Utils::Vector3d, Utils::Vector3d>
  get_local_domain() const {
    // We only have one block per mpi rank
    assert(++(m_blocks->begin()) == m_blocks->end());

    auto const ab = m_blocks->begin()->getAABB();
    return {to_vector3d(ab.min()), to_vector3d(ab.max())};
  };

  [[nodiscard]] bool node_in_local_domain(const Utils::Vector3i &node) const {
    // Note: Lattice constant =1, cell centers offset by .5
    return get_block_and_cell(node, false, m_blocks, get_ghost_layers()) !=
           boost::none;
  };
  [[nodiscard]] bool node_in_local_halo(const Utils::Vector3i &node) const {
    return get_block_and_cell(node, true, m_blocks, get_ghost_layers()) !=
           boost::none;
  };
  [[nodiscard]] bool pos_in_local_domain(const Utils::Vector3d &pos) const {
    return ::walberla::get_block(pos, false, m_blocks, get_ghost_layers()) !=
           nullptr;
  };
  [[nodiscard]] bool pos_in_local_halo(const Utils::Vector3d &pos) const {
    return ::walberla::get_block(pos, true, m_blocks, get_ghost_layers()) !=
           nullptr;
  };
};
} // namespace walberla

#endif // ESPRESSO_WALBERLABLOCKFOREST_HPP

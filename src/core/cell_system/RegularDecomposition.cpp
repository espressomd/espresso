/*
 * Copyright (C) 2010-2026 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#include "cell_system/RegularDecomposition.hpp"

#include "cell_system/Cell.hpp"
#include "ghosts/HaloPlanValidator.hpp"
#include "ghosts/mark_boundary_cells.hpp"

#include "communication.hpp"
#include "error_handling/RuntimeErrorStream.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>
#include <utils/mpi/cart_comm.hpp>
#include <utils/mpi/sendrecv.hpp>

#include <boost/container/flat_set.hpp>
#include <boost/mpi/collectives/all_reduce.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/mpi/request.hpp>
#include <boost/range/numeric.hpp>

#include <Kokkos_Core.hpp>

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <exception>
#include <functional>
#include <initializer_list>
#include <iterator>
#include <map>
#include <mutex>
#include <numeric>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

int RegularDecomposition::position_to_cell_index(
    Utils::Vector3d const &pos) const {
  Utils::Vector3i cpos;

  for (auto i = 0u; i < 3u; i++) {
    cpos[i] = static_cast<int>(std::floor(pos[i] * inv_cell_size[i])) + 1 -
              cell_offset[i];

    /* particles outside our box. Still take them if
       nonperiodic boundary. We also accept the particle if we are at
       the box boundary, and the particle is within the box. In this case
       the particle belongs here and could otherwise potentially be dismissed
       due to rounding errors. */
    if (cpos[i] < 1) {
      if ((!m_box.periodic(i) or (pos[i] >= m_box.length()[i])) and
          m_local_box.boundary()[2u * i])
        cpos[i] = 1;
      else
        return -1;
    } else if (cpos[i] > cell_grid[i]) {
      if ((!m_box.periodic(i) or (pos[i] < m_box.length()[i])) and
          m_local_box.boundary()[2u * i + 1u])
        cpos[i] = cell_grid[i];
      else
        return -1;
    }
  }

  return get_linear_index(cpos, ghost_cell_grid);
}

void RegularDecomposition::move_if_local(
    ParticleList &src, ParticleList &rest,
    std::vector<ParticleChange> &modified_cells) {
  for (auto &part : src) {
    auto target_cell = position_to_cell(part.pos());

    if (target_cell) {
      target_cell->particles().insert(std::move(part));
      modified_cells.emplace_back(ModifiedList{target_cell->particles()});
    } else {
      rest.insert(std::move(part));
    }
  }

  src.clear();
}

void RegularDecomposition::move_left_or_right(ParticleList &src,
                                              ParticleList &left,
                                              ParticleList &right,
                                              int dir) const {
  auto const is_open_boundary_left = m_local_box.boundary()[2 * dir] != 0;
  auto const is_open_boundary_right = m_local_box.boundary()[2 * dir + 1] != 0;
  auto const can_move_left = m_box.periodic(dir) or not is_open_boundary_left;
  auto const can_move_right = m_box.periodic(dir) or not is_open_boundary_right;
  auto const my_left = m_local_box.my_left()[dir];
  auto const my_right = m_local_box.my_right()[dir];
  for (auto it = src.begin(); it != src.end();) {
    auto const pos = it->pos()[dir];
    if (m_box.get_mi_coord(pos, my_right, dir) >= 0. and can_move_right) {
      right.insert(std::move(*it));
      it = src.erase(it);
    } else if (m_box.get_mi_coord(pos, my_left, dir) < 0. and can_move_left) {
      left.insert(std::move(*it));
      it = src.erase(it);
    } else {
      ++it;
    }
  }
}

void RegularDecomposition::exchange_neighbors(
    ParticleList &pl, std::vector<ParticleChange> &modified_cells) {
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(m_comm);
  static ParticleList send_buf_l, send_buf_r, recv_buf_l, recv_buf_r;

  for (int dir = 0; dir < 3; dir++) {
    /* Single node direction, no action needed. */
    if (Utils::Mpi::cart_get<3>(m_comm).dims[dir] == 1) {
      continue;
      /* In this (common) case left and right neighbors are
         the same, and we need only one communication */
    }
    if (Utils::Mpi::cart_get<3>(m_comm).dims[dir] == 2) {
      move_left_or_right(pl, send_buf_l, send_buf_l, dir);

      Utils::Mpi::sendrecv(m_comm, node_neighbors[2 * dir], 0, send_buf_l,
                           node_neighbors[2 * dir], 0, recv_buf_l);

      send_buf_l.clear();
    } else {
      using boost::mpi::request;
      using Utils::Mpi::isendrecv;

      move_left_or_right(pl, send_buf_l, send_buf_r, dir);

      auto req_l = isendrecv(m_comm, node_neighbors[2 * dir], 0, send_buf_l,
                             node_neighbors[2 * dir], 0, recv_buf_l);
      auto req_r = isendrecv(m_comm, node_neighbors[2 * dir + 1], 0, send_buf_r,
                             node_neighbors[2 * dir + 1], 0, recv_buf_r);

      std::array<request, 4> reqs{{req_l[0], req_l[1], req_r[0], req_r[1]}};
      boost::mpi::wait_all(reqs.begin(), reqs.end());

      send_buf_l.clear();
      send_buf_r.clear();
    }

    move_if_local(recv_buf_l, pl, modified_cells);
    move_if_local(recv_buf_r, pl, modified_cells);
  }
}

/**
 * @brief Fold coordinates to box and reset the old position.
 */
static void fold_and_reset(Particle &p, BoxGeometry const &box_geo) {
  box_geo.fold_position(p.pos(), p.image_box());

  p.pos_at_last_verlet_update() = p.pos();
}

void RegularDecomposition::resort(bool global,
                                  std::vector<ParticleChange> &diff) {
  ParticleList displaced_parts;

  auto const cells_span = local_cells();
  auto const n_cells = cells_span.size();

  /* Remove a misplaced particle from its cell and hand it to its target
   * cell (or the displaced list when it left the local domain), recording
   * the changes. Shared by the serial and the two-phase parallel sweep. */
  auto const apply_move = [&](ParticleList &parts, Particle &&p,
                              Cell *target_cell) {
    diff.emplace_back(ModifiedList{parts});
    /* Particle is not local */
    if (target_cell == nullptr) {
      diff.emplace_back(RemovedParticle{p.id()});
      displaced_parts.insert(std::move(p));
    }
    /* Particle belongs on this node but is in the wrong cell. */
    else {
      target_cell->particles().insert(std::move(p));
      diff.emplace_back(ModifiedList{target_cell->particles()});
    }
  };

  if (Kokkos::DefaultHostExecutionSpace().concurrency() == 1) {
    /* Single-threaded rank: one-pass sweep without the bookkeeping overhead
     * of the two-phase version below. */
    for (auto *const c : cells_span) {
      for (auto it = c->particles().begin(); it != c->particles().end();) {
        fold_and_reset(*it, m_box);

        auto *const target_cell = particle_to_cell(*it);

        /* Particle is in place */
        if (target_cell == c) {
          std::advance(it, 1);
          continue;
        }

        auto p = std::move(*it);
        it = c->particles().erase(it);
        apply_move(c->particles(), std::move(p), target_cell);
      }
    }
  } else {
    using exec_space = Kokkos::DefaultHostExecutionSpace;
    /* Phase 1 (parallel): fold every particle position and classify it
     * against its target cell. Only particle-local state and this cell's own
     * move list are written, so cells can be swept concurrently.
     * fold_and_reset() throws on image-box overflow; exceptions must not
     * escape the parallel region, so the first error is captured and
     * rethrown afterwards. */
    struct Move {
      unsigned index;
      Cell *target;
    };
    std::vector<std::vector<Move>> moves(n_cells);
    std::mutex fold_error_mutex;
    std::string fold_error_msg;
    Kokkos::RangePolicy<exec_space> policy(std::size_t{0}, n_cells);
    Kokkos::parallel_for(
        "RegularDecomposition::resort::classify", policy,
        [&](std::size_t const ci) {
          auto *cell = cells_span[ci];
          unsigned index = 0u;
          for (auto &p : cell->particles()) {
            try {
              fold_and_reset(p, m_box);
            } catch (std::exception const &err) {
              std::lock_guard<std::mutex> guard{fold_error_mutex};
              if (fold_error_msg.empty()) {
                fold_error_msg = err.what();
              }
              return;
            }
            if (auto *const target = particle_to_cell(p); target != cell) {
              moves[ci].emplace_back(index, target);
            }
            ++index;
          }
        });
    Kokkos::fence();
    if (not fold_error_msg.empty()) {
      throw std::runtime_error(fold_error_msg);
    }

    /* Phase 2 (serial): apply the moves. This replays the serial sweep
     * exactly: ParticleList::erase() swaps the last element into the erased
     * slot, so a slot-to-original-index map is maintained to look up the
     * phase-1 classification of swapped-in elements. Cells without moves
     * (the vast majority) are skipped entirely. */
    std::vector<int> slot;
    std::vector<Cell *> target_of;
    for (std::size_t ci = 0; ci < n_cells; ++ci) {
      if (moves[ci].empty()) {
        continue;
      }
      auto *const c = cells_span[ci];
      auto &parts = c->particles();
      auto const n = static_cast<int>(parts.size());
      target_of.assign(n, c); // target == own cell: particle is in place
      for (auto const &move : moves[ci]) {
        target_of[move.index] = move.target;
      }
      slot.resize(n);
      std::iota(slot.begin(), slot.end(), 0);
      int i = 0;
      int end = n;
      while (i < end) {
        auto *const target_cell = target_of[slot[i]];

        /* Particle is in place */
        if (target_cell == c) {
          ++i;
          continue;
        }

        auto p = std::move(*(parts.begin() + i));
        parts.erase(parts.begin() + i); // swaps the last element into slot i
        slot[i] = slot[--end];
        apply_move(parts, std::move(p), target_cell);
      }
    }
  }

  if (global) {
    auto const grid = Utils::Mpi::cart_get<3>(m_comm).dims;
    /* Worst case we need grid - 1 rounds per direction.
     * This correctly implies that if there is only one node,
     * no action should be taken. */
    int rounds_left = grid[0] + grid[1] + grid[2] - 3;
    for (; rounds_left > 0; rounds_left--) {
      exchange_neighbors(displaced_parts, diff);

      auto left_over = boost::mpi::all_reduce(m_comm, displaced_parts.size(),
                                              std::plus<std::size_t>());

      if (left_over == 0) {
        break;
      }
    }
  } else {
    exchange_neighbors(displaced_parts, diff);
  }

  if (not displaced_parts.empty()) {
    auto sort_cell = local_cells()[0];

    for (auto &part : displaced_parts) {
      runtimeErrorMsg() << "Particle " << part.id() << " moved more "
                        << "than one local box length in one timestep";
      sort_cell->particles().insert(std::move(part));

      diff.emplace_back(ModifiedList{sort_cell->particles()});
    }
  }
}

void RegularDecomposition::mark_cells() {
  m_local_cells.clear();
  m_ghost_cells.clear();

  int cnt_c = 0;
  for (int o = 0; o < ghost_cell_grid[2]; o++)
    for (int n = 0; n < ghost_cell_grid[1]; n++)
      for (int m = 0; m < ghost_cell_grid[0]; m++) {
        if ((m > 0 && m < ghost_cell_grid[0] - 1 && n > 0 &&
             n < ghost_cell_grid[1] - 1 && o > 0 && o < ghost_cell_grid[2] - 1))
          m_local_cells.push_back(&cells.at(cnt_c++));
        else
          m_ghost_cells.push_back(&cells.at(cnt_c++));
      }
}

Utils::Vector3d RegularDecomposition::max_cutoff() const {
  auto dir_max_range = [this](unsigned int i) {
    return std::min(0.5 * m_box.length()[i], m_local_box.length()[i]);
  };

  return {dir_max_range(0u), dir_max_range(1u), dir_max_range(2u)};
}

Utils::Vector3d RegularDecomposition::max_range() const { return cell_size; }
int RegularDecomposition::calc_processor_min_num_cells() const {
  /* the minimal number of cells can be lower if there are at least two nodes
     serving a direction,
     since this also ensures that the cell size is at most half the box
     length. However, if there is only one processor for a direction, there
     have to be at least two cells for this direction. */
  return boost::accumulate(Utils::Mpi::cart_get<3>(m_comm).dims, 1,
                           [](int n_cells, int grid) {
                             return (grid == 1) ? 2 * n_cells : n_cells;
                           });
}

void RegularDecomposition::create_cell_grid(double range) {
  auto const cart_info = Utils::Mpi::cart_get<3>(m_comm);

  int n_local_cells;
  auto cell_range = Utils::Vector3d::broadcast(range);
  auto const min_num_cells = calc_processor_min_num_cells();

  if (range <= 0.) {
    /* this is the non-interacting case */
    auto const cells_per_dir =
        static_cast<int>(std::ceil(std::cbrt(min_num_cells)));

    cell_grid = Utils::Vector3i::broadcast(cells_per_dir);
    n_local_cells = Utils::product(cell_grid);
  } else {
    /* Calculate initial cell grid */
    auto const &local_box_l = m_local_box.length();
    auto const volume = Utils::product(local_box_l);
    auto const scale = std::cbrt(RegularDecomposition::max_num_cells / volume);

    for (auto i = 0u; i < 3u; i++) {
      /* this is at least 1 */
      cell_grid[i] = static_cast<int>(std::ceil(local_box_l[i] * scale));
      cell_range[i] = local_box_l[i] / static_cast<double>(cell_grid[i]);

      if (cell_range[i] < range) {
        /* ok, too many cells for this direction, set to minimum */
        cell_grid[i] = static_cast<int>(std::floor(local_box_l[i] / range));
        if (cell_grid[i] < 1) {
          runtimeErrorMsg()
              << "interaction range " << range << " in direction " << i
              << " is larger than the local box size " << local_box_l[i];
          cell_grid[i] = 1;
        }
        cell_range[i] = local_box_l[i] / static_cast<double>(cell_grid[i]);
      }
    }

    /* It may be necessary to asymmetrically assign the scaling to the
       coordinates, which the above approach will not do.
       For a symmetric box, it gives a symmetric result. Here we correct that.
       */
    for (;;) {
      n_local_cells = Utils::product(cell_grid);

      /* done */
      if (n_local_cells <= RegularDecomposition::max_num_cells)
        break;

      /* find coordinate with the smallest cell range */
      auto min_ind = 0u;
      auto min_size = cell_range[0];

      for (auto i = 1u; i < 3u; ++i) {
        if (cell_grid[i] > 1 and cell_range[i] < min_size) {
          min_ind = i;
          min_size = cell_range[i];
        }
      }

      cell_grid[min_ind]--;
      cell_range[min_ind] = m_local_box.length()[min_ind] / cell_grid[min_ind];
    }

    /* sanity check */
    if (n_local_cells < min_num_cells) {
      runtimeErrorMsg() << "number of cells " << n_local_cells
                        << " is smaller than minimum " << min_num_cells
                        << ": either interaction range is too large for "
                        << "the current skin (range=" << range << ", "
                        << "half_local_box_l=[" << local_box_l / 2. << "]) "
                        << "or min_num_cells too large";
    }
  }

  if (n_local_cells > RegularDecomposition::max_num_cells) {
    runtimeErrorMsg() << "no suitable cell grid found";
  }

  auto const node_pos = cart_info.coords;

  /* now set all dependent variables */
  int new_cells = 1;
  for (auto i = 0u; i < 3u; i++) {
    ghost_cell_grid[i] = cell_grid[i] + 2;
    new_cells *= ghost_cell_grid[i];
    cell_size[i] = m_local_box.length()[i] / static_cast<double>(cell_grid[i]);
    inv_cell_size[i] = 1.0 / cell_size[i];
    cell_offset[i] = node_pos[i] * cell_grid[i];
  }

  /* allocate cell array and cell pointer arrays */
  cells.clear();
  cells.resize(static_cast<unsigned int>(new_cells));
  m_local_cells.resize(n_local_cells);
  m_ghost_cells.resize(new_cells - n_local_cells);
}

template <class K, class Comparator> auto make_flat_set(Comparator &&comp) {
  return boost::container::flat_set<K, std::remove_reference_t<Comparator>>(
      std::forward<Comparator>(comp));
}

void RegularDecomposition::init_cell_interactions() {

  // Note: the global index for physical cells is 0-based.
  // I.e., a global index of -1 refers to a ghost cell.
  auto const halo = Utils::Vector3i{1, 1, 1}; // number of ghost layers
  auto const cart_info = Utils::Mpi::cart_get<3>(m_comm);
  // 3D index of the MPI rank in the Cartesian grid of MPI ranks
  auto const &node_pos = cart_info.coords;
  // size of the Cartesian grid of MPI ranks
  auto const &node_grid = ::communicator.node_grid;
  auto const global_halo_offset = hadamard_product(node_pos, cell_grid) - halo;
  // MD cell index of lower halo layer on this MPI rank
  auto const global_size = hadamard_product(node_grid, cell_grid);

  // is a cell at the system boundary in the given coord
  auto const at_boundary = [&global_size](int coord, Utils::Vector3i cell_idx) {
    return (cell_idx[coord] == 0 or cell_idx[coord] == global_size[coord] - 1);
  };

  // For the fully connected feature (cells that don't share at least a corner)
  // only apply if one cell is a ghost cell (i.e. connections across the
  // periodic boundary.
  auto const fcb_is_inner_connection = [&global_size, this](Utils::Vector3i a,
                                                            Utils::Vector3i b) {
    if (fully_connected_boundary()) {
      auto const [fc_normal, fc_dir] = *fully_connected_boundary();
      auto const involves_ghost_cell =
          (a[fc_normal] == -1 or a[fc_normal] == global_size[fc_normal] or
           b[fc_normal] == -1 or b[fc_normal] == global_size[fc_normal]);
      if (not involves_ghost_cell) {
        // check if cells do not share at least a corner
        return std::abs((a - b)[fc_dir]) > 1;
      }
    }
    return false;
  };

  /* Translate a node local index (relative to the origin of the local grid)
   * to a global index. */
  auto global_index =
      [&](Utils::Vector3i const &local_index) -> Utils::Vector3i {
    return (global_halo_offset + local_index);
  };

  /* Linear index in the global cell grid. */
  auto folded_linear_index = [&](Utils::Vector3i const &global_index) {
    auto const folded_index = (global_index + global_size) % global_size;

    return get_linear_index(folded_index, global_size);
  };

  /* Translate a global index into a local one */
  auto local_index =
      [&](Utils::Vector3i const &global_index) -> Utils::Vector3i {
    return (global_index - global_halo_offset);
  };

  // sanity checks
  if (fully_connected_boundary()) {
    auto const [fc_normal, fc_dir] = *fully_connected_boundary();
    if (fc_normal == fc_dir) {
      throw std::domain_error("fully_connected_boundary normal and connection "
                              "coordinates need to differ.");
    }
    if (node_grid[fc_dir] != 1) {
      throw std::runtime_error(
          "The MPI nodegrid must be 1 in the fully connected direction.");
    }
    if (not m_box.periodic(fc_normal)) {
      throw std::runtime_error(
          "The fully connected boundary requires periodicity in the "
          "boundary normal direction.");
    }
  }

  /* We only consider local cells (e.g. not halo cells), which
   * span the range [(1,1,1), cell_grid) in local coordinates. */
  auto const start = global_index(Utils::Vector3i{1, 1, 1});
  auto const end = start + cell_grid;

  bool one_mpi_rank = m_comm.size() == 1;

  /* loop all local cells */
  for (int o = start[2]; o < end[2]; o++)
    for (int n = start[1]; n < end[1]; n++)
      for (int m = start[0]; m < end[0]; m++) {
        /* next-nearest neighbors in every direction */
        Utils::Vector3i lower_index = {m - 1, n - 1, o - 1};
        Utils::Vector3i upper_index = {m + 1, n + 1, o + 1};

        /* In the fully connected case, we consider all cells
        * in the direction as neighbors, not only the nearest ones.
        //         */
        if (fully_connected_boundary()) {
          auto const [fc_boundary, fc_direction] = *fully_connected_boundary();

          // Fully connected is only needed at the box surface
          if (at_boundary(fc_boundary, {m, n, o})) {
            lower_index[fc_direction] = -1;
            upper_index[fc_direction] = global_size[fc_direction];
          }
        }

        /* In non-periodic directions, the halo needs not
         * be considered. */
        for (auto i = 0u; i < 3u; i++) {
          if (not m_box.periodic(i)) {
            lower_index[i] = std::max(0, lower_index[i]);
            upper_index[i] = std::min(global_size[i] - 1, upper_index[i]);
          }
        }

        /* Unique set of neighbors, cells are compared by their linear
         * index in the global cell grid. */
        auto neighbors = make_flat_set<Utils::Vector3i>(
            [&](Utils::Vector3i const &a, Utils::Vector3i const &b) {
              return folded_linear_index(a) < folded_linear_index(b);
            });

        /* Collect neighbors */
        for (int p = lower_index[2]; p <= upper_index[2]; p++)
          for (int q = lower_index[1]; q <= upper_index[1]; q++)
            for (int r = lower_index[0]; r <= upper_index[0]; r++) {
              if (fully_connected_boundary()) {
                // Avoid fully connecting the boundary layer and the
                // next INNER layer
                if (fcb_is_inner_connection({m, n, o}, {r, q, p}))
                  continue;
              }
              neighbors.insert(Utils::Vector3i{r, q, p});
            }

        /* Red-black partition by global index. */
        auto const ind1 = folded_linear_index({m, n, o});

        std::vector<Cell *> red_neighbors;
        std::vector<Cell *> black_neighbors;

        /* If we are running on a single MPI rank, it is not necessary to use
         * ghost cells. Instead of adding a ghost cell as neighbor,
         * we directly connect to the corresponding
         * physical cell across the periodic boundary */
        for (auto &neighbor : neighbors) {
          if (one_mpi_rank) {
            for (auto coord : {0u, 1u, 2u}) {
              if (neighbor[coord] == -1) {
                neighbor[coord] += cell_grid[coord];
              } else if (neighbor[coord] == cell_grid[coord]) {
                neighbor[coord] -= cell_grid[coord];
              }
            }
          }
          auto const ind2 = folded_linear_index(neighbor);
          /* Exclude cell itself */
          if (ind1 == ind2)
            continue;

          auto cell = &cells.at(
              get_linear_index(local_index(neighbor), ghost_cell_grid));

          // Divide red and black neighbors
          if (ind2 > ind1) {
            red_neighbors.push_back(cell);
          } else {
            black_neighbors.push_back(cell);
          }
        }

        // Assign neighbors to the cell
        cells[get_linear_index(local_index({m, n, o}), ghost_cell_grid)]
            .m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);
      }
}

GhostComm::HaloPlan RegularDecomposition::make_halo_plan() {
  using GhostComm::HaloPlan;
  using GhostComm::LocalComm;
  using GhostComm::NeighborComm;
  using GhostComm::SendRegion;

  HaloPlan plan;
  plan.comm = m_comm;

  // Match the legacy communicator: on a single MPI rank there are no ghost
  // cells to fill. The cell neighbourships are set up (see
  // init_cell_interactions) so that cells across periodic boundaries are
  // connected directly, so the plan stays empty.
  if (m_comm.size() == 1)
    return plan;

  auto const cart_info = Utils::Mpi::cart_get<3>(m_comm);
  auto const &node_pos = cart_info.coords;
  auto const &node_grid = ::communicator.node_grid;
  // Total number of MD cells along each axis across all ranks.
  auto const global_size = hadamard_product(node_grid, cell_grid);
  // Global (0-based) cell index of this rank's first *local* cell, i.e. of
  // ghost-grid coordinate (1,1,1).
  auto const global_origin = hadamard_product(node_pos, cell_grid);

  // Cartesian rank owning the cell at (folded) global cell coordinate.
  // The Cartesian communicator is always fully periodic (see
  // Communicator::init_comm_cart), so this is well-defined for every offset.
  auto const owner_of = [&](Utils::Vector3i const &global_cell) {
    auto const owner_coords =
        hadamard_division((global_cell + global_size) % global_size, cell_grid);
    return Utils::Mpi::cart_rank<3>(m_comm, owner_coords);
  };

  // Deterministic ordering key shared by both ranks of a peer pair: the linear
  // index of the *real* cell (in the global cell grid) that a ghost mirrors.
  auto const global_key = [&](Utils::Vector3i const &global_cell) {
    auto const folded = (global_cell + global_size) % global_size;
    return Utils::get_linear_index(folded, global_size);
  };

  // Pointer to the particle list of the cell at ghost-grid coordinate c.
  auto const list_at = [this](Utils::Vector3i const &c) -> ParticleList * {
    return &cells
                .at(static_cast<std::size_t>(
                    Utils::get_linear_index(c, ghost_cell_grid)))
                .particles();
  };

  // Per-peer accumulators. A "recv" pair maps one of our ghost cells to the
  // global index of the real cell (on the peer) it mirrors. A "send" pair maps
  // one of our real cells to its own global index (the peer will receive it
  // into a ghost). Sorting both lists by their key makes recv[k] line up with
  // peer.send[k] without exchanging any index arrays.
  struct PeerBucket {
    std::vector<std::pair<int, ParticleList *>> recv; // (key, our ghost)
    std::vector<std::pair<int, ParticleList *>> send; // (key, our real cell)
  };
  std::map<int, PeerBucket> peers;
  std::vector<std::pair<int, LocalComm>> local; // (key, self-ghost copy)

  auto const this_rank = m_comm.rank();
  auto const one = Utils::Vector3i{1, 1, 1};

  // Enumerate every ghost cell exactly once (any ghost-grid coordinate with a
  // component in the halo, i.e. == 0 or == cell_grid+1). This guarantees each
  // ghost is a recv/dst target exactly once, regardless of how small the node
  // grid is (dims of 1 or 2 collapse several stencil directions onto the same
  // peer, so a local-cell x offset enumeration would double-count).
  //
  // For each ghost we record a *matched pair*:
  //   * recv: this ghost, keyed by the global index of the real cell (on the
  //     peer) it mirrors -- so it lines up with the peer's send of that cell.
  //   * send: our own boundary real cell that the peer mirrors as its ghost in
  //     the opposite direction, keyed by that real cell's own global index --
  //     so it lines up with the peer's recv. The recv<->send pairing is a
  //     bijection, which keeps send.size() == recv.size() per peer.
  for (int gz = 0; gz < ghost_cell_grid[2]; ++gz) {
    for (int gy = 0; gy < ghost_cell_grid[1]; ++gy) {
      for (int gx = 0; gx < ghost_cell_grid[0]; ++gx) {
        Utils::Vector3i const nc{gx, gy, gz};
        // Direction sign of the halo crossing (0 if interior in a dim).
        Utils::Vector3i side{};
        bool is_ghost = false;
        for (auto d = 0u; d < 3u; ++d) {
          if (nc[d] == 0) {
            side[d] = -1;
            is_ghost = true;
          } else if (nc[d] == cell_grid[d] + 1) {
            side[d] = +1;
            is_ghost = true;
          }
        }
        if (not is_ghost)
          continue; // interior (local) cell

        // Global cell mirrored by this ghost and its owning peer.
        auto const ghost_global = global_origin + (nc - one);
        auto const peer = owner_of(ghost_global);
        auto const recv_key = global_key(ghost_global);

        if (peer == this_rank) {
          // Periodic self-ghost (a node-grid dim equals 1): copy the matching
          // local cell straight into the ghost (replaces GHOST_LOCL).
          auto const src_coord = ((ghost_global + global_size) % global_size) -
                                 global_origin + one;
          local.emplace_back(recv_key,
                             LocalComm{list_at(src_coord), list_at(nc), {}});
          continue;
        }

        peers[peer].recv.emplace_back(recv_key, list_at(nc));

        // Dual send cell: our boundary real cell mirrored by the peer's ghost
        // in the opposite direction. Snap crossing dims to the near boundary,
        // keep tangential dims aligned with the ghost.
        auto mc = nc;
        for (auto d = 0u; d < 3u; ++d) {
          if (side[d] == -1)
            mc[d] = 1;
          else if (side[d] == +1)
            mc[d] = cell_grid[d];
        }
        auto const send_global = global_origin + (mc - one);
        peers[peer].send.emplace_back(global_key(send_global), list_at(mc));
      }
    }
  }

  // Emit one NeighborComm per peer, with send/recv sorted by their shared key.
  auto const by_key = [](auto const &a, auto const &b) {
    return a.first < b.first;
  };
  for (auto &[peer, bucket] : peers) {
    std::ranges::sort(bucket.recv, by_key);
    std::ranges::sort(bucket.send, by_key);
    NeighborComm nc;
    nc.peer = peer;
    nc.recv.reserve(bucket.recv.size());
    for (auto const &[key, cell] : bucket.recv)
      nc.recv.push_back(cell);
    nc.send.reserve(bucket.send.size());
    for (auto const &[key, cell] : bucket.send)
      nc.send.push_back(SendRegion{cell, {}});
    plan.neighbors.push_back(std::move(nc));
  }

  // Sort the self-copies deterministically too (not required, but keeps the
  // plan reproducible run to run).
  std::ranges::sort(
      local, [](auto const &a, auto const &b) { return a.first < b.first; });
  plan.local.reserve(local.size());
  for (auto &[key, lc] : local)
    plan.local.push_back(lc);

  return plan;
}

RegularDecomposition::RegularDecomposition(
    boost::mpi::communicator comm, double range, BoxGeometry const &box_geo,
    LocalBox const &local_geo,
    std::optional<std::pair<int, int>> fully_connected)
    : m_comm(std::move(comm)), m_box(box_geo), m_local_box(local_geo),
      m_fully_connected_boundary(std::move(fully_connected)) {

  /* set up new regular decomposition cell structure */
  create_cell_grid(range);

  /* setup cell neighbors */
  init_cell_interactions();

  /* mark local and ghost cells */
  mark_cells();

  /* build the topology-agnostic direct-neighbor halo plan */
  m_halo_plan = make_halo_plan();

  /* Classify local cells as interior or boundary.
   *
   * Rule (a): any neighbor that is a ghost cell -> boundary [base rule].
   * Rule (b): any neighbor relation that crosses a periodic box boundary
   *   must also make the cell boundary, even when both cells are local.
   *   This happens on a single MPI rank (node_grid[i]==1, periodic[i]):
   *   init_cell_interactions() wires the first and last local layer along
   *   axis i directly without going through ghost cells, so the base rule
   *   misses those periodic wrap-around neighbours.
   *
   * Implementation: precise pair-predicate that fires iff the two cells
   * sit in the "first layer ↔ last layer" pair along a wrap axis.
   * Ghost-grid coordinate of a cell: unpack the column-major linear index
   *   idx = a + G[0]*(b + G[1]*c)  (G = ghost_cell_grid).
   * First local layer along axis i: ghost coord == 1.
   * Last local layer along axis i: ghost coord == cell_grid[i].
   */
  auto const &node_grid = ::communicator.node_grid;
  auto const idx_of = [this](Cell const *c) {
    return static_cast<int>(c - cells.data());
  };
  auto const ghost_coord_of = [this](int idx) -> Utils::Vector3i {
    int const a = idx % ghost_cell_grid[0];
    int const bc = idx / ghost_cell_grid[0];
    int const b = bc % ghost_cell_grid[1];
    int const c = bc / ghost_cell_grid[1];
    return {a, b, c};
  };
  // Which axes wrap locally (the local domain spans the whole box along a
  // periodic axis)?  Precomputed by value so the predicate lambda does not
  // capture node_grid (AppleClang rejects that non-odr-use capture with
  // -Werror,-Wunused-lambda-capture).
  std::array<bool, 3> wrap_axis;
  for (int i = 0; i < 3; ++i)
    wrap_axis[i] = (node_grid[i] == 1) && m_box.periodic(i);
  // Build the predicate only when there is at least one wrap axis; otherwise
  // pass nullptr (no overhead in the inner loop of mark_boundary_cells).
  bool const has_wrap_axis = wrap_axis[0] || wrap_axis[1] || wrap_axis[2];
  std::function<bool(Cell const *, Cell const *)> wrap_pred;
  if (has_wrap_axis) {
    wrap_pred = [this, wrap_axis, idx_of, ghost_coord_of](
                    Cell const *a_cell, Cell const *b_cell) -> bool {
      auto const a_coord = ghost_coord_of(idx_of(a_cell));
      auto const b_coord = ghost_coord_of(idx_of(b_cell));
      for (int i = 0; i < 3; ++i) {
        if (wrap_axis[i]) {
          bool const a_first = (a_coord[i] == 1);
          bool const a_last = (a_coord[i] == cell_grid[i]);
          bool const b_first = (b_coord[i] == 1);
          bool const b_last = (b_coord[i] == cell_grid[i]);
          if ((a_first && b_last) || (a_last && b_first))
            return true;
        }
      }
      return false;
    };
  }
  GhostComm::mark_boundary_cells(local_cells(), ghost_cells(), wrap_pred);

  /* Degenerate case: node_grid[i]==1, periodic[i], cell_grid[i]==1.
   *
   * When cell_grid[i]==1 on a wrap axis, the single local cell layer along
   * that axis has itself as the only periodic neighbour (both ends fold to
   * the same global index).  init_cell_interactions() therefore excludes
   * the self-pair (ind1==ind2 guard at line ~552), leaving neighbors().all()
   * empty along that axis.  The wrap_predicate above is never called for
   * that cell, so mark_boundary_cells() leaves it interior — wrong.
   *
   * Correct interpretation: the cell interacts with itself across the
   * periodic boundary, so it is by definition wrap-adjacent and must be
   * boundary.  The simplest safe fix: if any wrap axis has cell_grid[i]==1,
   * every local cell is boundary (they all span that axis, so all are
   * wrap-adjacent).
   */
  for (int i = 0; i < 3; ++i) {
    if (wrap_axis[i] && cell_grid[i] == 1) {
      for (Cell *c : local_cells())
        c->m_is_boundary = true;
      break; // one degenerate axis is enough to force all-boundary
    }
  }

  /* Plan-membership pass (source 2): mark every local cell that the plan
   * exports as a send source.  The geometric rules above (source 1) miss
   * plan shapes such as Lees-Edwards fully-connected boundaries and ELC
   * periodicity-change paths, where boundary cells are determined by the
   * plan topology rather than ghost-cell adjacency.  Both sources are
   * complementary and must both be applied. */
  GhostComm::mark_plan_cells_boundary(m_halo_plan, local_cells());

#ifdef ESPRESSO_ADDITIONAL_CHECKS
  assert(GhostComm::report_violations(
      GhostComm::validate_halo_plan(m_halo_plan, local_cells(), ghost_cells()),
      "RegularDecomposition"));
  // NOTE: validate_halo_plan_symmetry is NOT called here.
  // During checkpoint loading, decompositions are transiently rebuilt while
  // maximal_cutoff is rank-divergent (ranks may have different cell grids for a
  // brief window before the next consistent rebuild).  The transient plan is
  // never used — it is immediately replaced — so the asymmetry is harmless.
  // A construction-time collective all_to_all inside a ctor is also dangerous:
  // if one rank aborts the others block forever in the collective.
  // Symmetry is instead validated at FIRST USE of the plan in
  // halo_exchange_start (see GhostComm::halo_exchange_start in
  // HaloExchange.cpp).
#endif
}

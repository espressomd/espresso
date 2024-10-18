// clang-format: off
/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#include "cell_system/RegularSingleNodeDecompositoin.hpp"

#include "cell_system/Cell.hpp"

#include "error_handling/RuntimeErrorStream.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/Vector.hpp>
#include <utils/index.hpp>
#include <utils/mpi/cart_comm.hpp>
#include <utils/mpi/sendrecv.hpp>

#include <boost/container/flat_set.hpp>
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/communicator.hpp>
#include <boost/range/algorithm/reverse.hpp>
#include <boost/range/numeric.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <functional>
#include <iterator>
#include <utility>
#include <vector>

/** Returns pointer to the cell which corresponds to the position if the
 *  position is in the node's spatial domain otherwise a nullptr.
 */
Cell *
RegularSingleNodeDecompositoin::position_to_cell(const Utils::Vector3d &pos) {
  Utils::Vector3i cpos;

  for (int i = 0; i < 3; i++) {
    cpos[i] = static_cast<int>(std::floor(pos[i] * inv_cell_size[i]));
    /* particles outside our box. Still take them if
       nonperiodic boundary. We also accept the particle if we are at
       the box boundary, and the particle is within the box. In this case
       the particle belongs here and could otherwise potentially be dismissed
       due to rounding errors. */
    if (cpos[i] < 0) {
      if ((!m_box.periodic(i) or (pos[i] >= m_box.length()[i])) &&
          m_local_box.boundary()[2 * i])
        cpos[i] = 0;
      else
        return nullptr;
    } else if (cpos[i] >= cell_grid[i]) {
      if ((!m_box.periodic(i) or (pos[i] < m_box.length()[i])) &&
          m_local_box.boundary()[2 * i])
        cpos[i] = cell_grid[i] - 1;
      else
        return nullptr;
    }
  }

  auto const ind = get_linear_index(cpos, cell_grid);
  return &(cells.at(ind));
}

void RegularSingleNodeDecompositoin::move_locally(
    ParticleList &src, std::vector<ParticleChange> &modified_cells) {
  for (auto &part : src) {
    auto target_cell = position_to_cell(part.pos());
    if (not target_cell) {
      thorw std::runtime_error("Target cell could not be determined during "
                               "resort. This shoudl not happen. Bug. ");
    }

    target_cell->particles().insert(std::move(part));
    modified_cells.emplace_back(ModifiedList{target_cell->particles()});
  }

  src.clear();
}

namespace {
/**
 * @brief Fold coordinates to box and reset the old position.
 */
void fold_and_reset(Particle &p, BoxGeometry const &box_geo) {
  fold_position(p.pos(), p.image_box(), box_geo);

  p.pos_at_last_verlet_update() = p.pos();
}
} // namespace

void RegularSingleNodeDecompositoin::resort(bool global,
                                            std::vector<ParticleChange> &diff) {
  ParticleList displaced_parts;

  for (auto &c : local_cells()) {
    for (auto it = c->particles().begin(); it != c->particles().end();) {
      fold_and_reset(*it, m_box);

      auto target_cell = particle_to_cell(*it);

      /* Particle is in place */
      if (target_cell == c) {
        std::advance(it, 1);
        continue;
      }

      auto p = std::move(*it);
      it = c->particles().erase(it);
      diff.emplace_back(ModifiedList{c->particles()});

      /* Particle is not local */
      if (target_cell == nullptr) {
        throw std::runtime_error("Target cell could not be determined. Tshis "
                                 "shoudl not happen. Bug.");
      }
      /* Particle belongs on this node but is in the wrong cell. */
      else if (target_cell != c) {
        target_cell->particles().insert(std::move(p));
        diff.emplace_back(ModifiedList{target_cell->particles()});
      }
    }
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

void RegularSingleNodeDecompositoin::mark_cells() {
  m_local_cells.clear();
  m_ghost_cells.clear();

  int cnt_c = 0;
  for (int o = 0; o < cell_grid[2]; o++)
    for (int n = 0; n < cell_grid[1]; n++)
      for (int m = 0; m < cell_grid[0]; m++) {
        m_local_cells.push_back(&cells.at(cnt_c++));
      }
}

// this can probabl go. TODO
void RegularSingleNodeDecompositoin::fill_comm_cell_lists(
    ParticleList **part_lists, Utils::Vector3i const &lc,
    Utils::Vector3i const &hc) {
  for (int o = lc[0]; o <= hc[0]; o++)
    for (int n = lc[1]; n <= hc[1]; n++)
      for (int m = lc[2]; m <= hc[2]; m++) {
        auto const i = Utils::get_linear_index(o, n, m, ghost_cell_grid);

        *part_lists++ = &(cells.at(i).particles());
      }
}
Utils::Vector3d RegularSingleNodeDecompositoin::max_cutoff() const {
  auto dir_max_range = [this](int i) {
    return std::min(0.5 * m_box.length()[i], m_local_box.length()[i]);
  };

  return {dir_max_range(0), dir_max_range(1), dir_max_range(2)};
}

Utils::Vector3d RegularSingleNodeDecompositoin::max_range() const {
  return cell_size;
}
int RegularSingleNodeDecompositoin::calc_processor_min_num_cells() const {
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

// TODO, I suspect, this can be simplified
void RegularSingleNodeDecompositoin::create_cell_grid(double range) {
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
    auto const &box_l = m_box.length();
    auto const volume = Utils::product(box_l);
    auto const scale =
        std::cbrt(RegularSingleNodeDecompositoin::max_num_cells / volume);

    for (int i = 0; i < 3; i++) {
      /* this is at least 1 */
      cell_grid[i] = static_cast<int>(std::ceil(box_l[i] * scale));
      cell_range[i] = box_l[i] / static_cast<double>(cell_grid[i]);

      if (cell_range[i] < range) {
        /* ok, too many cells for this direction, set to minimum */
        cell_grid[i] = static_cast<int>(std::floor(box_l[i] / range));
        if (cell_grid[i] < 1) {
          runtimeErrorMsg() << "interaction range " << range << " in direction "
                            << i << " is larger than the box size " << box_l[i];
          cell_grid[i] = 1;
        }
        cell_range[i] = box_l[i] / static_cast<double>(cell_grid[i]);
      }
    }

    /* It may be necessary to asymmetrically assign the scaling to the
       coordinates, which the above approach will not do.
       For a symmetric box, it gives a symmetric result. Here we correct that.
       */
    for (;;) {
      n_local_cells = Utils::product(cell_grid);

      /* done */
      if (n_local_cells <= RegularSingleNodeDecompositoin::max_num_cells)
        break;

      /* find coordinate with the smallest cell range */
      int min_ind = 0;
      double min_size = cell_range[0];

      for (int i = 1; i < 3; i++) {
        if (cell_grid[i] > 1 && cell_range[i] < min_size) {
          min_ind = i;
          min_size = cell_range[i];
        }
      }

      cell_grid[min_ind]--;
      cell_range[min_ind] = m_box.length()[min_ind] / cell_grid[min_ind];
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

  if (n_local_cells > RegularSingleNodeDecompositoin::max_num_cells) {
    runtimeErrorMsg() << "no suitable cell grid found";
  }

  /* now set all dependent variables */
  int new_cells = 1;
  for (int i = 0; i < 3; i++) {
    new_cells *= cell_grid[i];
    cell_size[i] = m_box.length()[i] / static_cast<double>(cell_grid[i]);
    inv_cell_size[i] = 1.0 / cell_size[i];
  }

  /* allocate cell array and cell pointer arrays */
  cells.clear();
  cells.resize(new_cells);
  m_local_cells.resize(n_local_cells);
  m_ghost_cells.clear();
}

template <class K, class Comparator> auto make_flat_set(Comparator &&comp) {
  return boost::container::flat_set<K, std::remove_reference_t<Comparator>>(
      std::forward<Comparator>(comp));
}

// this needs checking... TODO
void RegularSingleNodeDecompositoin::init_cell_interactions() {

  int n_cells = m_local_cells.size();
  /* Tanslate a node local index (relative to the origin of the local grid)
   * to a global index. */
  /* Linear index in the cell grid, respecting perioidic boudaires */
  auto folded_linear_index = [&](Utils::Vector3i const &index) {
    auto const folded_index = (index % cell_grid;

    return get_linear_index(folded_index, cell_grid);
  };

  auto const start = Vector3i{0,0,0});
  auto const end = cell_grid; // not including ight hand boundary in loop

  /* loop all local cells */
  for (int o = start[2]; o < end[2]; o++)
    for (int n = start[1]; n < end[1]; n++)
      for (int m = start[0]; m < end[0]; m++) {
        /* next-nearest neighbors in every direction */
        Utils::Vector3i lower_index = {m - 1, n - 1, o - 1};
        Utils::Vector3i upper_index = {m + 1, n + 1, o + 1};

        //        /* In the fully connected case, we consider all cells
        //         * in the direction as neighbors, not only the nearest ones.
        //         */
        //        for (int i = 0; i < 3; i++) {
        //          if (dd.fully_connected[i]) {
        //            // Fully connected is only needed at the box surface
        //            if (i==0 and (n!=start[1] or n!=end[1]-1) and (o!=start[2]
        //            or o!=end[2]-1)) continue; if (i==1 and (m!=start[0] or
        //            m!=end[0]-1) and (o!=start[2] or o!=end[2]-1)) continue;
        //            if (i==2 and (m!=start[0] or m!=end[0]-1) and (n!=start[1]
        //            or n!=end[1]-1)) continue; lower_index[i] = 0;
        //            upper_index[i] = global_size[i] - 1;
        //          }
        //        }

        // In non-periodic direction, cut-off at start and end. I.e. mo
        // neighbors
        for (int i = 0; i < 3; i++) {
          if (not box_geo.periodic(i)) {
            lower_index[i] = start[i];
            upper_index[i] = end[i];
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
              neighbors.insert(Utils::Vector3i{r, q, p});
            }

        /* Red-black partition by global index. */
        auto const ind1 = folded_linear_index({m, n, o});

        std::vector<Cell *> red_neighbors;
        std::vector<Cell *> black_neighbors;
        for (auto const &neighbor : neighbors) {
          auto const ind2 = folded_linear_index(neighbor);
          /* Exclude cell itself */
          if (ind1 == ind2)
            continue;

          auto cell = &cells.at(
              get_linear_index(local_index(neighbor), ghost_cell_grid));
          if (ind2 > ind1) {
            red_neighbors.push_back(cell);
          } else {
            black_neighbors.push_back(cell);
          }
        }

        cells[get_linear_index(index({m, n, o}), cell_grid)].m_neighbors =
            Neighbors<Cell *>(red_neighbors, black_neighbors);
      }
}

manamespace {
  /** Revert the order of a communicator: After calling this the
   *  communicator is working in reverted order with exchanged
   *  communication types GHOST_SEND <-> GHOST_RECV.
   */
  void revert_comm_order(GhostCommunicator & comm) {
    /* revert order */
    boost::reverse(comm.communications);

    /* exchange SEND/RECV */
    for (auto &c : comm.communications) {
      if (c.type == GHOST_SEND)
        c.type = GHOST_RECV;
      else if (c.type == GHOST_RECV)
        c.type = GHOST_SEND;
      else if (c.type == GHOST_LOCL) {
        boost::reverse(c.part_lists);
      }
    }
  }

  /** Of every two communication rounds, set the first receivers to prefetch and
   *  poststore
   */
  void assign_prefetches(GhostCommunicator & comm) {
    for (auto it = comm.communications.begin(); it != comm.communications.end();
         it += 2) {
      auto next = std::next(it);
      if (it->type == GHOST_RECV && next->type == GHOST_SEND) {
        it->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
        next->type |= GHOST_PREFETCH | GHOST_PSTSTORE;
      }
    }
  }

} // namespace

GhostCommunicator RegularSingleNodeDecompositoin::prepare_comm() {
  int dir, lr, i, cnt, n_comm_cells[3];
  Utils::Vector3i lc{}, hc{}, done{};

  auto const comm_info = Utils::Mpi::cart_get<3>(m_comm);
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(m_comm);

  /* calculate number of communications */
  std::size_t num = 0;
  for (dir = 0; dir < 3; dir++) {
    for (lr = 0; lr < 2; lr++) {
      /* No communication for border of non periodic direction */
      if (comm_info.dims[dir] == 1)
        num++;
      else
        num += 2;
    }
  }

  /* prepare communicator */
  auto ghost_comm = GhostCommunicator{m_comm, num};

  /* number of cells to communicate in a direction */
  n_comm_cells[0] = cell_grid[1] * cell_grid[2];
  n_comm_cells[1] = cell_grid[2] * ghost_cell_grid[0];
  n_comm_cells[2] = ghost_cell_grid[0] * ghost_cell_grid[1];

  cnt = 0;
  /* direction loop: x, y, z */
  for (dir = 0; dir < 3; dir++) {
    lc[(dir + 1) % 3] = 1 - done[(dir + 1) % 3];
    lc[(dir + 2) % 3] = 1 - done[(dir + 2) % 3];
    hc[(dir + 1) % 3] = cell_grid[(dir + 1) % 3] + done[(dir + 1) % 3];
    hc[(dir + 2) % 3] = cell_grid[(dir + 2) % 3] + done[(dir + 2) % 3];
    /* lr loop: left right */
    /* here we could in principle build in a one sided ghost
       communication, simply by taking the lr loop only over one
       value */
    for (lr = 0; lr < 2; lr++) {
      if (comm_info.dims[dir] == 1) {
        /* just copy cells on a single node */
        ghost_comm.communications[cnt].type = GHOST_LOCL;
        ghost_comm.communications[cnt].node = m_comm.rank();

        /* Buffer has to contain Send and Recv cells -> factor 2 */
        ghost_comm.communications[cnt].part_lists.resize(2 * n_comm_cells[dir]);

        /* fill send ghost_comm cells */
        lc[dir] = hc[dir] = 1 + lr * (cell_grid[dir] - 1);

        fill_comm_cell_lists(ghost_comm.communications[cnt].part_lists.data(),
                             lc, hc);

        /* fill recv ghost_comm cells */
        lc[dir] = hc[dir] = 0 + (1 - lr) * (cell_grid[dir] + 1);

        /* place receive cells after send cells */
        fill_comm_cell_lists(
            &ghost_comm.communications[cnt].part_lists[n_comm_cells[dir]], lc,
            hc);

        cnt++;
      } else {
        /* i: send/recv loop */
        for (i = 0; i < 2; i++) {
          if ((comm_info.coords[dir] + i) % 2 == 0) {
            ghost_comm.communications[cnt].type = GHOST_SEND;
            ghost_comm.communications[cnt].node = node_neighbors[2 * dir + lr];
            ghost_comm.communications[cnt].part_lists.resize(n_comm_cells[dir]);
            /* prepare folding of ghost positions */

            lc[dir] = hc[dir] = 1 + lr * (cell_grid[dir] - 1);

            fill_comm_cell_lists(
                ghost_comm.communications[cnt].part_lists.data(), lc, hc);
            cnt++;
          }
          if ((comm_info.coords[dir] + (1 - i)) % 2 == 0) {
            ghost_comm.communications[cnt].type = GHOST_RECV;
            ghost_comm.communications[cnt].node =
                node_neighbors[2 * dir + (1 - lr)];
            ghost_comm.communications[cnt].part_lists.resize(n_comm_cells[dir]);

            lc[dir] = hc[dir] = (1 - lr) * (cell_grid[dir] + 1);

            fill_comm_cell_lists(
                ghost_comm.communications[cnt].part_lists.data(), lc, hc);
            cnt++;
          }
        }
      }
      done[dir] = 1;
    }
  }

  return ghost_comm;
}

RegularSingleNodeDecompositoin::RegularSingleNodeDecompositoin(
    double range, BoxGeometry const &box_geo)
    : m_box(box_geo) {
  /* set up new regular decomposition cell structure */
  create_cell_grid(range);

  /* setup cell neighbors */
  init_cell_interactions();

  /* mark local and ghost cells */
  mark_cells();

  /* create communicators */
  m_exchange_ghosts_comm = {];
  m_collect_ghost_force_comm = {};
}

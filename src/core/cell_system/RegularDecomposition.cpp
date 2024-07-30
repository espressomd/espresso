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

#include "cell_system/RegularDecomposition.hpp"

#include "cell_system/Cell.hpp"

#include "communication.hpp"
#include "error_handling/RuntimeErrorStream.hpp"
#include "errorhandling.hpp"
#include "system/System.hpp"

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
        diff.emplace_back(RemovedParticle{p.id()});
        displaced_parts.insert(std::move(p));
      }
      /* Particle belongs on this node but is in the wrong cell. */
      else if (target_cell != c) {
        target_cell->particles().insert(std::move(p));
        diff.emplace_back(ModifiedList{target_cell->particles()});
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

void RegularDecomposition::fill_comm_cell_lists(ParticleList **part_lists,
                                                Utils::Vector3i const &lc,
                                                Utils::Vector3i const &hc) {
  for (int o = lc[0]; o <= hc[0]; o++)
    for (int n = lc[1]; n <= hc[1]; n++)
      for (int m = lc[2]; m <= hc[2]; m++) {
        auto const i = Utils::get_linear_index(o, n, m, ghost_cell_grid);

        *part_lists++ = &(cells.at(i).particles());
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

  auto const halo = Utils::Vector3i{1, 1, 1};
  auto const cart_info = Utils::Mpi::cart_get<3>(m_comm);
  auto const &node_pos = cart_info.coords;
  auto const &node_grid = ::communicator.node_grid;
  auto const global_halo_offset = hadamard_product(node_pos, cell_grid) - halo;
  auto const global_size = hadamard_product(node_grid, cell_grid);
  auto const at_boundary = [&global_size](int coord, Utils::Vector3i cell_idx) {
    return (cell_idx[coord] == 0 or cell_idx[coord] == global_size[coord]);
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
  }

  /* We only consider local cells (e.g. not halo cells), which
   * span the range [(1,1,1), cell_grid) in local coordinates. */
  auto const start = global_index(Utils::Vector3i{1, 1, 1});
  auto const end = start + cell_grid;

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
            upper_index[fc_direction] = global_size[fc_boundary];
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

        cells[get_linear_index(local_index({m, n, o}), ghost_cell_grid)]
            .m_neighbors = Neighbors<Cell *>(red_neighbors, black_neighbors);
      }
}

namespace {
/** Revert the order of a communicator: After calling this the
 *  communicator is working in reverted order with exchanged
 *  communication types GHOST_SEND <-> GHOST_RECV.
 */
void revert_comm_order(GhostCommunicator &comm) {
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
void assign_prefetches(GhostCommunicator &comm) {
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

GhostCommunicator RegularDecomposition::prepare_comm() {
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

  /* create communicators */
  m_exchange_ghosts_comm = prepare_comm();
  m_collect_ghost_force_comm = prepare_comm();

  /* collect forces has to be done in reverted order! */
  revert_comm_order(m_collect_ghost_force_comm);

  assign_prefetches(m_exchange_ghosts_comm);
  assign_prefetches(m_collect_ghost_force_comm);
}

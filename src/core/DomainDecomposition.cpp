#include "DomainDecomposition.hpp"

#include "errorhandling.hpp"

#include <boost/mpi/collectives.hpp>
#include <utils/mpi/sendrecv.hpp>

/** Returns pointer to the cell which corresponds to the position if the
 *  position is in the nodes spatial domain otherwise a nullptr pointer.
 */
Cell *DomainDecomposition::position_to_cell(const Utils::Vector3d &pos) {
  Utils::Vector3i cpos;

  for (int i = 0; i < 3; i++) {
    cpos[i] = static_cast<int>(std::floor(pos[i] * inv_cell_size[i])) + 1 -
              cell_offset[i];

    /* particles outside our box. Still take them if
       nonperiodic boundary. We also accept the particle if we are at
       the box boundary, and the particle is within the box. In this case
       the particle belongs here and could otherwise potentially be dismissed
       due to rounding errors. */
    if (cpos[i] < 1) {
      if ((!box_geo.periodic(i) or (pos[i] >= box_geo.length()[i])) &&
          local_geo.boundary()[2 * i])
        cpos[i] = 1;
      else
        return nullptr;
    } else if (cpos[i] > cell_grid[i]) {
      if ((!box_geo.periodic(i) or (pos[i] < box_geo.length()[i])) &&
          local_geo.boundary()[2 * i + 1])
        cpos[i] = cell_grid[i];
      else
        return nullptr;
    }
  }

  auto const ind = get_linear_index(cpos, ghost_cell_grid);
  return &(cells.at(ind));
}

void DomainDecomposition::move_if_local(ParticleList &src, ParticleList &rest,
                                        std::vector<Cell *> &modified_cells) {
  for (auto &part : src) {
    auto target_cell = position_to_cell(part.r.p);

    if (target_cell) {
      target_cell->particles().insert(std::move(part));
      modified_cells.push_back(target_cell);
    } else {
      rest.insert(std::move(part));
    }
  }

  src.clear();
}

/**
 * @brief Split particle list by direction.
 *
 * Moves all particles from src into left
 * and right depending if they belong to
 * the left or right side from local node
 * in direction dir.
 *
 * @param src Particles to sort.
 * @param left Particles that should go to the left
 * @param right Particles that should go to the right
 * @param dir Direction to consider.
 */
void DomainDecomposition::move_left_or_right(ParticleList &src,
                                             ParticleList &left,
                                             ParticleList &right,
                                             int dir) const {
  for (auto it = src.begin(); it != src.end();) {
    if ((get_mi_coord(it->r.p[dir], local_geo.my_left()[dir],
                      box_geo.length()[dir], box_geo.periodic(dir)) < 0.0) and
        (box_geo.periodic(dir) || (local_geo.boundary()[2 * dir] == 0))) {
      left.insert(std::move(*it));
      it = src.erase(it);
    } else if ((get_mi_coord(it->r.p[dir], local_geo.my_right()[dir],
                             box_geo.length()[dir],
                             box_geo.periodic(dir)) >= 0.0) and
               (box_geo.periodic(dir) ||
                (local_geo.boundary()[2 * dir + 1] == 0))) {
      right.insert(std::move(*it));
      it = src.erase(it);
    } else {
      ++it;
    }
  }
}

void DomainDecomposition::exchange_neighbors(
    ParticleList &pl, std::vector<Cell *> &modified_cells) {
  auto const node_neighbors = Utils::Mpi::cart_neighbors<3>(comm);
  static ParticleList send_buf_l, send_buf_r, recv_buf_l, recv_buf_r;

  for (int dir = 0; dir < 3; dir++) {
    /* Single node direction, no action needed. */
    if (Utils::Mpi::cart_get<3>(comm).dims[dir] == 1) {
      continue;
      /* In this (common) case left and right neighbors are
         the same, and we need only one communication */
    }
    if (Utils::Mpi::cart_get<3>(comm).dims[dir] == 2) {
      move_left_or_right(pl, send_buf_l, send_buf_l, dir);

      Utils::Mpi::sendrecv(comm, node_neighbors[2 * dir], 0, send_buf_l,
                           node_neighbors[2 * dir], 0, recv_buf_l);

      send_buf_l.clear();
    } else {
      using boost::mpi::request;
      using Utils::Mpi::isendrecv;

      move_left_or_right(pl, send_buf_l, send_buf_r, dir);

      auto req_l = isendrecv(comm, node_neighbors[2 * dir], 0, send_buf_l,
                             node_neighbors[2 * dir], 0, recv_buf_l);
      auto req_r = isendrecv(comm, node_neighbors[2 * dir + 1], 0, send_buf_r,
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

void DomainDecomposition::resort(bool global, ParticleList &pl,
                                 std::vector<Cell *> &modified_cells) {
  if (global) {
    auto const grid = Utils::Mpi::cart_get<3>(comm).dims;
    /* Worst case we need grid - 1 rounds per direction.
     * This correctly implies that if there is only one node,
     * no action should be taken. */
    int rounds_left = grid[0] + grid[1] + grid[2] - 3;
    for (; rounds_left > 0; rounds_left--) {
      exchange_neighbors(pl, modified_cells);

      auto left_over =
          boost::mpi::all_reduce(comm, pl.size(), std::plus<size_t>());

      if (left_over == 0) {
        break;
      }
    }
  } else {
    exchange_neighbors(pl, modified_cells);
  }
}
void DomainDecomposition::mark_cells() {
  int cnt_c = 0;

  m_local_cells.clear();
  m_ghost_cells.clear();

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
void DomainDecomposition::fill_comm_cell_lists(ParticleList **part_lists,
                                               const Utils::Vector3i &lc,
                                               const Utils::Vector3i &hc) {
  for (int o = lc[0]; o <= hc[0]; o++)
    for (int n = lc[1]; n <= hc[1]; n++)
      for (int m = lc[2]; m <= hc[2]; m++) {
        auto const i = Utils::get_linear_index(o, n, m, ghost_cell_grid);

        *part_lists++ = &(cells.at(i).particles());
      }
}
Utils::Vector3d DomainDecomposition::max_range() const {
  auto dir_max_range = [this](int i) {
    if (fully_connected[i]) {
      return std::numeric_limits<double>::infinity();
    }

    return std::min(0.5 * box_geo.length()[i], local_geo.length()[i]);
  };

  return {dir_max_range(0), dir_max_range(1), dir_max_range(2)};
}
int DomainDecomposition::calc_processor_min_num_cells() const {
  /* the minimal number of cells can be lower if there are at least two nodes
     serving a direction,
     since this also ensures that the cell size is at most half the box
     length. However, if there is only one processor for a direction, there
     have to be at least two cells for this direction. */
  return boost::accumulate(Utils::Mpi::cart_get<3>(comm).dims, 1,
                           [](int n_cells, int grid) {
                             return (grid == 1) ? 2 * n_cells : n_cells;
                           });
}

void DomainDecomposition::create_cell_grid(double range) {
  auto const cart_info = Utils::Mpi::cart_get<3>(comm);

  int i, n_local_cells, new_cells;
  double cell_range[3];

  /* initialize */
  cell_range[0] = cell_range[1] = cell_range[2] = range;

  /* Min num cells can not be smaller than calc_processor_min_num_cells. */
  int min_num_cells = calc_processor_min_num_cells();

  if (range <= 0.) {
    /* this is the non-interacting case */
    const int cells_per_dir = std::ceil(std::pow(min_num_cells, 1. / 3.));

    cell_grid[0] = cells_per_dir;
    cell_grid[1] = cells_per_dir;
    cell_grid[2] = cells_per_dir;

    n_local_cells = cell_grid[0] * cell_grid[1] * cell_grid[2];
  } else {
    /* Calculate initial cell grid */
    double volume = local_geo.length()[0];
    for (i = 1; i < 3; i++)
      volume *= local_geo.length()[i];
    double scale = pow(DomainDecomposition::max_num_cells / volume, 1. / 3.);
    for (i = 0; i < 3; i++) {
      /* this is at least 1 */
      cell_grid[i] = (int)ceil(local_geo.length()[i] * scale);
      cell_range[i] = local_geo.length()[i] / cell_grid[i];

      if (cell_range[i] < range) {
        /* ok, too many cells for this direction, set to minimum */
        cell_grid[i] = (int)floor(local_geo.length()[i] / range);
        if (cell_grid[i] < 1) {
          runtimeErrorMsg()
              << "interaction range " << range << " in direction " << i
              << " is larger than the local box size " << local_geo.length()[i];
          cell_grid[i] = 1;
        }
        cell_range[i] = local_geo.length()[i] / cell_grid[i];
      }
    }

    /* It may be necessary to asymmetrically assign the scaling to the
       coordinates, which the above approach will not do.
       For a symmetric box, it gives a symmetric result. Here we correct that.
       */
    for (;;) {
      n_local_cells = cell_grid[0] * cell_grid[1] * cell_grid[2];

      /* done */
      if (n_local_cells <= DomainDecomposition::max_num_cells)
        break;

      /* find coordinate with the smallest cell range */
      int min_ind = 0;
      double min_size = cell_range[0];

      for (i = 1; i < 3; i++) {
        if (cell_grid[i] > 1 && cell_range[i] < min_size) {
          min_ind = i;
          min_size = cell_range[i];
        }
      }

      cell_grid[min_ind]--;
      cell_range[min_ind] = local_geo.length()[min_ind] / cell_grid[min_ind];
    }

    /* sanity check */
    if (n_local_cells < min_num_cells) {
      runtimeErrorMsg()
          << "number of cells " << n_local_cells << " is smaller than minimum "
          << min_num_cells
          << " (interaction range too large or min_num_cells too large)";
    }
  }

  /* quit program if unsuccessful */
  if (n_local_cells > DomainDecomposition::max_num_cells) {
    runtimeErrorMsg() << "no suitable cell grid found ";
  }

  auto const node_pos = cart_info.coords;

  /* now set all dependent variables */
  new_cells = 1;
  for (i = 0; i < 3; i++) {
    ghost_cell_grid[i] = cell_grid[i] + 2;
    new_cells *= ghost_cell_grid[i];
    cell_size[i] = local_geo.length()[i] / (double)cell_grid[i];
    inv_cell_size[i] = 1.0 / cell_size[i];
    cell_offset[i] = node_pos[i] * cell_grid[i];
  }

  /* allocate cell array and cell pointer arrays */
  cells.clear();
  cells.resize(new_cells);
  m_local_cells.resize(n_local_cells);
  m_ghost_cells.resize(new_cells - n_local_cells);
}

void DomainDecomposition::init_cell_interactions() {
  /* loop all local cells */
  for (int o = 1; o < cell_grid[2] + 1; o++)
    for (int n = 1; n < cell_grid[1] + 1; n++)
      for (int m = 1; m < cell_grid[0] + 1; m++) {

        auto const ind1 = get_linear_index(m, n, o, ghost_cell_grid);

        std::vector<Cell *> red_neighbors;
        std::vector<Cell *> black_neighbors;

        /* loop all neighbor cells */
        int lower_index[3] = {m - 1, n - 1, o - 1};
        int upper_index[3] = {m + 1, n + 1, o + 1};

        for (int i = 0; i < 3; i++) {
          if (fully_connected[i]) {
            lower_index[i] = 0;
            upper_index[i] = ghost_cell_grid[i] - 1;
          }
        }

        for (int p = lower_index[2]; p <= upper_index[2]; p++)
          for (int q = lower_index[1]; q <= upper_index[1]; q++)
            for (int r = lower_index[0]; r <= upper_index[0]; r++) {
              auto const ind2 = get_linear_index(r, q, p, ghost_cell_grid);
              if (ind2 > ind1) {
                red_neighbors.push_back(&cells.at(ind2));
              } else {
                black_neighbors.push_back(&cells.at(ind2));
              }
            }
        cells.at(ind1).m_neighbors =
            Neighbors<Cell *>(red_neighbors, black_neighbors);
      }
}

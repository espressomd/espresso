/*
 * Copyright (C) 2010-2020 The ESPResSo project
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
#include "generic_dd.hpp"

#ifdef HAVE_REPA
#include <boost/mpi/collectives.hpp>
#include <boost/mpi/nonblocking.hpp>
#include <boost/range/adaptor/uniqued.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <memory>
#include <repa/repa.hpp>

#include "communication.hpp" // comm_cart
#include "errorhandling.hpp"
#include "grid.hpp" // fold_position()
#include <utils/Vector.hpp>
#include <utils/mpi/waitany.hpp>

namespace {
/** Wraps the underlying grid and type.
 */
struct GridImplementationWrapper {
  std::unique_ptr<repa::grids::ParallelLCGrid> pargrid;
  /** Gridtype that is supposed to be created the next time a new pargrid
   * is instantiated.
   */
  repa::GridType new_gridtype;

  GridImplementationWrapper()
      : pargrid(nullptr), new_gridtype(repa::GridType::NONE) {}

  repa::grids::ParallelLCGrid *operator->() { return pargrid.get(); }
  const repa::grids::ParallelLCGrid *operator->() const {
    return pargrid.get();
  }
};

/** Returns the one GridImplementationWrapper singleton.
 */
GridImplementationWrapper &grid_instance() {
  static GridImplementationWrapper g;
  return g;
}
} // namespace

namespace generic_dd {

namespace impl {

/** Convert an ESPResSo vec to a repa vec.
 */
static inline repa::Vec3d to_repa_vec(const Utils::Vector3d &v) {
  return {v[0], v[1], v[2]};
}

/** Convert a repa vec to an ESPResSo vec
 */
static inline Utils::Vector3d to_espresso_vec(const repa::Vec3d &v) {
  return {v[0], v[1], v[2]};
}

/** Folds back a position into the primary simulation box.
 */
static Utils::Vector3d folded_pos(const Utils::Vector3d &pos) {
  auto pfold = pos;
  Utils::Vector3i im = {0, 0, 0};
  fold_position(pfold, im, box_geo);
  return pfold;
}

/** Returns the rank of the process which is responsible for the cell
 * "pos" falls into.
 */
static int position_to_node(const Utils::Vector3d &pos) {
  return grid_instance()->position_to_rank(to_repa_vec(folded_pos(pos)));
}

/** Folded version of Pargrid::position_to_neighidx.
 */
static int position_to_neighidx(const Utils::Vector3d &pos) {
  return grid_instance()->position_to_neighidx(to_repa_vec(folded_pos(pos)));
}

/** Maps a position "pos" to a cell, local or ghost.
 */
static Cell *position_to_cell(const Utils::Vector3d &pos) {
  try {
    const auto i =
        grid_instance()->position_to_cell_index(to_repa_vec(folded_pos(pos)));
    return cell_structure.m_local_cells[i];
  } catch (const std::domain_error &e) {
    return nullptr;
  }
}

/** Takes the ownership of a ParticleList and moves the particles to their
 * corresponding local cells. If the particle does not belong on this node,
 * aborts.
 * Also keeps track of the cells the function appends particles to and
 * stores pointers to these in "modified_cells".
 */
static void insert_particles(ParticleList &&recvbuf,
                             std::vector<Cell *> &modified_cells) {
  for (auto &p : recvbuf) {
    // fold_position(p.r.p, p.l.i, box_geo);
    Cell *c = position_to_cell(p.r.p);

    if (!c)
      runtimeErrorMsg()
          << "[" << this_node
          << "] Insertion: Particle does not belong on this node but on: "
          << position_to_node(p.r.p);

    c->particles().insert(std::move(p));
    modified_cells.push_back(c);
  }
  recvbuf.clear();
}

/** Fills the particle sendbuffers "sendbuf" for global or local resorting with
 * the particles from "pl". The ownership of all particles in "pl" is
 * transferred and pl is emptied.
 */
static void fill_sendbufs(bool global, ParticleList *pl,
                          std::vector<ParticleList> &sendbuf) {
  for (auto &part : *pl) {
    const int neigh =
        (global ? position_to_node : position_to_neighidx)(part.r.p);
    sendbuf[neigh].insert(std::move(part));
  }
  pl->clear();
}

namespace exchg {
static bool in_progress = false;
static int nneigh;
static std::vector<ParticleList> sendbuf, recvbuf;
static std::vector<boost::mpi::request> rreq, sreq;
} // namespace exchg

/** Starts the transfer of all particles from "pl" to their corresponding
 * owners.
 * If global == true, communicates globally with with all processes to
 * exchange particles from "pl".
 * If global == false, communicates only with the neighboring nodes. If
 * a particle is encountered that cannot be resolved to any neighboring node,
 * this procedure will abort the program.
 */
static void exchange_particles_begin(bool global, ParticleList *pl) {
  assert(!exchg::in_progress);

  // Some number that can be recognized if they appear in an error message
  static constexpr int GDD_EX_TAG = 11011;
  // On a global exchange, communicate with 1..n. On local exchange, only with
  // the neighbors.
  auto neighbor_rank = [global](int neigh) {
    return global ? neigh : grid_instance()->neighbor_rank(neigh);
  };

  exchg::nneigh = global ? comm_cart.size() : grid_instance()->n_neighbors();

  // Send and receive data
  exchg::sendbuf.resize(exchg::nneigh);
  exchg::recvbuf.resize(exchg::nneigh);

  exchg::rreq.reserve(exchg::nneigh);
  for (int i = 0; i < exchg::nneigh; i++) {
    exchg::rreq.push_back(
        comm_cart.irecv(neighbor_rank(i), GDD_EX_TAG, exchg::recvbuf[i]));
  }

  fill_sendbufs(global, pl, exchg::sendbuf);
  assert(pl->size() == 0);

  // Send particle data
  exchg::sreq.reserve(exchg::nneigh);
  for (int i = 0; i < exchg::nneigh; i++) {
    exchg::sreq.push_back(
        comm_cart.isend(neighbor_rank(i), GDD_EX_TAG, exchg::sendbuf[i]));
  }

  exchg::in_progress = true;
}

/** Waits for the end of all communications still in progress and
 * keeps track of the cells that particles were newly added to.
 */
static void exchange_particles_end(std::vector<Cell *> &modified_cells) {
  assert(exchg::in_progress);

  // Wait for receives and immediately process them
  for (int i = 0; i < exchg::nneigh; i++) {
    const auto p = Utils::Mpi::wait_any<Utils::Mpi::Status::Ignore>(
                       std::begin(exchg::rreq), std::end(exchg::rreq))
                       .iterator;
    const auto idx = std::distance(std::begin(exchg::rreq), p);
    insert_particles(std::move(exchg::recvbuf[idx]), modified_cells);
  }

  boost::mpi::wait_all(std::begin(exchg::sreq), std::end(exchg::sreq));
  for (auto &pl : exchg::sendbuf) {
    pl.clear();
  }

  exchg::in_progress = false;
  exchg::sendbuf.clear();
  exchg::recvbuf.clear();
  exchg::rreq.clear();
  exchg::sreq.clear();
}

/** Fills the neighbor lists of all local cells.
 */
static void fill_cell_inter_lists() {
  for (int i = 0; i < cell_structure.m_local_cells.size(); i++) {
    std::vector<Cell *> red_neighbors;
    std::vector<Cell *> black_neighbors;

    for (int n = 0; n < 27; ++n) {
      (n >= 1 && n < 14 ? red_neighbors : black_neighbors)
          .push_back(&cells[grid_instance()->cell_neighbor_index(i, n)]);
    }
    cell_structure.m_local_cells[i]->m_neighbors =
        Neighbors<Cell *>(red_neighbors, black_neighbors);
  }
}

/** Fills a GhostCommunicator object with the necessary information
 * for boundary/ghost layer exchanges.
 */
static void fill_communicator(GhostCommunicator *comm) {
  const auto gexds = grid_instance()->get_boundary_info();
  const size_t ncomm = gexds.size();
  prepare_comm(comm, static_cast<int>(2 * ncomm));
  comm->async = true;
  for (size_t i = 0; i < ncomm; ++i) {
    const auto &gexd = gexds[i];
    comm->comm[i].type = GHOST_SEND;
    comm->comm[i].node = gexd.dest;
    comm->comm[i].part_lists.resize(gexd.send.size());
    for (int j = 0; j < gexd.send.size(); ++j)
      comm->comm[i].part_lists[j] = &cells[gexd.send[j]];

    comm->comm[ncomm + i].type = GHOST_RECV;
    comm->comm[ncomm + i].node = gexd.dest;
    comm->comm[ncomm + i].part_lists.resize(gexd.recv.size());
    for (int j = 0; j < gexd.recv.size(); ++j)
      comm->comm[ncomm + i].part_lists[j] = &cells[gexd.recv[j]];
  }
}

/** Turns all send GhostCommunications into receive ones and vice versa.
 */
static void revert_communicator(GhostCommunicator *comm) {
  std::for_each(comm->comm.begin(), comm->comm.end(),
                [](GhostCommunication &gc) {
                  if (gc.type == GHOST_SEND)
                    gc.type = GHOST_RECV;
                  else if (gc.type == GHOST_RECV)
                    gc.type = GHOST_SEND;
                });
}

/** Fills "local_cells" and "ghost_cells" with the respective pointers to
 * "cells".
 */
static void mark_cells() {
  int c = 0;
  for (int i = 0; i < grid_instance()->n_local_cells(); i++)
    cell_structure.m_local_cells[i] = &cells[c++];
  for (int i = 0; i < grid_instance()->n_ghost_cells(); i++)
    cell_structure.m_ghost_cells[i] = &cells[c++];
}

/** Returns a ParticleList containing all out of box particles
 * from an range of Cell pointers.
 * Out of box particles are moved from.
 */
template <typename It> ParticleList filter_oob_particles(It first, It last) {
  ParticleList oobparts;
  auto extract_oob = [&oobparts](ParticleList &src) {
    for (auto it = src.begin(); it != src.end();) {
      if (position_to_node(folded_pos(it->r.p)) != this_node) {
        oobparts.insert(std::move(*it));
        it = src.erase(it);
      } else {
        ++it;
      }
    }
  };
  for (; first != last; ++first) {
    extract_oob((*first)->particles());
  }
  return oobparts;
}

/** Returns the center of mass of this subdomain.
 */
static repa::Vec3d center_of_mass() {
  int npart = 0;
  repa::Vec3d c = {{0., 0., 0.}};

  for (const auto &p : Cells::particles(cell_structure.local_cells())) {
    npart++;
    for (int d = 0; d < 3; ++d)
      c[d] += p.r.p[d];
  }

  for (int d = 0; d < 3; ++d)
    c[d] /= npart;

  return c;
}

/** Updates the particle index for cells modified by fn.
 * @param fn Function void(std::vector<Cell*>&) which stores modified cells in
 * parameter handed to it.
 */
template <typename Func> void call_with_modified_cells_update(Func &&fn) {
  static std::vector<Cell *> modified_cells;
  modified_cells.clear();
  fn(modified_cells);
  boost::sort(modified_cells);
  for (auto cell : modified_cells | boost::adaptors::uniqued) {
    cell_structure.update_particle_index(cell->particles());
  }
}

/** Ensures that each process has at least one cell.
 * If this cannot be ensured, aborts the program.
 */
static void ensure_at_least_one_cell() {
  // An alternative to the following is to simply abort if
  // grid_instance()->n_local_cells() == 0. Note as below, do not use
  // runtimeErrorMsg to abort in this case.
  std::vector<int> ncells_per_proc;
  boost::mpi::all_gather(comm_cart, grid_instance()->n_local_cells(),
                         ncells_per_proc);
  const auto is_zero = [](int i) { return i == 0; };

  if (std::any_of(ncells_per_proc.begin(), ncells_per_proc.end(), is_zero)) {
    // Try repartitioning, all cells equally weighted.
    if (grid_instance()->n_local_cells() == 0)
      std::cerr << "WARN: A process has 0 local cells. "
                << "Trying to repartition as remedy." << std::endl;
    grid_instance()->repartition(
        []() {
          return std::vector<double>(grid_instance()->n_local_cells(), 1.0);
        },
        [](int i, int j) { return 1.0; }, []() {});

    if (grid_instance()->n_local_cells() == 0) {
      // Don't use runtimeErrorMsg. This process will segfault on its way to the
      // next runtime error message collection.
      std::cerr << "ERROR: A process still has 0 local cells. "
                << "System too small for this number of processes."
                << std::endl;
      MPI_Abort(comm_cart, 1);
    }
  }
}

} // namespace impl

void on_geometry_change(int flags, double range) {
  auto grid_size = grid_instance()->grid_size();
  // Don't change cell system if it is still sufficient
  if (flags & CELL_FLAG_GRIDCHANGED ||
      grid_size[0] != std::max<int>(box_geo.length()[0] / range, 1) ||
      grid_size[1] != std::max<int>(box_geo.length()[1] / range, 1) ||
      grid_size[2] != std::max<int>(box_geo.length()[2] / range, 1)) {
    cells_re_init(CELL_STRUCTURE_CURRENT, range);
  }
}

void topology_init(double range, bool is_repart) {
  if (!box_geo.periodic(0) || !box_geo.periodic(1) || !box_geo.periodic(2)) {
    runtimeErrorMsg() << "Generic_dd does only support fully periodic domains.";
  }

  cell_structure.type = CELL_STRUCTURE_GENERIC_DD;
  cell_structure.particle_to_cell = [](const Particle &p) {
    return impl::position_to_cell(p.r.p);
  };

  // Create new grid only if not repartitioning.
  // This lets us effectively call cells_re_init after lb.
  if (!is_repart) {
    auto ep = repa::ExtraParams{};
    ep.subdomain_midpoint = impl::center_of_mass;
    grid_instance().pargrid =
        repa::make_pargrid(grid_instance().new_gridtype, comm_cart,
                           impl::to_repa_vec(box_geo.length()), range, ep);
  }

  impl::ensure_at_least_one_cell();

  cells.resize(grid_instance()->n_local_cells() +
               grid_instance()->n_ghost_cells());
  cell_structure.m_local_cells.resize(grid_instance()->n_local_cells());
  cell_structure.m_ghost_cells.resize(grid_instance()->n_ghost_cells());
  cell_structure.max_range =
      impl::to_espresso_vec(grid_instance()->cell_size());

  impl::mark_cells();
  impl::fill_cell_inter_lists();

  impl::fill_communicator(&cell_structure.exchange_ghosts_comm);
  impl::fill_communicator(&cell_structure.collect_ghost_force_comm);
  impl::revert_communicator(&cell_structure.collect_ghost_force_comm);

  if (is_repart) {
    impl::call_with_modified_cells_update(impl::exchange_particles_end);
  }
}

void topology_release() {
  free_comm(&cell_structure.exchange_ghosts_comm);
  free_comm(&cell_structure.collect_ghost_force_comm);

  // Free ghost cell storage
  cell_structure.m_ghost_cells.resize(0);
}

void exchange_particles(int global_flag, ParticleList *displaced_particles,
                        std::vector<Cell *> &modified_cells) {
  impl::exchange_particles_begin(global_flag, displaced_particles);
  impl::exchange_particles_end(modified_cells);
}

void set_grid(const std::string &desc) {
  grid_instance().new_gridtype = repa::parse_grid_type(desc);
}

void repartition(const repart::Metric &m) {
  auto exchange_start_callback = []() {
    auto oobparts =
        impl::filter_oob_particles(cell_structure.m_local_cells.begin(),
                                   cell_structure.m_local_cells.end());
    impl::exchange_particles_begin(true, &oobparts);
  };

  using namespace std::placeholders;
  auto ccm = std::bind(&repart::Metric::cell_cell_weight, m, _1, _2);

  if (grid_instance()->repartition(m, ccm, exchange_start_callback))
    cells_re_init(CELL_STRUCTURE_CURRENT, cell_structure.min_range, true);
}

void command(const std::string &cmd) { grid_instance()->command(cmd); }

std::vector<std::string> librepa_supported_grid_types() {
  auto s = repa::supported_grid_types();
  auto v = std::vector<std::string>{};
  std::transform(std::begin(s), std::end(s), std::back_inserter(v),
                 repa::grid_type_to_string);
  return v;
}

} // namespace generic_dd

#else // HAVE_REPA

#include "communication.hpp" // comm_cart
#include "errorhandling.hpp"

namespace generic_dd {

std::vector<std::string> librepa_supported_grid_types() { return {}; }

/** Prints an error message that indicates missing features and aborts the
 * program execution.
 *
 * Note, we cannot use ESPResSo-internal error handline (runtimeErrorMsg, etc.)
 * Using it would make this function return. If it did,
 * it would leave ESPResSo in an undefined state (e.g.
 * no local cells are created).
 */
[[noreturn]] static void err_not_compiled_in() {
  // 1 error message is enough
  if (comm_cart.rank() == 0) {
    std::cerr << "Generic_dd is only available if ESPResSo has been "
                 "compiled with 'repa'."
              << std::endl;
  }
  comm_cart
      .barrier(); // No guarantee at all but increases the hope that no process
                  // will call abort before rank 0 has flushed its stderr buffer
  errexit();
  // Not reached.
}

/** Defines a function "fn_name" returning "retval" and following arguments
 * to call "err_not_compiled_in".
 */
#define errdef(retval, fn_name, ...)                                           \
  retval fn_name(__VA_ARGS__) { err_not_compiled_in(); }

/** Error defines for all exported functionality.
 * This helps eliminating #ifdefs at the inclusion and caller site.
 */
// clang-format off
errdef(void, on_geometry_change, int, double)
errdef(void, topology_init, double, bool)
errdef(void, topology_release, void)
errdef(void, exchange_particles, int, ParticleList *, std::vector<Cell *>&)
errdef(void, set_grid, const std::string &)
errdef(void, repartition, const repart::Metric &)
errdef(void, command, const std::string &) // clang-format on

} // namespace generic_dd

#endif

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

#include "particle_node.hpp"

#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "partCfg_global.hpp"
#include "particle_data.hpp"

#include <utils/Cache.hpp>
#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/keys.hpp>
#include <utils/mpi/gatherv.hpp>

#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/collectives/reduce.hpp>
#include <boost/mpi/collectives/scatter.hpp>
#include <boost/optional.hpp>
#include <boost/range/algorithm/sort.hpp>
#include <boost/range/numeric.hpp>

#include <algorithm>
#include <cmath>
#include <functional>
#include <iterator>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

constexpr auto some_tag = 42;

/** @brief Enable particle type tracking in @ref particle_type_map. */
static bool type_list_enable;

/** @brief Mapping particle types to lists of particle ids. */
static std::unordered_map<int, std::unordered_set<int>> particle_type_map;

/** @brief Mapping particle ids to MPI ranks. */
static std::unordered_map<int, int> particle_node;

/**
 * @brief Keep track of the largest particle id.
 * This book-keeping variable is necessary to make particle insertion run
 * in constant time. Traversing the @ref particle_node to find the largest
 * particle id scales with O(N) and traversing the local cells in parallel
 * followed by a reduction scales with O(N^2).
 */
static int max_seen_pid = -1;

void init_type_map(int type) {
  type_list_enable = true;
  if (type < 0)
    throw std::runtime_error("Types may not be negative");

  auto &map_for_type = particle_type_map[type];
  map_for_type.clear();
  for (auto const &p : partCfg()) {
    if (p.type() == type)
      map_for_type.insert(p.id());
  }
}

static void remove_id_from_map(int p_id, int type) {
  auto it = particle_type_map.find(type);
  if (it != particle_type_map.end())
    it->second.erase(p_id);
}

static void add_id_to_type_map(int p_id, int type) {
  auto it = particle_type_map.find(type);
  if (it != particle_type_map.end())
    it->second.insert(p_id);
}

void on_particle_type_change(int p_id, int type) {
  if (type_list_enable) {
    // check if the particle exists already and the type is changed, then remove
    // it from the list which contains it
    auto const &cur_par = get_particle_data(p_id);
    int prev_type = cur_par.type();
    if (prev_type != type) {
      // particle existed before so delete it from the list
      remove_id_from_map(p_id, prev_type);
    }
    add_id_to_type_map(p_id, type);
  }
}

namespace {
/* Limit cache to 100 MiB */
std::size_t const max_cache_size = (100ul * 1048576ul) / sizeof(Particle);
Utils::Cache<int, Particle> particle_fetch_cache(max_cache_size);
} // namespace

void invalidate_fetch_cache() { particle_fetch_cache.invalidate(); }
std::size_t fetch_cache_max_size() { return particle_fetch_cache.max_size(); }

static boost::optional<const Particle &> get_particle_data_local(int p_id) {
  auto p = cell_structure.get_local_particle(p_id);

  if (p and (not p->is_ghost())) {
    return *p;
  }

  return {};
}

REGISTER_CALLBACK_ONE_RANK(get_particle_data_local)

const Particle &get_particle_data(int p_id) {
  auto const pnode = get_particle_node(p_id);

  if (pnode == this_node) {
    assert(cell_structure.get_local_particle(p_id));
    return *cell_structure.get_local_particle(p_id);
  }

  /* Query the cache */
  auto const p_ptr = particle_fetch_cache.get(p_id);
  if (p_ptr) {
    return *p_ptr;
  }

  /* Cache miss, fetch the particle,
   * put it into the cache and return a pointer into the cache. */
  auto const cache_ptr = particle_fetch_cache.put(
      p_id, Communication::mpiCallbacks().call(Communication::Result::one_rank,
                                               get_particle_data_local, p_id));
  return *cache_ptr;
}

static void mpi_get_particles_local() {
  std::vector<int> ids;
  boost::mpi::scatter(comm_cart, ids, 0);

  std::vector<Particle> parts(ids.size());
  std::transform(ids.begin(), ids.end(), parts.begin(), [](int id) {
    assert(cell_structure.get_local_particle(id));
    return *cell_structure.get_local_particle(id);
  });

  Utils::Mpi::gatherv(comm_cart, parts.data(), static_cast<int>(parts.size()),
                      0);
}

REGISTER_CALLBACK(mpi_get_particles_local)

/**
 * @brief Get multiple particles at once.
 *
 * *WARNING* Particles are returned in an arbitrary order.
 *
 * @param ids The ids of the particles that should be returned.
 *
 * @returns The particle list.
 */
static std::vector<Particle> mpi_get_particles(Utils::Span<const int> ids) {
  mpi_call(mpi_get_particles_local);
  /* Return value */
  std::vector<Particle> parts(ids.size());

  /* Group ids per node */
  static std::vector<std::vector<int>> node_ids(comm_cart.size());
  for (auto &per_node : node_ids) {
    per_node.clear();
  }

  for (auto const &p_id : ids) {
    auto const p_node = get_particle_node(p_id);

    if (p_node != this_node)
      node_ids[p_node].push_back(p_id);
  }

  /* Distributed ids to the nodes */
  {
    static std::vector<int> ignore;
    boost::mpi::scatter(comm_cart, node_ids, ignore, 0);
  }

  /* Copy local particles */
  std::transform(node_ids[this_node].cbegin(), node_ids[this_node].cend(),
                 parts.begin(), [](int p_id) {
                   assert(cell_structure.get_local_particle(p_id));
                   return *cell_structure.get_local_particle(p_id);
                 });

  static std::vector<int> node_sizes(comm_cart.size());
  std::transform(
      node_ids.cbegin(), node_ids.cend(), node_sizes.begin(),
      [](std::vector<int> const &ids) { return static_cast<int>(ids.size()); });

  Utils::Mpi::gatherv(comm_cart, parts.data(), static_cast<int>(parts.size()),
                      parts.data(), node_sizes.data(), 0);

  return parts;
}

void prefetch_particle_data(Utils::Span<const int> in_ids) {
  /* Nothing to do on a single node. */
  // NOLINTNEXTLINE(clang-analyzer-core.NonNullParamChecker)
  if (comm_cart.size() == 1)
    return;

  static std::vector<int> ids;
  ids.clear();
  auto out_ids = std::back_inserter(ids);

  std::copy_if(in_ids.begin(), in_ids.end(), out_ids, [](int id) {
    return (get_particle_node(id) != this_node) && particle_fetch_cache.has(id);
  });

  /* Don't prefetch more particles than fit the cache. */
  if (ids.size() > particle_fetch_cache.max_size())
    ids.resize(particle_fetch_cache.max_size());

  /* Fetch the particles... */
  for (auto &p : mpi_get_particles(ids)) {
    auto id = p.id();
    particle_fetch_cache.put(id, std::move(p));
  }
}

static void mpi_who_has_local() {
  static std::vector<int> sendbuf;

  auto local_particles = cell_structure.local_particles();
  auto const n_part = static_cast<int>(local_particles.size());
  boost::mpi::gather(comm_cart, n_part, 0);

  if (n_part == 0)
    return;

  sendbuf.resize(n_part);

  std::transform(local_particles.begin(), local_particles.end(),
                 sendbuf.begin(), [](Particle const &p) { return p.id(); });

  comm_cart.send(0, some_tag, sendbuf);
}

REGISTER_CALLBACK(mpi_who_has_local)

static void mpi_who_has() {
  mpi_call(mpi_who_has_local);

  auto local_particles = cell_structure.local_particles();

  static std::vector<int> n_parts;
  boost::mpi::gather(comm_cart, static_cast<int>(local_particles.size()),
                     n_parts, 0);

  static std::vector<int> pdata;
  max_seen_pid = -1;

  /* then fetch particle locations */
  for (int pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      for (auto const &p : local_particles) {
        particle_node[p.id()] = this_node;
        max_seen_pid = std::max(max_seen_pid, p.id());
      }
    } else if (n_parts[pnode] > 0) {
      pdata.resize(n_parts[pnode]);
      comm_cart.recv(pnode, some_tag, pdata);
      for (int i = 0; i < n_parts[pnode]; i++) {
        particle_node[pdata[i]] = pnode;
        max_seen_pid = std::max(max_seen_pid, pdata[i]);
      }
    }
  }
}

/**
 * @brief Rebuild the particle index.
 */
static void build_particle_node() { mpi_who_has(); }

int get_particle_node(int p_id) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }

  if (particle_node.empty())
    build_particle_node();

  auto const needle = particle_node.find(p_id);

  // Check if particle has a node, if not, we assume it does not exist.
  if (needle == particle_node.end()) {
    throw std::runtime_error("Particle node for id " + std::to_string(p_id) +
                             " not found!");
  }
  return needle->second;
}

void clear_particle_node() { particle_node.clear(); }

/**
 * @brief Calculate the largest particle id.
 * Traversing the @ref particle_node to find the largest particle id
 * scales with O(N). Consider using the cached value in @ref max_seen_pid
 * if possible. This function is only necessary when the cached value is
 * invalidated, for example when removing the particle which has the
 * largest id.
 */
static int calculate_max_seen_id() {
  return boost::accumulate(particle_node, -1,
                           [](int max, const std::pair<int, int> &kv) {
                             return std::max(max, kv.first);
                           });
}

/**
 * @brief Create a new particle and attach it to a cell.
 * @param p_id  The identity of the particle to create.
 * @param pos   The particle position.
 * @return Whether the particle was created on that node.
 */
static bool maybe_insert_particle(int p_id, Utils::Vector3d const &pos) {
  auto folded_pos = pos;
  auto image_box = Utils::Vector3i{};
  fold_position(folded_pos, image_box, box_geo);

  Particle new_part;
  new_part.id() = p_id;
  new_part.pos() = folded_pos;
  new_part.image_box() = image_box;

  return ::cell_structure.add_local_particle(std::move(new_part)) != nullptr;
}

/**
 * @brief Move particle to a new position.
 * @param p_id  The identity of the particle to move.
 * @param pos   The new particle position.
 * @return Whether the particle was moved from that node.
 */
static bool maybe_move_particle(int p_id, Utils::Vector3d const &pos) {
  auto pt = ::cell_structure.get_local_particle(p_id);
  if (pt == nullptr) {
    return false;
  }
  auto folded_pos = pos;
  auto image_box = Utils::Vector3i{};
  fold_position(folded_pos, image_box, box_geo);
  pt->pos() = folded_pos;
  pt->image_box() = image_box;
  return true;
}

static void particle_checks(int p_id, Utils::Vector3d const &pos) {
  if (p_id < 0) {
    throw std::domain_error("Invalid particle id: " + std::to_string(p_id));
  }
#ifndef __FAST_MATH__
  if (std::isnan(pos[0]) or std::isnan(pos[1]) or std::isnan(pos[2]) or
      std::isinf(pos[0]) or std::isinf(pos[1]) or std::isinf(pos[2])) {
    throw std::domain_error("Particle position must be finite");
  }
#endif // __FAST_MATH__
}

static void mpi_remove_particle_local(int p_id) {
  cell_structure.remove_particle(p_id);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_remove_particle_local)

static void mpi_remove_all_particles_local() {
  cell_structure.remove_all_particles();
  on_particle_change();
}

REGISTER_CALLBACK(mpi_remove_all_particles_local)

void remove_all_particles() {
  mpi_call_all(mpi_remove_all_particles_local);
  clear_particle_node();
}

void remove_particle(int p_id) {
  auto const &cur_par = get_particle_data(p_id);
  if (type_list_enable) {
    // remove particle from its current type_list
    int type = cur_par.type();
    remove_id_from_map(p_id, type);
  }

  particle_node[p_id] = -1;
  mpi_call_all(mpi_remove_particle_local, p_id);
  particle_node.erase(p_id);
  if (p_id == max_seen_pid) {
    --max_seen_pid;
    // if there is a gap (i.e. there is no particle with id max_seen_pid - 1,
    // then the cached value is invalidated and has to be recomputed (slow)
    if (particle_node.count(max_seen_pid) == 0 or
        particle_node[max_seen_pid] == -1) {
      max_seen_pid = calculate_max_seen_id();
    }
  }
}

static void mpi_make_new_particle_local(int p_id, Utils::Vector3d const &pos) {
  auto const has_created = maybe_insert_particle(p_id, pos);
  on_particle_change();

  auto node = -1;
  auto const node_local = (has_created) ? ::comm_cart.rank() : 0;
  boost::mpi::reduce(::comm_cart, node_local, node, std::plus<int>{}, 0);
  if (::this_node == 0) {
    particle_node[p_id] = node;
    max_seen_pid = std::max(max_seen_pid, p_id);
    assert(not has_created or node == 0);
  }
}

REGISTER_CALLBACK(mpi_make_new_particle_local)

void mpi_make_new_particle(int p_id, Utils::Vector3d const &pos) {
  particle_checks(p_id, pos);
  mpi_call_all(mpi_make_new_particle_local, p_id, pos);
}

static void mpi_set_particle_pos_local(int p_id, Utils::Vector3d const &pos) {
  auto const has_moved = maybe_move_particle(p_id, pos);
  ::cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();

  auto success = false;
  boost::mpi::reduce(::comm_cart, has_moved, success, std::plus<bool>{}, 0);
  if (::this_node == 0 and !success) {
    throw std::runtime_error("Particle node for id " + std::to_string(p_id) +
                             " not found!");
  }
}

REGISTER_CALLBACK(mpi_set_particle_pos_local)

void mpi_set_particle_pos(int p_id, Utils::Vector3d const &pos) {
  particle_checks(p_id, pos);
  mpi_call_all(mpi_set_particle_pos_local, p_id, pos);
}

int get_random_p_id(int type, int random_index_in_type_map) {
  auto it = particle_type_map.find(type);
  if (it == particle_type_map.end()) {
    throw std::runtime_error("The provided particle type " +
                             std::to_string(type) +
                             " is currently not tracked by the system.");
  }

  if (random_index_in_type_map + 1 > it->second.size())
    throw std::runtime_error("The provided index exceeds the number of "
                             "particle types listed in the particle_type_map");
  return *std::next(it->second.begin(), random_index_in_type_map);
}

int number_of_particles_with_type(int type) {
  auto it = particle_type_map.find(type);
  if (it == particle_type_map.end()) {
    throw std::runtime_error("The provided particle type " +
                             std::to_string(type) +
                             " is currently not tracked by the system.");
  }

  return static_cast<int>(it->second.size());
}

bool particle_exists(int p_id) {
  if (particle_node.empty())
    build_particle_node();
  return particle_node.count(p_id);
}

std::vector<int> get_particle_ids() {
  if (particle_node.empty())
    build_particle_node();

  auto ids = Utils::keys(particle_node);
  boost::sort(ids);

  return ids;
}

int get_maximal_particle_id() {
  if (particle_node.empty())
    build_particle_node();

  return max_seen_pid;
}

int get_n_part() {
  if (particle_node.empty())
    build_particle_node();

  return static_cast<int>(particle_node.size());
}

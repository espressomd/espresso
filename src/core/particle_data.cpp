/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

  This file is part of ESPResSo.

  ESPResSo is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 3 of the License, or
  (at your option) any later version.

  ESPResSo is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
/** \file particle_data.cpp
    This file contains everything related to particle storage. If you want to
   add a new property to the particles, it is probably a good idea to modify
   \ref Particle to give scripts access to that property. You always have to
   modify two positions: first the print section, where you should add your new
   data at the end, and second the read section where you have to find a nice
   and short name for your property to appear in the Tcl code. Then you just
   parse your part out of argc and argv.

    The corresponding header file is particle_data.hpp.
*/
#include "particle_data.hpp"
#include "PartCfg.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "interaction_data.hpp"
#include "partCfg_global.hpp"
#include "rotation.hpp"
#include "virtual_sites.hpp"

#include "utils.hpp"
#include "utils/Cache.hpp"
#include "utils/make_unique.hpp"
#include "utils/mpi/gatherv.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <unordered_map>
#include <unordered_set>
/************************************************
 * defines
 ************************************************/

/** granularity of the particle buffers in particles */
#define PART_INCREMENT 8

/** my magic MPI code for send/recv_particles */
#define REQ_SNDRCV_PART 0xaa

/************************************************
 * variables
 ************************************************/
// List of particles for grandcanonical simulations
int number_of_type_lists;
bool type_list_enable;
std::unordered_map<int, std::unordered_set<int>> particle_type_map{};
void remove_id_from_map(int part_id, int type);
void add_id_to_type_map(int part_id, int type);

int max_seen_particle = -1;
int n_part = 0;
/**
 * @brief id -> rank
 */
std::unordered_map<int, int> particle_node;

int max_local_particles = 0;
Particle **local_particles = nullptr;

/************************************************
 * local functions
 ************************************************/

/** Remove bond from particle if possible */
int try_delete_bond(Particle *part, int *bond);

/** Remove exclusion from particle if possible */
void try_delete_exclusion(Particle *part, int part2);

/** Insert an exclusion if not already set */
void try_add_exclusion(Particle *part, int part2);

/** Automatically add the next \<distance\> neighbors in each molecule to the
   exclusion list. This uses the bond topology obtained directly from the
   particles, since only this contains the full topology, in contrast to \ref
   topology::topology. To easily setup the bonds, all data should be on a single
   node, therefore the \ref partCfg array is used. With large amounts of
   particles, you should avoid this function and setup exclusions manually. */
void auto_exclusion(int distance);

/************************************************
 * particle initialization functions
 ************************************************/

/** Deallocate the dynamic storage of a particle. */
void free_particle(Particle *part) { part->~Particle(); }

void mpi_who_has_slave(int node, int param) {
  static int *sendbuf;
  int n_part;

  n_part = cells_get_n_particles();
  MPI_Gather(&n_part, 1, MPI_INT, nullptr, 0, MPI_INT, 0, comm_cart);
  if (n_part == 0)
    return;

  sendbuf = Utils::realloc(sendbuf, sizeof(int) * n_part);

  auto end = std::transform(local_cells.particles().begin(),
                            local_cells.particles().end(), sendbuf,
                            [](Particle const &p) { return p.p.identity; });

  auto npart = std::distance(sendbuf, end);
  MPI_Send(sendbuf, npart, MPI_INT, 0, SOME_TAG, comm_cart);
}

void mpi_who_has() {
  static int *sizes = new int[n_nodes];
  int *pdata = nullptr;
  int pdata_s = 0;

  mpi_call(mpi_who_has_slave, -1, 0);

  int n_part = cells_get_n_particles();
  /* first collect number of particles on each node */
  MPI_Gather(&n_part, 1, MPI_INT, sizes, 1, MPI_INT, 0, comm_cart);

  /* then fetch particle locations */
  for (int pnode = 0; pnode < n_nodes; pnode++) {
    COMM_TRACE(
        fprintf(stderr, "node %d reports %d particles\n", pnode, sizes[pnode]));
    if (pnode == this_node) {
      for (auto const &p : local_cells.particles())
        particle_node[p.p.identity] = this_node;

    } else if (sizes[pnode] > 0) {
      if (pdata_s < sizes[pnode]) {
        pdata_s = sizes[pnode];
        pdata = Utils::realloc(pdata, sizeof(int) * pdata_s);
      }
      MPI_Recv(pdata, sizes[pnode], MPI_INT, pnode, SOME_TAG, comm_cart,
               MPI_STATUS_IGNORE);
      for (int i = 0; i < sizes[pnode]; i++)
        particle_node[pdata[i]] = pnode;
    }
  }
  free(pdata);
}

/**
 * @brief Rebuild the particle index.
 */
void build_particle_node() { mpi_who_has(); }

/**
 *  @brief Get the mpi rank which owns the particle with id.
*/
int get_particle_node(int id) {
  if ((id < 0) or (id > max_seen_particle))
    throw std::runtime_error("Invalid particle id!");

  if (particle_node.empty())
    build_particle_node();

  auto const needle = particle_node.find(id);

  // Check if particle has a node, if not, we assume it does not exist.
  if (needle == particle_node.end()) {
    throw std::runtime_error("Particle node not found!");
  } else {
    return needle->second;
  }
}

void clear_particle_node() { particle_node.clear(); }

/************************************************
 * organizational functions
 ************************************************/

/** resize \ref local_particles.
    \param part the highest existing particle
*/
void realloc_local_particles(int part) {
  if (part >= max_local_particles) {
    /* round up part + 1 in granularity PART_INCREMENT */
    max_local_particles =
        PART_INCREMENT * ((part + PART_INCREMENT) / PART_INCREMENT);
    local_particles = Utils::realloc(local_particles,
                                     sizeof(Particle *) * max_local_particles);

    /* Set new memory to 0 */
    for (int i = (max_seen_particle + 1); i < max_local_particles; i++)
      local_particles[i] = nullptr;
  }
}

void init_particlelist(ParticleList *pList) {
  pList->n = 0;
  pList->max = 0;
  pList->part = nullptr;
}

int realloc_particlelist(ParticleList *l, int size) {
  int old_max = l->max;
  Particle *old_start = l->part;

  PART_TRACE(fprintf(stderr, "%d: realloc_particlelist %p: %d/%d->%d\n",
                     this_node, (void *)l, l->n, l->max, size));

  if (size < l->max) {
    if (size == 0)
      /* to be able to free an array again */
      l->max = 0;
    else
      /* shrink not as fast, just lose half, rounded up */
      l->max =
          PART_INCREMENT *
          (((l->max + size + 1) / 2 + PART_INCREMENT - 1) / PART_INCREMENT);
  } else
    /* round up */
    l->max = PART_INCREMENT * ((size + PART_INCREMENT - 1) / PART_INCREMENT);
  if (l->max != old_max)
    l->part = Utils::realloc(l->part, sizeof(Particle) * l->max);
  return l->part != old_start;
}

void update_local_particles(ParticleList *pl) {
  Particle *p = pl->part;
  int n = pl->n, i;
  for (i = 0; i < n; i++)
    local_particles[p[i].p.identity] = &p[i];
}

Particle *got_particle(ParticleList *l, int id) {
  int i;

  for (i = 0; i < l->n; i++)
    if (l->part[i].p.identity == id)
      break;
  if (i == l->n)
    return nullptr;
  return &(l->part[i]);
}

void append_unindexed_particle(ParticleList *l, Particle &&part) {
  realloc_particlelist(l, ++l->n);
  new (&(l->part[l->n - 1])) Particle(std::move(part));
}

Particle *append_indexed_particle(ParticleList *l, Particle &&part) {
  auto const re = realloc_particlelist(l, ++l->n);
  auto p = new (&(l->part[l->n - 1])) Particle(std::move(part));

  if (re)
    update_local_particles(l);
  else
    local_particles[p->p.identity] = p;
  return p;
}

Particle *move_unindexed_particle(ParticleList *dl, ParticleList *sl, int i) {
  realloc_particlelist(dl, ++dl->n);
  auto dst = &dl->part[dl->n - 1];
  auto src = &sl->part[i];
  auto end = &sl->part[sl->n - 1];

  new (dst) Particle(std::move(*src));
  if (src != end) {
    new (src) Particle(std::move(*end));
  }

  sl->n -= 1;
  realloc_particlelist(sl, sl->n);
  return dst;
}

Particle *move_indexed_particle(ParticleList *dl, ParticleList *sl, int i) {
  int re = realloc_particlelist(dl, ++dl->n);
  Particle *dst = &dl->part[dl->n - 1];
  Particle *src = &sl->part[i];
  Particle *end = &sl->part[sl->n - 1];

  new (dst) Particle(std::move(*src));
  if (re) {
    update_local_particles(dl);
  } else {
    local_particles[dst->p.identity] = dst;
  }
  if (src != end) {
    new (src) Particle(std::move(*end));
  }
  if (realloc_particlelist(sl, --sl->n)) {
    update_local_particles(sl);
  } else if (src != end) {
    local_particles[src->p.identity] = src;
  }
  return dst;
}

namespace {
/* Limit cache to 100 MiB */
std::size_t const max_cache_size = (100ul * 1048576ul) / sizeof(Particle);
Utils::Cache<int, Particle> particle_fetch_cache(max_cache_size);
}

void invalidate_fetch_cache() { particle_fetch_cache.invalidate(); }

const Particle &get_particle_data(int part) {
  auto const pnode = get_particle_node(part);

  if (pnode == this_node) {
    assert(local_particles[part]);
    return *local_particles[part];
  }

  /* Query the cache */
  auto const p_ptr = particle_fetch_cache.get(part);
  if (p_ptr) {
    return *p_ptr;
  }

  /* Cache miss, fetch the particle,
  * put it into the cache and return a pointer into the cache. */
  auto const cache_ptr =
      particle_fetch_cache.put(part, mpi_recv_part(pnode, part));
  return *cache_ptr;
}

const Particle *get_particle_data_ptr(int part) {
  return &get_particle_data(part);
}

void mpi_get_particles_slave(int, int) {
  std::vector<int> ids;
  boost::mpi::scatter(comm_cart, ids, 0);

  std::vector<Particle> parts(ids.size());
  std::transform(ids.begin(), ids.end(), parts.begin(), [](int id) {
    assert(local_particles[id]);
    return *local_particles[id];
  });

  Utils::Mpi::gatherv(comm_cart, parts.data(), parts.size(), 0);
}

/**
 * @brief Get multiple particles at once.
 *
 * *WARNING* Particles are returned in an arbitrary order.
 *
 * @param ids The ids of the particles that should be returned.
 *
 * @returns The particle data.
 */
std::vector<Particle> mpi_get_particles(std::vector<int> const &ids) {
  mpi_call(mpi_get_particles_slave, 0, 0);
  /* Return value */
  std::vector<Particle> parts(ids.size());

  /* Group ids per node */
  std::vector<std::vector<int>> node_ids(comm_cart.size());
  for (auto const &id : ids) {
    auto const pnode = get_particle_node(id);

    node_ids[pnode].push_back(id);
  }

  /* Distributed ids to the nodes */
  {
    std::vector<int> ignore;
    boost::mpi::scatter(comm_cart, node_ids, ignore, 0);
  }

  /* Copy local particles */
  std::transform(node_ids[this_node].cbegin(), node_ids[this_node].cend(),
                 parts.begin(), [](int id) {
                   assert(id);
                   return *local_particles[id];
                 });

  std::vector<int> node_sizes(comm_cart.size());
  std::transform(
      node_ids.cbegin(), node_ids.cend(), node_sizes.begin(),
      [](std::vector<int> const &ids) { return static_cast<int>(ids.size()); });

  Utils::Mpi::gatherv(comm_cart, parts.data(), parts.size(), parts.data(),
                      node_sizes.data(), 0);

  return parts;
}

void prefetch_particle_data(std::vector<int> ids) {
  /* Nothing to do on a single node. */
  if (comm_cart.size() == 1)
    return;

  /* Remove local, already cached and non-existent particles from the list. */
  ids.erase(std::remove_if(ids.begin(), ids.end(),
                           [](int id) {
                             if (not particle_exists(id)) {
                               return true;
                             } else {
                               auto const pnode = get_particle_node(id);
                               return (pnode == this_node) ||
                                      particle_fetch_cache.has(id);
                             }
                           }),
            ids.end());

  /* Don't prefetch more particles than fit the cache. */
  if (ids.size() > particle_fetch_cache.max_size())
    ids.resize(particle_fetch_cache.max_size());

  /* Fetch the particles... */
  auto parts = mpi_get_particles(ids);

  /* mpi_get_particles does not return the parts in the correct
     order, so the ids need to be updated. */
  std::transform(parts.cbegin(), parts.cend(), ids.begin(),
                 [](Particle const &p) { return p.identity(); });

  /* ... and put them into the cache. */
  particle_fetch_cache.put(ids.cbegin(), ids.cend(),
                           std::make_move_iterator(parts.begin()));
}

int place_particle(int part, double p[3]) {
  int retcode = ES_PART_OK;

  int pnode;
  if (particle_exists(part)) {
    pnode = get_particle_node(part);
    mpi_place_particle(pnode, part, p);
  } else {
    /* new particle, node by spatial position */
    pnode = cell_structure.position_to_node(p);

    /* master node specific stuff */
    particle_node[part] = pnode;

    retcode = ES_PART_CREATED;

    mpi_place_new_particle(pnode, part, p);
  }

  return retcode;
}

int set_particle_v(int part, double v[3]) {
  auto const pnode = get_particle_node(part);

  mpi_send_v(pnode, part, v);
  return ES_OK;
}

#ifdef ENGINE
int set_particle_swimming(int part, ParticleParametersSwimming swim) {
  auto const pnode = get_particle_node(part);

  mpi_send_swimming(pnode, part, swim);
  return ES_OK;
}
#endif

int set_particle_f(int part, const Vector3d &F) {
  auto const pnode = get_particle_node(part);

  mpi_send_f(pnode, part, F);
  return ES_OK;
}

#ifdef SHANCHEN
int set_particle_solvation(int part, double *solvation) {
  auto const pnode = get_particle_node(part);

  mpi_send_solvation(pnode, part, solvation);
  return ES_OK;
}

#endif

#if defined(MASS) || defined(LB_BOUNDARIES_GPU)
int set_particle_mass(int part, double mass) {
  auto const pnode = get_particle_node(part);

  mpi_send_mass(pnode, part, mass);
  return ES_OK;
}
#else
constexpr double ParticleProperties::mass;
#endif

#ifdef ROTATIONAL_INERTIA
int set_particle_rotational_inertia(int part, double rinertia[3]) {
  auto const pnode = get_particle_node(part);

  mpi_send_rotational_inertia(pnode, part, rinertia);
  return ES_OK;
}
#endif
#ifdef ROTATION
int set_particle_rotation(int part, int rot) {
  auto const pnode = get_particle_node(part);

  mpi_send_rotation(pnode, part, rot);
  return ES_OK;
}
#endif
#ifdef ROTATION
int rotate_particle(int part, double axis[3], double angle) {
  auto const pnode = get_particle_node(part);

  mpi_rotate_particle(pnode, part, axis, angle);
  return ES_OK;
}
#endif

#ifdef AFFINITY
int set_particle_affinity(int part, double bond_site[3]) {
  auto const pnode = get_particle_node(part);

  mpi_send_affinity(pnode, part, bond_site);
  return ES_OK;
}
#endif

#ifdef MEMBRANE_COLLISION
int set_particle_out_direction(int part, double out_direction[3]) {
  auto const pnode = get_particle_node(part);

  mpi_send_out_direction(pnode, part, out_direction);
  return ES_OK;
}
#endif

#ifdef DIPOLES
int set_particle_dipm(int part, double dipm) {
  auto const pnode = get_particle_node(part);

  mpi_send_dipm(pnode, part, dipm);
  return ES_OK;
}

int set_particle_dip(int part, double dip[3]) {
  auto const pnode = get_particle_node(part);

  mpi_send_dip(pnode, part, dip);

  return ES_OK;
}

#endif

#ifdef VIRTUAL_SITES
int set_particle_virtual(int part, int is_virtual) {
  auto const pnode = get_particle_node(part);

  if (pnode == -1)
    return ES_ERROR;
  mpi_send_virtual(pnode, part, is_virtual);
  return ES_OK;
}
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void set_particle_vs_quat(int part, double *vs_quat) {
  auto const pnode = get_particle_node(part);
  mpi_send_vs_quat(pnode, part, vs_quat);
}

int set_particle_vs_relative(int part, int vs_relative_to, double vs_distance,
                             double *rel_ori) {
  auto const pnode = get_particle_node(part);

  // Send the stuff
  mpi_send_vs_relative(pnode, part, vs_relative_to, vs_distance, rel_ori);
  return ES_OK;
}
#endif

int set_particle_q(int part, double q) {
  auto const pnode = get_particle_node(part);

  mpi_send_q(pnode, part, q);
  return ES_OK;
}

#ifdef LB_ELECTROHYDRODYNAMICS
int set_particle_mu_E(int part, double mu_E[3]) {
  auto const pnode = get_particle_node(part);

  mpi_send_mu_E(pnode, part, mu_E);
  return ES_OK;
}

void get_particle_mu_E(int part, double (&mu_E)[3]) {
  auto const &p = get_particle_data(part);

  for (int i = 0; i < 3; i++) {
    mu_E[i] = p.p.mu_E[i];
  }
}
#endif

int set_particle_type(int p_id, int type) {
  auto const pnode = get_particle_node(p_id);
  make_particle_type_exist(type);

  if (type_list_enable) {
    // check if the particle exists already and the type is changed, then remove
    // it from the list which contains it
    auto const &cur_par = get_particle_data(p_id);
    int prev_type = cur_par.p.type;
    if (prev_type != type ) {
      // particle existed before so delete it from the list
      remove_id_from_map(p_id, prev_type);
    }
    add_id_to_type_map(p_id, type);
  }

  mpi_send_type(pnode, p_id, type);

  return ES_OK;
}

int set_particle_mol_id(int part, int mid) {
  auto const pnode = get_particle_node(part);

  mpi_send_mol_id(pnode, part, mid);
  return ES_OK;
}

#ifdef ROTATION
int set_particle_quat(int part, double quat[4]) {
  auto const pnode = get_particle_node(part);

  mpi_send_quat(pnode, part, quat);
  return ES_OK;
}

int set_particle_omega_lab(int part, double omega_lab[3]) {
  auto const &particle = get_particle_data(part);

  /* Internal functions require the body coordinates
     so we need to convert to these from the lab frame */

  double A[9];
  double omega[3];

  define_rotation_matrix(particle, A);

  omega[0] = A[0 + 3 * 0] * omega_lab[0] + A[0 + 3 * 1] * omega_lab[1] +
             A[0 + 3 * 2] * omega_lab[2];
  omega[1] = A[1 + 3 * 0] * omega_lab[0] + A[1 + 3 * 1] * omega_lab[1] +
             A[1 + 3 * 2] * omega_lab[2];
  omega[2] = A[2 + 3 * 0] * omega_lab[0] + A[2 + 3 * 1] * omega_lab[1] +
             A[2 + 3 * 2] * omega_lab[2];

  auto const pnode = get_particle_node(part);
  mpi_send_omega(pnode, part, omega);
  return ES_OK;
}

int set_particle_omega_body(int part, double omega[3]) {
  auto const pnode = get_particle_node(part);
  mpi_send_omega(pnode, part, omega);
  return ES_OK;
}

int set_particle_torque_lab(int part, double torque_lab[3]) {
  auto const &particle = get_particle_data(part);

  /* Internal functions require the body coordinates
     so we need to convert to these from the lab frame */

  double A[9];
  double torque[3];

  define_rotation_matrix(particle, A);

  torque[0] = A[0 + 3 * 0] * torque_lab[0] + A[0 + 3 * 1] * torque_lab[1] +
              A[0 + 3 * 2] * torque_lab[2];
  torque[1] = A[1 + 3 * 0] * torque_lab[0] + A[1 + 3 * 1] * torque_lab[1] +
              A[1 + 3 * 2] * torque_lab[2];
  torque[2] = A[2 + 3 * 0] * torque_lab[0] + A[2 + 3 * 1] * torque_lab[1] +
              A[2 + 3 * 2] * torque_lab[2];

  auto const pnode = get_particle_node(part);
  mpi_send_torque(pnode, part, torque);
  return ES_OK;
}

int set_particle_torque_body(int part, double torque[3]) {
  auto const pnode = get_particle_node(part);

  /* Nothing to be done but pass, since the coordinates
     are already in the proper frame */

  mpi_send_torque(pnode, part, torque);
  return ES_OK;
}

#endif

#ifdef LANGEVIN_PER_PARTICLE
int set_particle_temperature(int part, double T) {
  auto const pnode = get_particle_node(part);

  mpi_set_particle_temperature(pnode, part, T);
  return ES_OK;
}

#ifndef PARTICLE_ANISOTROPY
int set_particle_gamma(int part, double gamma)
#else
int set_particle_gamma(int part, Vector3d gamma)
#endif // PARTICLE_ANISOTROPY
{
  auto const pnode = get_particle_node(part);

  mpi_set_particle_gamma(pnode, part, gamma);
  return ES_OK;
}
#ifdef ROTATION
#ifndef PARTICLE_ANISOTROPY
int set_particle_gamma_rot(int part, double gamma_rot)
#else
int set_particle_gamma_rot(int part, Vector3d gamma_rot)
#endif // PARTICLE_ANISOTROPY
{
  auto const pnode = get_particle_node(part);

  mpi_set_particle_gamma_rot(pnode, part, gamma_rot);
  return ES_OK;
}
#endif // ROTATION
#endif // LANGEVIN_PER_PARTICLE

#ifdef EXTERNAL_FORCES
#ifdef ROTATION
int set_particle_ext_torque(int part, int flag, double torque[3]) {
  auto const pnode = get_particle_node(part);

  mpi_send_ext_torque(pnode, part, flag, PARTICLE_EXT_TORQUE, torque);
  return ES_OK;
}
#endif

int set_particle_ext_force(int part, int flag, double force[3]) {
  auto const pnode = get_particle_node(part);

  mpi_send_ext_force(pnode, part, flag, PARTICLE_EXT_FORCE, force);
  return ES_OK;
}

int set_particle_fix(int part, int flag) {
  auto const pnode = get_particle_node(part);

  mpi_send_ext_force(pnode, part, flag, COORDS_FIX_MASK, nullptr);
  return ES_OK;
}

#endif

int change_particle_bond(int part, int *bond, int _delete) {
  auto const pnode = get_particle_node(part);

  if (_delete != 0 || bond == nullptr)
    _delete = 1;

  if (bond != nullptr) {
    if (bond[0] < 0 || bond[0] >= bonded_ia_params.size()) {
      runtimeErrorMsg() << "invalid/unknown bonded interaction type "
                        << bond[0];
      return ES_ERROR;
    }
  }
  return mpi_send_bond(pnode, part, bond, _delete);
}

void remove_all_particles() {
  mpi_remove_particle(-1, -1);
  clear_particle_node();
}

int remove_particle(int p_id) {
  auto const &cur_par = get_particle_data(p_id);
  if (type_list_enable == true) {
    // remove particle from its current type_list
    int type = cur_par.p.type;
    remove_id_from_map(p_id, type);
  }

  auto const pnode = get_particle_node(p_id);

  particle_node[p_id] = -1;
  mpi_remove_particle(pnode, p_id);

  particle_node.erase(p_id);

  if (p_id == max_seen_particle) {
    max_seen_particle--;
    mpi_bcast_parameter(FIELD_MAXPART);
  }
  return ES_OK;
}

void local_remove_particle(int part) {
  int ind, c;
  Particle *p = local_particles[part];
  ParticleList *pl = nullptr, *tmp;

  /* the tricky - say ugly - part: determine
     the cell the particle is located in by checking
     wether the particle address is inside the array */
  for (c = 0; c < local_cells.n; c++) {
    tmp = local_cells.cell[c];
    ind = p - tmp->part;
    if (ind >= 0 && ind < tmp->n) {
      pl = tmp;
      break;
    }
  }
  if (!pl) {
    fprintf(stderr,
            "%d: INTERNAL ERROR: could not find cell of particle %d, exiting\n",
            this_node, part);
    errexit();
  }

  free_particle(p);

  /* remove local_particles entry */
  local_particles[p->p.identity] = nullptr;

  if (&pl->part[pl->n - 1] != p) {
    /* move last particle to free position */
    *p = pl->part[pl->n - 1];

    /* update the local_particles array for the moved particle */
    local_particles[p->p.identity] = p;
  }
  pl->n--;
}

void local_place_particle(int part, const double p[3], int _new) {
  Cell *cell;
  double pp[3];
  int i[3], rl;
  Particle *pt;

  i[0] = 0;
  i[1] = 0;
  i[2] = 0;
  pp[0] = p[0];
  pp[1] = p[1];
  pp[2] = p[2];

  double vv[3] = {0., 0., 0.};
  fold_position(pp, vv, i);

  if (_new) {
    /* allocate particle anew */
    cell = cell_structure.position_to_cell(pp);
    if (!cell) {
      fprintf(stderr, "%d: INTERNAL ERROR: particle %d at %f(%f) %f(%f) %f(%f) "
                      "does not belong on this node\n",
              this_node, part, p[0], pp[0], p[1], pp[1], p[2], pp[2]);
      errexit();
    }
    rl = realloc_particlelist(cell, ++cell->n);
    pt = new (&cell->part[cell->n - 1]) Particle;

    pt->p.identity = part;
    if (rl)
      update_local_particles(cell);
    else
      local_particles[pt->p.identity] = pt;
  } else
    pt = local_particles[part];

  PART_TRACE(fprintf(
      stderr, "%d: local_place_particle: got particle id=%d @ %f %f %f\n",
      this_node, part, p[0], p[1], p[2]));

  memmove(pt->r.p.data(), pp, 3 * sizeof(double));
  memmove(pt->l.i.data(), i, 3 * sizeof(int));
#ifdef BOND_CONSTRAINT
  memmove(pt->r.p_old.data(), pp, 3 * sizeof(double));
#endif
}

void local_remove_all_particles() {
  Cell *cell;
  int c;
  n_part = 0;
  max_seen_particle = -1;
  std::fill(local_particles, local_particles + max_local_particles, nullptr);

  for (c = 0; c < local_cells.n; c++) {
    Particle *p;
    int i, np;
    cell = local_cells.cell[c];
    p = cell->part;
    np = cell->n;
    for (i = 0; i < np; i++)
      free_particle(&p[i]);
    cell->n = 0;
  }
}

void local_rescale_particles(int dir, double scale) {
  for (auto &p : local_cells.particles()) {
    if (dir < 3)
      p.r.p[dir] *= scale;
    else {
      p.r.p[0] *= scale;
      p.r.p[1] *= scale;
      p.r.p[2] *= scale;
    }
  }
}

void added_particle(int part) {
  n_part++;

  if (part > max_seen_particle) {
    realloc_local_particles(part);

    max_seen_particle = part;
  }
}

int local_change_bond(int part, int *bond, int _delete) {
  auto p = local_particles[part];
  if (_delete)
    return try_delete_bond(p, bond);

  auto const bond_size = bonded_ia_params[bond[0]].num + 1;

  std::copy_n(bond, bond_size, std::back_inserter(p->bl));

  return ES_OK;
}

int try_delete_bond(Particle *part, int *bond) {
  IntList *bl = &part->bl;
  int i, j, type, partners;

  // Empty bond means: delete all bonds
  if (!bond) {
    bl->clear();

    return ES_OK;
  }

  // Go over the bond list to find the bond to delete
  for (i = 0; i < bl->n;) {
    type = bl->e[i];
    partners = bonded_ia_params[type].num;

    // If the bond type does not match the one, we want to delete, skip
    if (type != bond[0])
      i += 1 + partners;
    else {
      // Go over the bond partners
      for (j = 1; j <= partners; j++) {
        // Leave the loop early, if the bond to delete and the bond with in the
        // particle don't match
        if (bond[j] != bl->e[i + j])
          break;
      }
      // If we did not exit from the loop early, all parameters matched
      // and we go on with deleting
      if (j > partners) {
        // New length of bond list
        bl->erase(bl->begin() + i, bl->begin() + i + 1 + partners);

        return ES_OK;
      }
      i += 1 + partners;
    }
  }
  return ES_ERROR;
}

void remove_all_bonds_to(int identity) {
  for (auto &p : local_cells.particles()) {
    IntList *bl = &p.bl;
    int i, j, partners;

    for (i = 0; i < bl->n;) {
      partners = bonded_ia_params[bl->e[i]].num;
      for (j = 1; j <= partners; j++)
        if (bl->e[i + j] == identity)
          break;
      if (j <= partners) {
        bl->erase(bl->begin() + i, bl->begin() + i + 1 + partners);
      } else
        i += 1 + partners;
    }
    if (i != bl->n) {
      fprintf(stderr, "%d: INTERNAL ERROR: bond information corrupt for "
                      "particle %d, exiting...\n",
              this_node, p.p.identity);
      errexit();
    }
  }
}

#ifdef EXCLUSIONS
void local_change_exclusion(int part1, int part2, int _delete) {
  if (part1 == -1 && part2 == -1) {
    for (auto &p : local_cells.particles()) {
      p.el.clear();
    }

    return;
  }

  /* part1, if here */
  auto part = local_particles[part1];
  if (part) {
    if (_delete)
      try_delete_exclusion(part, part2);
    else
      try_add_exclusion(part, part2);
  }

  /* part2, if here */
  part = local_particles[part2];
  if (part) {
    if (_delete)
      try_delete_exclusion(part, part1);
    else
      try_add_exclusion(part, part1);
  }
}

void try_add_exclusion(Particle *part, int part2) {
  for (int i = 0; i < part->el.n; i++)
    if (part->el.e[i] == part2)
      return;

  part->el.push_back(part2);
}

void try_delete_exclusion(Particle *part, int part2) {
  IntList &el = part->el;

  el.erase(std::remove(el.begin(), el.end(), part2), el.end());
}
#endif

#include "utils/serialization/ParticleList.hpp"

void send_particles(ParticleList *particles, int node) {
  PART_TRACE(fprintf(stderr, "%d: send_particles %d to %d\n", this_node,
                     particles->n, node));

  comm_cart.send(node, REQ_SNDRCV_PART, *particles);

  /* remove particles from this nodes local list and free data */
  for (int pc = 0; pc < particles->n; pc++) {
    local_particles[particles->part[pc].p.identity] = nullptr;
    free_particle(&particles->part[pc]);
  }

  realloc_particlelist(particles, particles->n = 0);
}

void recv_particles(ParticleList *particles, int node) {
  PART_TRACE(fprintf(stderr, "%d: recv_particles from %d\n", this_node, node));
  comm_cart.recv(node, REQ_SNDRCV_PART, *particles);

  update_local_particles(particles);
}

#ifdef EXCLUSIONS

namespace {
/* keep a unique list for particle i. Particle j is only added if it is not i
   and not already in the list. */
void add_partner(IntList *il, int i, int j, int distance) {
  int k;
  if (j == i)
    return;
  for (k = 0; k < il->n; k += 2)
    if (il->e[k] == j)
      return;

  il->push_back(j);
  il->push_back(distance);
}
}

int change_exclusion(int part1, int part2, int _delete) {
  if (particle_exists(part1) && particle_exists(part2)) {
    mpi_send_exclusion(part1, part2, _delete);
    return ES_OK;
  } else {
    return ES_ERROR;
  }
}

void remove_all_exclusions() { mpi_send_exclusion(-1, -1, 1); }

void auto_exclusions(int distance) {
  int count, p, i, j, p1, p2, p3, dist1, dist2;
  Bonded_ia_parameters *ia_params;

  /* partners is a list containing the currently found excluded particles for
     each particle, and their distance, as a interleaved list */
  std::unordered_map<int, IntList> partners;

  /* We need bond information */
  partCfg().update_bonds();

  /* determine initial connectivity */
  for (auto const &part1 : partCfg()) {
    p1 = part1.p.identity;
    for (i = 0; i < part1.bl.n;) {
      ia_params = &bonded_ia_params[part1.bl.e[i++]];
      if (ia_params->num == 1) {
        p2 = part1.bl.e[i++];
        /* you never know what the user does, may bond a particle to itself...?
         */
        if (p2 != p1) {
          add_partner(&partners[p1], p1, p2, 1);
          add_partner(&partners[p2], p2, p1, 1);
        }
      } else
        i += ia_params->num;
    }
  }

  /* calculate transient connectivity. For each of the current neighbors,
     also exclude their close enough neighbors.
  */
  for (count = 1; count < distance; count++) {
    for (p1 = 0; p1 <= max_seen_particle; p1++) {
      for (i = 0; i < partners[p1].n; i += 2) {
        p2 = partners[p1].e[i];
        dist1 = partners[p1].e[i + 1];
        if (dist1 > distance)
          continue;
        /* loop over all partners of the partner */
        for (j = 0; j < partners[p2].n; j += 2) {
          p3 = partners[p2].e[j];
          dist2 = dist1 + partners[p2].e[j + 1];
          if (dist2 > distance)
            continue;
          add_partner(&partners[p1], p1, p3, dist2);
          add_partner(&partners[p3], p3, p1, dist2);
        }
      }
    }
  }

  /* setup the exclusions and clear the arrays. We do not setup the exclusions
     up there, since on_part_change clears the partCfg, so that we would have to
     restore it continously. Of course this could be optimized by bundling the
     exclusions, but this is only done once and the overhead is as much as for
     setting the bonds, which the user apparently accepted.
  */
  for (p = 0; p <= max_seen_particle; p++) {
    for (j = 0; j < partners[p].n; j++)
      if (p < partners[p].e[j])
        change_exclusion(p, partners[p].e[j], 0);
  }
}

#endif

void init_type_map(int type) {
  type_list_enable = true;
  if (type < 0)
    throw std::runtime_error("Types may not be negative");

  // fill particle map
  if (particle_type_map.count(type) == 0)
    particle_type_map[type] = std::unordered_set<int>();

  for (auto const &p : partCfg()) {
    if (p.p.type == type)
      particle_type_map.at(type).insert(p.p.identity);
  }
}

void remove_id_from_map(int part_id, int type) {
  if(particle_type_map.find(type)!=particle_type_map.end())
    particle_type_map.at(type).erase(part_id);
}

int get_random_p_id(int type) {
  if (particle_type_map.at(type).size() == 0)
    throw std::runtime_error("No particles of given type could be found");
  int rand_index = i_random(particle_type_map.at(type).size());
  return *std::next(particle_type_map[type].begin(), rand_index);
}

void add_id_to_type_map(int part_id, int type) {
  if(particle_type_map.find(type)!=particle_type_map.end())
    particle_type_map.at(type).insert(part_id);
}

int number_of_particles_with_type(int type) {
  return static_cast<int>(particle_type_map.at(type).size());
}

// The following functions are used by the python interface to obtain
// properties of a particle, which are only compiled in in some configurations
// This is needed, because cython does not support conditional compilation
// within a ctypedef definition

#ifdef ROTATION
void pointer_to_omega_body(Particle const *p, double const *&res) {
  res = p->m.omega.data();
}

void pointer_to_torque_lab(Particle const *p, double const *&res) {
  res = p->f.torque.data();
}

void pointer_to_quat(Particle const *p, double const *&res) { res = p->r.quat.data(); }

void pointer_to_quatu(Particle const *p, double const *&res) {
  res = p->r.quatu.data();
}
#endif

#ifdef ELECTROSTATICS
void pointer_to_q(Particle const *p, double const *&res) { res = &(p->p.q); }
#endif

#ifdef VIRTUAL_SITES
void pointer_to_virtual(Particle const *p, int const *&res) {
  res = &(p->p.is_virtual);
}
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void pointer_to_vs_quat(Particle const *p, double const *&res) {
  res = (p->p.vs_quat);
}

void pointer_to_vs_relative(Particle const *p, int const *&res1,
                            double const *&res2, double const *&res3) {
  res1 = &(p->p.vs_relative_to_particle_id);
  res2 = &(p->p.vs_relative_distance);
  res3 = (p->p.vs_relative_rel_orientation);
}
#endif

#ifdef DIPOLES
void pointer_to_dip(Particle const *p, double const *&res) { res = p->r.dip.data(); }

void pointer_to_dipm(Particle const *p, double const *&res) {
  res = &(p->p.dipm);
}
#endif

#ifdef EXTERNAL_FORCES
void pointer_to_ext_force(Particle const *p, int const *&res1,
                          double const *&res2) {
  res1 = &(p->p.ext_flag);
  res2 = p->p.ext_force.data();
}
#ifdef ROTATION
void pointer_to_ext_torque(Particle const *p, int const *&res1,
                           double const *&res2) {
  res1 = &(p->p.ext_flag);
  res2 = p->p.ext_torque.data();
}
#endif
void pointer_to_fix(Particle const *p, int const *&res) {
  res = &(p->p.ext_flag);
}
#endif

#ifdef LANGEVIN_PER_PARTICLE
void pointer_to_gamma(Particle const *p, double const *&res) {
#ifndef PARTICLE_ANISOTROPY
  res = &(p->p.gamma);
#else
  res = p->p.gamma.data(); // array [3]
#endif // PARTICLE_ANISTROPY
}

#ifdef ROTATION
void pointer_to_gamma_rot(Particle const *p, double const *&res) {
#ifndef PARTICLE_ANISOTROPY
  res = &(p->p.gamma_rot);
#else
  res = p->p.gamma_rot.data(); // array [3]
#endif // ROTATIONAL_INERTIA
}
#endif // ROTATION

void pointer_to_temperature(Particle const *p, double const *&res) {
  res = &(p->p.T);
}
#endif // LANGEVIN_PER_PARTICLE

void pointer_to_rotation(Particle const *p, short int const *&res) {
  res = &(p->p.rotation);
}

#ifdef ENGINE
void pointer_to_swimming(Particle const *p,
                         ParticleParametersSwimming const *&swim) {
  swim = &(p->swim);
}
#endif

#ifdef ROTATIONAL_INERTIA
void pointer_to_rotational_inertia(Particle const *p, double const *&res) {
  res = p->p.rinertia.data();
}
#endif

#ifdef AFFINITY
void pointer_to_bond_site(Particle const* p, double const*& res) {
  res =p->p.bond_site.data();
}
#endif

#ifdef MEMBRANE_COLLISION
void pointer_to_out_direction(const Particle* p, const double*& res) {
 res = p->p.out_direction.data();
}
#endif



bool particle_exists(int part_id) {
  if (particle_node.empty())
    build_particle_node();
  return particle_node.count(part_id);
}

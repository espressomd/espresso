/*
  Copyright (C) 2010-2018 The ESPResSo project
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
/** \file
 *  Particles and particle lists.
 *
 *  The corresponding header file is particle_data.hpp.
 */

#include "particle_data.hpp"

#include "PartCfg.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "debug.hpp"
#include "event.hpp"
#include "global.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "partCfg_global.hpp"
#include "random.hpp"
#include "rotation.hpp"
#include "serialization/ParticleList.hpp"
#include "virtual_sites.hpp"

#include <utils/Cache.hpp>
#include <utils/constants.hpp>
#include <utils/mpi/gatherv.hpp>

#include <boost/range/algorithm.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>

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

namespace {
/**
 * @brief A generic particle update.
 *
 * Here the sub-struct struture of Particle is
 * used: the specification of the data member to update
 * consists of two parts, the pointer to the substruct @p s
 * and a pointer to a member of that substruct @p m.
 *
 * @tparam S Substruct type of Particle
 * @tparam s Pointer to a member of Particle
 * @tparam T Type of the member to update, must be serializable
 * @tparam m Pointer to the member.
 */
template <typename S, S Particle::*s, typename T, T S::*m>
struct UpdateParticle {
  T value;

  void operator()(Particle &p) const { (p.*s).*m = value; }

  template <class Archive> void serialize(Archive &ar, long int) { ar &value; }
};

template <typename T, T ParticleProperties::*m>
using UpdateProperty = UpdateParticle<ParticleProperties, &Particle::p, T, m>;
template <typename T, T ParticlePosition ::*m>
using UpdatePosition = UpdateParticle<ParticlePosition, &Particle::r, T, m>;
template <typename T, T ParticleMomentum ::*m>
using UpdateMomentum = UpdateParticle<ParticleMomentum, &Particle::m, T, m>;
template <typename T, T ParticleForce ::*m>
using UpdateForce = UpdateParticle<ParticleForce, &Particle::f, T, m>;

#ifdef EXTERNAL_FORCES
/**
 * @brief Special updater for the external flags.
 *
 * These need to be treated specialy as they are
 * updated masked and not overwritten.
 */
struct UpdateExternalFlag {
  /* The bits to update */
  int mask;
  /* The actual values for the update */
  int flag;
  void operator()(Particle &p) const {
    /* mask out old flags */
    p.p.ext_flag &= ~mask;
    /* set new values */
    p.p.ext_flag |= (mask & flag);
  }

  template <class Archive> void serialize(Archive &ar, long int) {
    ar &mask;
    ar &flag;
  }
};
#endif

using Prop = ParticleProperties;

// clang-format off
using UpdatePropertyMessage = boost::variant
        < UpdateProperty<int, &Prop::type>
        , UpdateProperty<int, &Prop::mol_id>
#ifdef MASS
        , UpdateProperty<double, &Prop::mass>
#endif
#ifdef SHANCHEN
        , UpdateProperty<std::array<double, 2 * LB_COMPONENTS>, &Prop::solvation>
#endif
#ifdef ROTATIONAL_INERTIA
        , UpdateProperty<Utils::Vector3d, &Prop::rinertia>
#endif
#ifdef AFFINITY
        , UpdateProperty<Utils::Vector3d, &Prop::bond_site>
#endif
#ifdef MEMBRANE_COLLISION
        , UpdateProperty<Utils::Vector3d, &Prop::out_direction>
#endif
        , UpdateProperty<int, &Prop::rotation>
#ifdef ELECTROSTATICS
        , UpdateProperty<double, &Prop::q>
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
        , UpdateProperty<Utils::Vector3d, &Prop::mu_E>
#endif
#ifdef DIPOLES
        , UpdateProperty<double, &Prop::dipm>
#endif
#ifdef VIRTUAL_SITES
        , UpdateProperty<int, &Prop::is_virtual>
#ifdef VIRTUAL_SITES_RELATIVE
        , UpdateProperty<ParticleProperties::VirtualSitesRelativeParameteres,
                         &Prop::vs_relative>
#endif
#endif
#ifdef LANGEVIN_PER_PARTICLE
        , UpdateProperty<double, &Prop::T>
#ifndef PARTICLE_ANISOTROPY
        , UpdateProperty<double, &Prop::gamma>
#else
        , UpdateProperty<Utils::Vector3d, &Prop::gamma>
#endif // PARTICLE_ANISOTROPY
#ifdef ROTATION
#ifndef PARTICLE_ANISOTROPY
        , UpdateProperty<double, &Prop::gamma_rot>
#else
        , UpdateProperty<Utils::Vector3d, &Prop::gamma_rot>
#endif // ROTATIONAL_INERTIA
#endif // ROTATION
#endif // LANGEVIN_PER_PARTICLE
#ifdef EXTERNAL_FORCES
        , UpdateExternalFlag
        , UpdateProperty<Utils::Vector3d, &Prop::ext_force>
#ifdef ROTATION
        , UpdateProperty<Utils::Vector3d, &Prop::ext_torque>
#endif
#endif
        >;

using UpdatePositionMessage = boost::variant
        < UpdatePosition<Utils::Vector3d, &ParticlePosition::p>
#ifdef ROTATION
        , UpdatePosition<Utils::Vector4d, &ParticlePosition::quat>
#endif
        >;

using UpdateMomentumMessage = boost::variant
      < UpdateMomentum<Utils::Vector3d, &ParticleMomentum::v>
#ifdef ROTATION
      , UpdateMomentum<Utils::Vector3d, &ParticleMomentum::omega>
#endif
      >;

using UpdateForceMessage = boost::variant
      < UpdateForce<Utils::Vector3d, &ParticleForce::f>
#ifdef ROTATION
      , UpdateForce<Utils::Vector3d, &ParticleForce::torque>
#endif
      >;

/**
 * @brief Delete specific bond.
 */
struct RemoveBond {
    std::vector<int> bond;

    void operator()(Particle &p) const {
      try_delete_bond(&p, bond.data());
    }

    template<class Archive>
            void serialize(Archive &ar, long int) {
        ar & bond;
    }
};


/**
 * @brief Delete all bonds.
 */
struct RemoveBonds {
    void operator()(Particle &p) const {
      p.bl.clear();
    }

    template<class Archive>
    void serialize(Archive &ar, long int) {
    }
};

struct AddBond {
    std::vector<int> bond;

    void operator()(Particle &p) const {
        local_add_particle_bond(p, bond);
    }

    template<class Archive>
    void serialize(Archive &ar, long int) {
        ar & bond;
    }
};

using UpdateBondMessage = boost::variant
        < RemoveBond
        , RemoveBonds
        , AddBond
        >;

#ifdef ENGINE
struct UpdateSwim {
    ParticleParametersSwimming swim;

    void operator()(Particle &p) const {
      p.swim = swim;
    }

    template<class Archive>
    void serialize(Archive &ar, long int) {
      ar & swim;
    }
};
#endif

#ifdef ROTATION
struct UpdateOrientation {
    Utils::Vector3d axis;
    double angle;

    void operator()(Particle &p) const {
        local_rotate_particle(p, axis, angle);
    }

    template<class Archive>
    void serialize(Archive &ar, long int) {
        ar & axis & angle;
    }
};
#endif

/**
 * @brief Top-level message.
 *
 * A message is either updates a property,
 * or a position, or ...
 * New messages can be added here, if they
 * fullfill the type requirements, namely:
 * They either have an integer member id indicating
 * the particle that should be updated an operator()(const Particle&)
 * that is called with the particle, or a tree of
 * variants with leafs that have such a operator() and member.
 */
using UpdateMessage = boost::variant
        < UpdatePropertyMessage
        , UpdatePositionMessage
        , UpdateMomentumMessage
        , UpdateForceMessage
        , UpdateBondMessage
#ifdef ENGINE
        , UpdateSwim
#endif
#ifdef ROTATION
        , UpdateOrientation
#endif
        >;
// clang-format on

/**
 * @brief Meta-function to detect the message type from
 *        the particle substruct.
 */
template <typename S, S Particle::*s> struct message_type;

template <> struct message_type<ParticleProperties, &Particle::p> {
  using type = UpdatePropertyMessage;
};

template <> struct message_type<ParticlePosition, &Particle::r> {
  using type = UpdatePositionMessage;
};

template <> struct message_type<ParticleMomentum, &Particle::m> {
  using type = UpdateMomentumMessage;
};

template <> struct message_type<ParticleForce, &Particle::f> {
  using type = UpdateForceMessage;
};

template <typename S, S Particle::*s>
using message_type_t = typename message_type<S, s>::type;

/**
 * @brief Visitor for message evaluation.
 *
 * This visitor either recurses into the active type
 * if it is a variant type, or otherwise callse
 * operator()(Particle&) on the active type with
 * the particle that ought to be updated. This construction
 * allows to nest message variants to allow for sub-
 * categories. Those are mostly used here to differentiate
 * the updates for the substructs of Particle.
 */
struct UpdateVisitor : public boost::static_visitor<void> {
  explicit UpdateVisitor(int id) : id(id) {}

  const int id;

  /* Recurse into sub-variants */
  template <class... Message>
  void operator()(const boost::variant<Message...> &msg) const {
    boost::apply_visitor(*this, msg);
  }
  /* Plain messages are just called. */
  template <typename Message> void operator()(const Message &msg) const {
    assert(local_particles[id]);
    msg(*local_particles[id]);
  }
};
} // namespace

/**
 * @brief Callback for \ref mpi_send_update_message.
 */
void mpi_update_particle_slave(int node, int id) {
  if (node == comm_cart.rank()) {
    UpdateMessage msg{};
    comm_cart.recv(0, SOME_TAG, msg);
    boost::apply_visitor(UpdateVisitor{id}, msg);
  }

  on_particle_change();
}

/**
 * @brief Send a particle update message.
 *
 * This sends the message to the node that is responsible for the particle,
 * where
 * @p msg is called with the particle as argument. The message then performs the
 * change to the particle that is encoded in it. The mechanism to call a functor
 * based on the active type of a variant is called visitation. Here we can use
 * @c UpdateVisitor with a single templated call operator on the visitor because
 * all message (eventually) provide the same interface. Overall this is
 * logically equivalent to nested switch statements over the message types,
 * where the case statements play the role of the call of the messages in our
 * case. A general introduction can be found in the documentation of
 * boost::variant.
 *
 * @param id Id of the particle to update
 * @param msg The message
 */
void mpi_send_update_message(int id, const UpdateMessage &msg) {
  auto const pnode = get_particle_node(id);

  mpi_call(mpi_update_particle_slave, pnode, id);

  /* If the particle is remote, send the
   * message to the target, otherwise we
   * can just apply the update directly. */
  if (pnode != comm_cart.rank()) {
    comm_cart.send(pnode, SOME_TAG, msg);
  } else {
    boost::apply_visitor(UpdateVisitor{id}, msg);
  }

  on_particle_change();
}

template <typename S, S Particle::*s, typename T, T S::*m>
void mpi_update_particle(int id, const T &value) {
  using MessageType = message_type_t<S, s>;
  MessageType msg = UpdateParticle<S, s, T, m>{value};
  mpi_send_update_message(id, msg);
}

template <typename T, T ParticleProperties::*m>
void mpi_update_particle_property(int id, const T &value) {
  mpi_update_particle<ParticleProperties, &Particle::p, T, m>(id, value);
}

/************************************************
 * variables
 ************************************************/
bool type_list_enable;
std::unordered_map<int, std::unordered_set<int>> particle_type_map{};
void remove_id_from_map(int part_id, int type);
void add_id_to_type_map(int part_id, int type);

int max_seen_particle = -1;
int n_part = 0;
bool swimming_particles_exist = false;
/**
 * @brief id -> rank
 */
std::unordered_map<int, int> particle_node;

int max_local_particles = 0;
Particle **local_particles = nullptr;

/************************************************
 * local functions
 ************************************************/

int try_delete_bond(Particle *part, const int *bond);

void try_delete_exclusion(Particle *part, int part2);

void try_add_exclusion(Particle *part, int part2);

void auto_exclusion(int distance);

/************************************************
 * particle initialization functions
 ************************************************/

void free_particle(Particle *part) { part->~Particle(); }

void mpi_who_has_slave(int, int) {
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
  static auto *sizes = new int[n_nodes];
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
    throw std::runtime_error("Particle node for id " + std::to_string(id) +
                             " not found!");
  }
  return needle->second;
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
  assert(size >= 0);
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

void append_unindexed_particle(ParticleList *l, Particle &&part) {
  realloc_particlelist(l, ++l->n);
  new (&(l->part[l->n - 1])) Particle(std::move(part));
}

Particle *append_indexed_particle(ParticleList *l, Particle &&part) {
  auto const re = realloc_particlelist(l, ++l->n);
  auto p = new (&(l->part[l->n - 1])) Particle(std::move(part));

  assert(p->p.identity <= max_seen_particle);

  if (re)
    update_local_particles(l);
  else
    local_particles[p->p.identity] = p;
  return p;
}

Particle *move_unindexed_particle(ParticleList *dl, ParticleList *sl, int i) {
  assert(sl->n > 0);
  assert(i < sl->n);

  realloc_particlelist(dl, ++dl->n);
  auto dst = &dl->part[dl->n - 1];
  auto src = &sl->part[i];
  auto end = &sl->part[sl->n - 1];

  new (dst) Particle(std::move(*src));
  if (src != end) {
    new (src) Particle(std::move(*end));
  }

  realloc_particlelist(sl, --sl->n);
  return dst;
}

Particle *move_indexed_particle(ParticleList *dl, ParticleList *sl, int i) {
  assert(sl->n > 0);
  assert(i < sl->n);
  int re = realloc_particlelist(dl, ++dl->n);
  Particle *dst = &dl->part[dl->n - 1];
  Particle *src = &sl->part[i];
  Particle *end = &sl->part[sl->n - 1];

  new (dst) Particle(std::move(*src));

  assert(dst->p.identity <= max_seen_particle);

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

/**
 * @brief Extract an indexed particle from a list.
 *
 * Removes a particle from a particle list and
 * from the particle index.
 *
 * @param i Index of particle to remove,
 *          needs to be valid.
 * @param sl List to remove the particle from,
 *           needs to be non-empty.
 * @return The extracted particle.
 */
Particle extract_indexed_particle(ParticleList *sl, int i) {
  assert(sl->n > 0);
  assert(i < sl->n);
  Particle *src = &sl->part[i];
  Particle *end = &sl->part[sl->n - 1];

  Particle p = std::move(*src);

  assert(p.p.identity <= max_seen_particle);
  local_particles[p.p.identity] = nullptr;

  if (src != end) {
    new (src) Particle(std::move(*end));
  }

  if (realloc_particlelist(sl, --sl->n)) {
    update_local_particles(sl);
  } else if (src != end) {
    local_particles[src->p.identity] = src;
  }
  return p;
}

namespace {
/* Limit cache to 100 MiB */
std::size_t const max_cache_size = (100ul * 1048576ul) / sizeof(Particle);
Utils::Cache<int, Particle> particle_fetch_cache(max_cache_size);
} // namespace

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
                   assert(local_particles[id]);
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
                             }
                             auto const pnode = get_particle_node(id);
                             return (pnode == this_node) ||
                                    particle_fetch_cache.has(id);
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
    pnode = cell_structure.position_to_node(Utils::Vector3d{p, p + 3});

    /* master node specific stuff */
    particle_node[part] = pnode;

    retcode = ES_PART_CREATED;

    mpi_place_new_particle(pnode, part, p);
  }

  return retcode;
}

void set_particle_v(int part, double *v) {
  mpi_update_particle<ParticleMomentum, &Particle::m, Utils::Vector3d,
                      &ParticleMomentum::v>(part, Utils::Vector3d(v, v + 3));
}

#ifdef ENGINE
void set_particle_swimming(int part, ParticleParametersSwimming swim) {
  mpi_send_update_message(part, UpdateSwim{swim});
}
#endif

void set_particle_f(int part, const Utils::Vector3d &F) {
  mpi_update_particle<ParticleForce, &Particle::f, Utils::Vector3d,
                      &ParticleForce::f>(part, F);
}

#if defined(MASS)
void set_particle_mass(int part, double mass) {
  mpi_update_particle_property<double, &ParticleProperties::mass>(part, mass);
}
#else
const constexpr double ParticleProperties::mass;
#endif

#ifdef ROTATIONAL_INERTIA
void set_particle_rotational_inertia(int part, double *rinertia) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::rinertia>(
      part, Utils::Vector3d(rinertia, rinertia + 3));
}
#else
constexpr Utils::Vector3d ParticleProperties::rinertia;
#endif
#ifdef ROTATION
void set_particle_rotation(int part, int rot) {
  mpi_update_particle_property<int, &ParticleProperties::rotation>(part, rot);
}
#endif
#ifdef ROTATION
void rotate_particle(int part, const Utils::Vector3d &axis, double angle) {
  mpi_send_update_message(part, UpdateOrientation{axis, angle});
}
#endif

#ifdef AFFINITY
void set_particle_affinity(int part, double *bond_site) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::bond_site>(
      part, Utils::Vector3d(bond_site, bond_site + 3));
}
#endif

#ifdef MEMBRANE_COLLISION
void set_particle_out_direction(int part, double *out_direction) {
  mpi_update_particle_property<Utils::Vector3d,
                               &ParticleProperties::out_direction>(
      part, Utils::Vector3d(out_direction, out_direction + 3));
}
#endif

#ifdef DIPOLES
void set_particle_dipm(int part, double dipm) {
  mpi_update_particle_property<double, &ParticleProperties::dipm>(part, dipm);
}

void set_particle_dip(int part, double const *const dip) {
  Utils::Vector4d quat;
  double dipm;
  std::tie(quat, dipm) =
      convert_dip_to_quat(Utils::Vector3d({dip[0], dip[1], dip[2]}));

  set_particle_dipm(part, dipm);
  set_particle_quat(part, quat.data());
}

#endif

#ifdef VIRTUAL_SITES
void set_particle_virtual(int part, int is_virtual) {
  mpi_update_particle_property<int, &ParticleProperties::is_virtual>(
      part, is_virtual);
}
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void set_particle_vs_quat(int part, double *vs_relative_quat) {
  auto vs_relative = get_particle_data(part).p.vs_relative;
  vs_relative.quat = Utils::Vector4d(vs_relative_quat, vs_relative_quat + 4);

  mpi_update_particle_property<
      ParticleProperties::VirtualSitesRelativeParameteres,
      &ParticleProperties::vs_relative>(part, vs_relative);
}

void set_particle_vs_relative(int part, int vs_relative_to, double vs_distance,
                              double *rel_ori) {
  ParticleProperties::VirtualSitesRelativeParameteres vs_relative;
  vs_relative.distance = vs_distance;
  vs_relative.to_particle_id = vs_relative_to;
  vs_relative.rel_orientation = {rel_ori, rel_ori + 4};

  mpi_update_particle_property<
      ParticleProperties::VirtualSitesRelativeParameteres,
      &ParticleProperties::vs_relative>(part, vs_relative);
}
#endif

void set_particle_q(int part, double q) {
#ifdef ELECTROSTATICS
  mpi_update_particle_property<double, &ParticleProperties::q>(part, q);
#endif
}

#ifndef ELECTROSTATICS
const constexpr double ParticleProperties::q;
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
void set_particle_mu_E(int part, double *mu_E) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::mu_E>(
      part, Utils::Vector3d(mu_E, mu_E + 3));
}

void get_particle_mu_E(int part, double (&mu_E)[3]) {
  auto const &p = get_particle_data(part);

  for (int i = 0; i < 3; i++) {
    mu_E[i] = p.p.mu_E[i];
  }
}
#endif

void set_particle_type(int p_id, int type) {
  make_particle_type_exist(type);

  if (type_list_enable) {
    // check if the particle exists already and the type is changed, then remove
    // it from the list which contains it
    auto const &cur_par = get_particle_data(p_id);
    int prev_type = cur_par.p.type;
    if (prev_type != type) {
      // particle existed before so delete it from the list
      remove_id_from_map(p_id, prev_type);
    }
    add_id_to_type_map(p_id, type);
  }

  mpi_update_particle_property<int, &ParticleProperties::type>(p_id, type);
}

void set_particle_mol_id(int part, int mid) {
  mpi_update_particle_property<int, &ParticleProperties::mol_id>(part, mid);
}

#ifdef ROTATION
void set_particle_quat(int part, double *quat) {
  mpi_update_particle<ParticlePosition, &Particle::r, Utils::Vector4d,
                      &ParticlePosition::quat>(part,
                                               Utils::Vector4d(quat, quat + 4));
}

void set_particle_omega_lab(int part, const Utils::Vector3d &omega_lab) {
  auto const &particle = get_particle_data(part);

  mpi_update_particle<ParticleMomentum, &Particle::m, Utils::Vector3d,
                      &ParticleMomentum::omega>(
      part, convert_vector_space_to_body(particle, omega_lab));
}

void set_particle_omega_body(int part, const Utils::Vector3d &omega) {
  mpi_update_particle<ParticleMomentum, &Particle::m, Utils::Vector3d,
                      &ParticleMomentum::omega>(part, omega);
}

void set_particle_torque_lab(int part, const Utils::Vector3d &torque_lab) {
  auto const &particle = get_particle_data(part);

  mpi_update_particle<ParticleForce, &Particle::f, Utils::Vector3d,
                      &ParticleForce::torque>(
      part, convert_vector_space_to_body(particle, torque_lab));
}
#endif

#ifdef LANGEVIN_PER_PARTICLE
void set_particle_temperature(int part, double T) {
  mpi_update_particle_property<double, &ParticleProperties::T>(part, T);
}

#ifndef PARTICLE_ANISOTROPY
void set_particle_gamma(int part, double gamma) {
  mpi_update_particle_property<double, &ParticleProperties::gamma>(part, gamma);
}
#else
void set_particle_gamma(int part, Utils::Vector3d gamma) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::gamma>(
      part, gamma);
}
#endif // PARTICLE_ANISOTROPY

#ifdef ROTATION
#ifndef PARTICLE_ANISOTROPY
void set_particle_gamma_rot(int part, double gamma_rot) {
  mpi_update_particle_property<double, &ParticleProperties::gamma_rot>(
      part, gamma_rot);
}
#else
void set_particle_gamma_rot(int part, Utils::Vector3d gamma_rot) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::gamma_rot>(
      part, gamma_rot);
}
#endif // PARTICLE_ANISOTROPY
#endif // ROTATION
#endif // LANGEVIN_PER_PARTICLE

#ifdef EXTERNAL_FORCES
#ifdef ROTATION
void set_particle_ext_torque(int part, const Utils::Vector3d &torque) {
  auto const flag = (!torque.empty()) ? PARTICLE_EXT_TORQUE : 0;
  if (flag) {
    mpi_update_particle_property<Utils::Vector3d,
                                 &ParticleProperties::ext_torque>(part, torque);
  }
  mpi_send_update_message(part, UpdatePropertyMessage(UpdateExternalFlag{
                                    PARTICLE_EXT_TORQUE, flag}));
}
#endif

void set_particle_ext_force(int part, const Utils::Vector3d &force) {
  auto const flag = (!force.empty()) ? PARTICLE_EXT_FORCE : 0;
  if (flag) {
    mpi_update_particle_property<Utils::Vector3d,
                                 &ParticleProperties::ext_force>(part, force);
  }
  mpi_send_update_message(part, UpdatePropertyMessage(UpdateExternalFlag{
                                    PARTICLE_EXT_FORCE, flag}));
}

void set_particle_fix(int part, int flag) {
  mpi_send_update_message(
      part, UpdatePropertyMessage(UpdateExternalFlag{COORDS_FIX_MASK, flag}));
}
#endif

void delete_particle_bond(int part, Utils::Span<const int> bond) {
  mpi_send_update_message(
      part, UpdateBondMessage{RemoveBond{{bond.begin(), bond.end()}}});
}

void delete_particle_bonds(int part) {
  mpi_send_update_message(part, UpdateBondMessage{RemoveBonds{}});
}

void add_particle_bond(int part, Utils::Span<const int> bond) {
  mpi_send_update_message(
      part, UpdateBondMessage{AddBond{{bond.begin(), bond.end()}}});
}

void remove_all_particles() {
  mpi_remove_particle(-1, -1);
  clear_particle_node();
}

int remove_particle(int p_id) {
  auto const &cur_par = get_particle_data(p_id);
  if (type_list_enable) {
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

namespace {
std::pair<Cell *, size_t> find_particle(Particle *p, Cell *c) {
  for (int i = 0; i < c->n; ++i) {
    if ((c->part + i) == p) {
      return {c, i};
    }
  }
  return {nullptr, 0};
}

std::pair<Cell *, size_t> find_particle(Particle *p, CellPList cells) {
  for (auto &c : cells) {
    auto res = find_particle(p, c);
    if (res.first) {
      return res;
    }
  }

  return {nullptr, 0};
}
} // namespace

void local_remove_particle(int part) {
  Particle *p = local_particles[part];
  assert(p);
  assert(not p->l.ghost);

  /* If the particles are sorted we can use the
   * cell system to find the cell containing the
   * particle. Otherwise we do a brute force search
   * of the cells. */
  Cell *cell = nullptr;
  size_t n = 0;
  if (Cells::RESORT_NONE == get_resort_particles()) {
    std::tie(cell, n) = find_particle(p, find_current_cell(*p));
  }

  if (not cell) {
    std::tie(cell, n) = find_particle(p, local_cells);
  }

  assert(cell && cell->part && (n < cell->n) && ((cell->part + n) == p));

  Particle p_destroy = extract_indexed_particle(cell, n);
}

Particle *local_place_particle(int part, const double p[3], int _new) {
  Particle *pt;

  Utils::Vector3i i{};
  Utils::Vector3d pp = {p[0], p[1], p[2]};
  fold_position(pp, i);

  if (_new) {
    /* allocate particle anew */
    auto cell = cell_structure.position_to_cell(pp);
    if (!cell) {
      fprintf(stderr,
              "%d: INTERNAL ERROR: particle %d at %f(%f) %f(%f) %f(%f) "
              "does not belong on this node\n",
              this_node, part, p[0], pp[0], p[1], pp[1], p[2], pp[2]);
      errexit();
    }
    auto rl = realloc_particlelist(cell, ++cell->n);
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

  pt->r.p = pp;
  pt->l.i = i;
#ifdef BOND_CONSTRAINT
  pt->r.p_old = pp;
#endif

  assert(local_particles[part] == pt);

  return pt;
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
      p.r.p *= scale;
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

void local_add_particle_bond(Particle &p, Utils::Span<const int> bond) {
  boost::copy(bond, std::back_inserter(p.bl));
}

int try_delete_bond(Particle *part, const int *bond) {
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
      fprintf(stderr,
              "%d: INTERNAL ERROR: bond information corrupt for "
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

  if (!el.empty()) {
    el.erase(std::remove(el.begin(), el.end(), part2), el.end());
  };
}
#endif

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
} // namespace

int change_exclusion(int part1, int part2, int _delete) {
  if (particle_exists(part1) && particle_exists(part2)) {
    mpi_send_exclusion(part1, part2, _delete);
    return ES_OK;
  }
  return ES_ERROR;
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
     restore it continuously. Of course this could be optimized by bundling the
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
  if (particle_type_map.find(type) != particle_type_map.end())
    particle_type_map.at(type).erase(part_id);
}

int get_random_p_id(int type) {
  if (particle_type_map.at(type).empty())
    throw std::runtime_error("No particles of given type could be found");
  int rand_index = i_random(particle_type_map.at(type).size());
  return *std::next(particle_type_map[type].begin(), rand_index);
}

void add_id_to_type_map(int part_id, int type) {
  if (particle_type_map.find(type) != particle_type_map.end())
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

void pointer_to_quat(Particle const *p, double const *&res) {
  res = p->r.quat.data();
}

#endif

void pointer_to_q(Particle const *p, double const *&res) { res = &(p->p.q); }

#ifdef VIRTUAL_SITES
void pointer_to_virtual(Particle const *p, int const *&res) {
  res = &(p->p.is_virtual);
}
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void pointer_to_vs_quat(Particle const *p, double const *&res) {
  res = (p->p.vs_relative.quat.data());
}

void pointer_to_vs_relative(Particle const *p, int const *&res1,
                            double const *&res2, double const *&res3) {
  res1 = &(p->p.vs_relative.to_particle_id);
  res2 = &(p->p.vs_relative.distance);
  res3 = p->p.vs_relative.rel_orientation.data();
}
#endif

#ifdef DIPOLES

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

void pointer_to_rotation(Particle const *p, int const *&res) {
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
void pointer_to_bond_site(Particle const *p, double const *&res) {
  res = p->p.bond_site.data();
}
#endif

#ifdef MEMBRANE_COLLISION
void pointer_to_out_direction(const Particle *p, const double *&res) {
  res = p->p.out_direction.data();
}
#endif

bool particle_exists(int part_id) {
  if (particle_node.empty())
    build_particle_node();
  return particle_node.count(part_id);
}

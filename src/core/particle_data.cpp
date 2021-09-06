/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
/** \file
 *  Particles and particle lists.
 *
 *  The corresponding header file is particle_data.hpp.
 */
#include "particle_data.hpp"

#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "grid.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "partCfg_global.hpp"
#include "rotation.hpp"

#include <string>
#include <utils/Cache.hpp>
#include <utils/constants.hpp>
#include <utils/keys.hpp>
#include <utils/mpi/gatherv.hpp>

#include <boost/algorithm/cxx11/copy_if.hpp>
#include <boost/mpi/collectives/gather.hpp>
#include <boost/mpi/collectives/scatter.hpp>
#include <boost/optional.hpp>
#include <boost/range/algorithm.hpp>
#include <boost/range/numeric.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace {
/**
 * @brief A generic particle update.
 *
 * Here the sub-struct structure of Particle is
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

using Prop = ParticleProperties;

// clang-format off
using UpdatePropertyMessage = boost::variant
        < UpdateProperty<int, &Prop::type>
        , UpdateProperty<int, &Prop::mol_id>
#ifdef MASS
        , UpdateProperty<double, &Prop::mass>
#endif
#ifdef ROTATIONAL_INERTIA
        , UpdateProperty<Utils::Vector3d, &Prop::rinertia>
#endif
#ifdef ROTATION
        , UpdateProperty<uint8_t, &Prop::rotation>
#endif
#ifdef ELECTROSTATICS
        , UpdateProperty<double, &Prop::q>
#endif
#ifdef LB_ELECTROHYDRODYNAMICS
        , UpdateProperty<Utils::Vector3d, &Prop::mu_E>
#endif
#ifdef ENGINE
        , UpdateProperty<ParticleParametersSwimming, &Prop::swim>
#endif
#ifdef DIPOLES
        , UpdateProperty<double, &Prop::dipm>
#endif
#ifdef VIRTUAL_SITES
        , UpdateProperty<bool, &Prop::is_virtual>
#ifdef VIRTUAL_SITES_RELATIVE
        , UpdateProperty<ParticleProperties::VirtualSitesRelativeParameters,
                         &Prop::vs_relative>
#endif
#endif
#ifdef THERMOSTAT_PER_PARTICLE
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
#endif // PARTICLE_ANISOTROPY
#endif // ROTATION
#endif // THERMOSTAT_PER_PARTICLE
#ifdef EXTERNAL_FORCES
        , UpdateProperty<uint8_t, &Prop::ext_flag>
        , UpdateProperty<Utils::Vector3d, &Prop::ext_force>
#ifdef ROTATION
        , UpdateProperty<Utils::Vector3d, &Prop::ext_torque>
#endif
#endif
        >;

using UpdatePositionMessage = boost::variant
        < UpdatePosition<Utils::Vector3d, &ParticlePosition::p>
#ifdef ROTATION
        , UpdatePosition<Utils::Quaternion<double>, &ParticlePosition::quat>
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
      assert(not bond.empty());
      auto const view = BondView(bond.front(), {bond.data() + 1, bond.size() - 1});
      auto it = boost::find(p.bonds(), view);

      if(it != p.bonds().end()) {
       p.bonds().erase(it);
      }
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
      p.bonds().clear();
    }

    template<class Archive>
    void serialize(Archive &ar, long int) {
    }
};

struct AddBond {
    std::vector<int> bond;

    void operator()(Particle &p) const {
      auto const view = BondView(bond.at(0), {bond.data() + 1, bond.size() - 1});

      p.bonds().insert(view);
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
 * fulfill the type requirements, namely:
 * They either have an integer member id indicating the particle
 * that should be updated with an <tt>operator()(const Particle&)</tt>
 * that is called with the particle, or a tree of
 * variants with leaves that have such an <tt>operator()</tt> member.
 */
using UpdateMessage = boost::variant
        < UpdatePropertyMessage
        , UpdatePositionMessage
        , UpdateMomentumMessage
        , UpdateForceMessage
        , UpdateBondMessage
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
 * if it is a variant type, or otherwise calls
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
    assert(cell_structure.get_local_particle(id));
    msg(*cell_structure.get_local_particle(id));
  }
};
} // namespace

void mpi_send_update_message_local(int node, int id) {
  if (node == comm_cart.rank()) {
    UpdateMessage msg{};
    comm_cart.recv(0, SOME_TAG, msg);
    boost::apply_visitor(UpdateVisitor{id}, msg);
  }

  on_particle_change();
}

REGISTER_CALLBACK(mpi_send_update_message_local)

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

  mpi_call(mpi_send_update_message_local, pnode, id);

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

/**
 * @brief id -> rank
 */
std::unordered_map<int, int> particle_node;

void delete_exclusion(Particle *part, int part2);

void add_exclusion(Particle *part, int part2);

void auto_exclusion(int distance);

void mpi_who_has_local(int, int) {
  static std::vector<int> sendbuf;

  auto local_particles = cell_structure.local_particles();
  auto const n_part = static_cast<int>(local_particles.size());
  boost::mpi::gather(comm_cart, n_part, 0);

  if (n_part == 0)
    return;

  sendbuf.resize(n_part);

  std::transform(local_particles.begin(), local_particles.end(),
                 sendbuf.begin(),
                 [](Particle const &p) { return p.p.identity; });

  MPI_Send(sendbuf.data(), n_part, MPI_INT, 0, SOME_TAG, comm_cart);
}

REGISTER_CALLBACK(mpi_who_has_local)

void mpi_who_has() {
  mpi_call(mpi_who_has_local, -1, 0);

  auto local_particles = cell_structure.local_particles();

  static std::vector<int> n_parts;
  boost::mpi::gather(comm_cart, static_cast<int>(local_particles.size()),
                     n_parts, 0);

  static std::vector<int> pdata;

  /* then fetch particle locations */
  for (int pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      for (auto const &p : local_particles)
        particle_node[p.p.identity] = this_node;

    } else if (n_parts[pnode] > 0) {
      pdata.resize(n_parts[pnode]);
      MPI_Recv(pdata.data(), n_parts[pnode], MPI_INT, pnode, SOME_TAG,
               comm_cart, MPI_STATUS_IGNORE);
      for (int i = 0; i < n_parts[pnode]; i++)
        particle_node[pdata[i]] = pnode;
    }
  }
}

/**
 * @brief Rebuild the particle index.
 */
void build_particle_node() { mpi_who_has(); }

/**
 *  @brief Get the mpi rank which owns the particle with id.
 */
int get_particle_node(int id) {
  if (id < 0)
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

namespace {
/* Limit cache to 100 MiB */
std::size_t const max_cache_size = (100ul * 1048576ul) / sizeof(Particle);
Utils::Cache<int, Particle> particle_fetch_cache(max_cache_size);
} // namespace

void invalidate_fetch_cache() { particle_fetch_cache.invalidate(); }
std::size_t fetch_cache_max_size() { return particle_fetch_cache.max_size(); }

boost::optional<const Particle &> get_particle_data_local(int id) {
  auto p = cell_structure.get_local_particle(id);

  if (p and (not p->l.ghost)) {
    return *p;
  }

  return {};
}

REGISTER_CALLBACK_ONE_RANK(get_particle_data_local)

const Particle &get_particle_data(int part) {
  auto const pnode = get_particle_node(part);

  if (pnode == this_node) {
    assert(cell_structure.get_local_particle(part));
    return *cell_structure.get_local_particle(part);
  }

  /* Query the cache */
  auto const p_ptr = particle_fetch_cache.get(part);
  if (p_ptr) {
    return *p_ptr;
  }

  /* Cache miss, fetch the particle,
   * put it into the cache and return a pointer into the cache. */
  auto const cache_ptr = particle_fetch_cache.put(
      part, Communication::mpiCallbacks().call(Communication::Result::one_rank,
                                               get_particle_data_local, part));
  return *cache_ptr;
}

void mpi_get_particles_local(int, int) {
  std::vector<int> ids;
  boost::mpi::scatter(comm_cart, ids, 0);

  std::vector<Particle> parts(ids.size());
  std::transform(ids.begin(), ids.end(), parts.begin(), [](int id) {
    assert(cell_structure.get_local_particle(id));
    return *cell_structure.get_local_particle(id);
  });

  Utils::Mpi::gatherv(comm_cart, parts.data(), parts.size(), 0);
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
std::vector<Particle> mpi_get_particles(Utils::Span<const int> ids) {
  mpi_call(mpi_get_particles_local, 0, 0);
  /* Return value */
  std::vector<Particle> parts(ids.size());

  /* Group ids per node */
  static std::vector<std::vector<int>> node_ids(comm_cart.size());
  for (auto &per_node : node_ids) {
    per_node.clear();
  }

  for (auto const &id : ids) {
    auto const pnode = get_particle_node(id);

    if (pnode != this_node)
      node_ids[pnode].push_back(id);
  }

  /* Distributed ids to the nodes */
  {
    static std::vector<int> ignore;
    boost::mpi::scatter(comm_cart, node_ids, ignore, 0);
  }

  /* Copy local particles */
  std::transform(node_ids[this_node].cbegin(), node_ids[this_node].cend(),
                 parts.begin(), [](int id) {
                   assert(cell_structure.get_local_particle(id));
                   return *cell_structure.get_local_particle(id);
                 });

  static std::vector<int> node_sizes(comm_cart.size());
  std::transform(
      node_ids.cbegin(), node_ids.cend(), node_sizes.begin(),
      [](std::vector<int> const &ids) { return static_cast<int>(ids.size()); });

  Utils::Mpi::gatherv(comm_cart, parts.data(), parts.size(), parts.data(),
                      node_sizes.data(), 0);

  return parts;
}

void prefetch_particle_data(Utils::Span<const int> in_ids) {
  /* Nothing to do on a single node. */
  // NOLINTNEXTLINE(clang-analyzer-core.NonNullParamChecker)
  if (comm_cart.size() == 1)
    return;

  static std::vector<int> ids;
  ids.clear();

  boost::algorithm::copy_if(in_ids, std::back_inserter(ids), [](int id) {
    return (get_particle_node(id) != this_node) && particle_fetch_cache.has(id);
  });

  /* Don't prefetch more particles than fit the cache. */
  if (ids.size() > particle_fetch_cache.max_size())
    ids.resize(particle_fetch_cache.max_size());

  /* Fetch the particles... */
  for (auto &p : mpi_get_particles(ids)) {
    auto id = p.identity();
    particle_fetch_cache.put(id, std::move(p));
  }
}

/** Move a particle to a new position. If it does not exist, it is created.
 *  The position must be on the local node!
 *
 *  @param id    the identity of the particle to move
 *  @param pos   its new position
 *  @param _new  if true, the particle is allocated, else has to exists already
 *
 *  @return Pointer to the particle.
 */
Particle *local_place_particle(int id, const Utils::Vector3d &pos, int _new) {
  auto pp = Utils::Vector3d{pos[0], pos[1], pos[2]};
  auto i = Utils::Vector3i{};
  fold_position(pp, i, box_geo);

  if (_new) {
    Particle new_part;
    new_part.p.identity = id;
    new_part.r.p = pp;
    new_part.l.i = i;

    return cell_structure.add_local_particle(std::move(new_part));
  }

  auto pt = cell_structure.get_local_particle(id);
  pt->r.p = pp;
  pt->l.i = i;

  return pt;
}

boost::optional<int> mpi_place_new_particle_local(int p_id,
                                                  Utils::Vector3d const &pos) {
  auto p = local_place_particle(p_id, pos, 1);
  on_particle_change();
  if (p) {
    return comm_cart.rank();
  }
  return {};
}

REGISTER_CALLBACK_ONE_RANK(mpi_place_new_particle_local)

/** Create particle at a position on a node.
 *  Also calls \ref on_particle_change.
 *  \param p_id  the particle to create.
 *  \param pos   the particles position.
 */
int mpi_place_new_particle(int p_id, const Utils::Vector3d &pos) {
  return mpi_call(Communication::Result::one_rank, mpi_place_new_particle_local,
                  p_id, pos);
}

void mpi_place_particle_local(int pnode, int p_id) {
  if (pnode == this_node) {
    Utils::Vector3d pos;
    comm_cart.recv(0, SOME_TAG, pos);
    local_place_particle(p_id, pos, 0);
  }

  cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_place_particle_local)

/** Move particle to a position on a node.
 *  Also calls \ref on_particle_change.
 *  \param node  the node to attach it to.
 *  \param p_id  the particle to move.
 *  \param pos   the particles position.
 */
void mpi_place_particle(int node, int p_id, const Utils::Vector3d &pos) {
  mpi_call(mpi_place_particle_local, node, p_id);

  if (node == this_node)
    local_place_particle(p_id, pos, 0);
  else {
    comm_cart.send(node, SOME_TAG, pos);
  }

  cell_structure.set_resort_particles(Cells::RESORT_GLOBAL);
  on_particle_change();
}

int place_particle(int p_id, Utils::Vector3d const &pos) {
  if (particle_exists(p_id)) {
    mpi_place_particle(get_particle_node(p_id), p_id, pos);

    return ES_PART_OK;
  }
  particle_node[p_id] = mpi_place_new_particle(p_id, pos);

  return ES_PART_CREATED;
}

void set_particle_v(int part, Utils::Vector3d const &v) {
  mpi_update_particle<ParticleMomentum, &Particle::m, Utils::Vector3d,
                      &ParticleMomentum::v>(part, v);
}

#ifdef ENGINE
void set_particle_swimming(int part, ParticleParametersSwimming swim) {
  mpi_update_particle_property<ParticleParametersSwimming,
                               &ParticleProperties::swim>(part, swim);
}
#endif

void set_particle_f(int part, const Utils::Vector3d &f) {
  mpi_update_particle<ParticleForce, &Particle::f, Utils::Vector3d,
                      &ParticleForce::f>(part, f);
}

#if defined(MASS)
void set_particle_mass(int part, double mass) {
  mpi_update_particle_property<double, &ParticleProperties::mass>(part, mass);
}
#else
const constexpr double ParticleProperties::mass;
#endif

#ifdef ROTATIONAL_INERTIA
void set_particle_rotational_inertia(int part,
                                     Utils::Vector3d const &rinertia) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::rinertia>(
      part, rinertia);
}
#else
constexpr Utils::Vector3d ParticleProperties::rinertia;
#endif

#ifdef ROTATION
void set_particle_rotation(int part, int rot) {
  mpi_update_particle_property<uint8_t, &ParticleProperties::rotation>(part,
                                                                       rot);
}

void rotate_particle(int part, const Utils::Vector3d &axis, double angle) {
  mpi_send_update_message(part, UpdateOrientation{axis, angle});
}
#endif

#ifdef DIPOLES
void set_particle_dipm(int part, double dipm) {
  mpi_update_particle_property<double, &ParticleProperties::dipm>(part, dipm);
}

void set_particle_dip(int part, Utils::Vector3d const &dip) {
  Utils::Quaternion<double> quat;
  double dipm;
  std::tie(quat, dipm) = convert_dip_to_quat(dip);

  set_particle_dipm(part, dipm);
  set_particle_quat(part, quat.data());
}
#endif

#ifdef VIRTUAL_SITES
void set_particle_virtual(int part, bool is_virtual) {
  mpi_update_particle_property<bool, &ParticleProperties::is_virtual>(
      part, is_virtual);
}
#endif

#ifdef VIRTUAL_SITES_RELATIVE
void set_particle_vs_quat(int part,
                          Utils::Quaternion<double> const &vs_relative_quat) {
  auto vs_relative = get_particle_data(part).p.vs_relative;
  vs_relative.quat = vs_relative_quat;

  mpi_update_particle_property<
      ParticleProperties::VirtualSitesRelativeParameters,
      &ParticleProperties::vs_relative>(part, vs_relative);
}

void set_particle_vs_relative(int part, int vs_relative_to, double vs_distance,
                              Utils::Quaternion<double> const &rel_ori) {
  ParticleProperties::VirtualSitesRelativeParameters vs_relative;
  vs_relative.distance = vs_distance;
  vs_relative.to_particle_id = vs_relative_to;
  vs_relative.rel_orientation = rel_ori;

  mpi_update_particle_property<
      ParticleProperties::VirtualSitesRelativeParameters,
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
void set_particle_mu_E(int part, Utils::Vector3d const &mu_E) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::mu_E>(
      part, mu_E);
}

Utils::Vector3d get_particle_mu_E(int part) {
  auto const &p = get_particle_data(part);
  return p.p.mu_E;
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
void set_particle_quat(int part, double *const quat) {
  mpi_update_particle<ParticlePosition, &Particle::r, Utils::Quaternion<double>,
                      &ParticlePosition::quat>(
      part, Utils::Quaternion<double>{quat[0], quat[1], quat[2], quat[3]});
}

void set_particle_director(int part, const Utils::Vector3d &director) {
  Utils::Quaternion<double> quat =
      convert_director_to_quaternion(director.normalized());
  set_particle_quat(part, quat.data());
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
#endif // ROTATION

#ifdef THERMOSTAT_PER_PARTICLE
#ifndef PARTICLE_ANISOTROPY
void set_particle_gamma(int part, double gamma) {
  mpi_update_particle_property<double, &ParticleProperties::gamma>(part, gamma);
}
#else
void set_particle_gamma(int part, Utils::Vector3d const &gamma) {
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
void set_particle_gamma_rot(int part, Utils::Vector3d const &gamma_rot) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::gamma_rot>(
      part, gamma_rot);
}
#endif // PARTICLE_ANISOTROPY
#endif // ROTATION
#endif // THERMOSTAT_PER_PARTICLE

#ifdef EXTERNAL_FORCES
#ifdef ROTATION
void set_particle_ext_torque(int part, const Utils::Vector3d &torque) {
  mpi_update_particle_property<Utils::Vector3d,
                               &ParticleProperties::ext_torque>(part, torque);
}
#endif // ROTATION

void set_particle_ext_force(int part, const Utils::Vector3d &force) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::ext_force>(
      part, force);
}

void set_particle_fix(int part, uint8_t flag) {
  mpi_update_particle_property<uint8_t, &ParticleProperties::ext_flag>(part,
                                                                       flag);
}
#endif // EXTERNAL_FORCES

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

const std::vector<BondView> &get_particle_bonds(int part) {
  static std::vector<BondView> ret;
  ret.clear();

  boost::copy(get_particle_data(part).bonds(), std::back_inserter(ret));

  return ret;
}

void mpi_remove_particle_local(int, int part) {
  if (part != -1) {
    cell_structure.remove_particle(part);
  } else {
    cell_structure.remove_all_particles();
  }
  on_particle_change();
}

REGISTER_CALLBACK(mpi_remove_particle_local)

/** Remove a particle.
 *  Also calls \ref on_particle_change.
 *  \param p_id  the particle to remove, use -1 to remove all particles.
 */
void mpi_remove_particle(int, int p_id) {
  mpi_call_all(mpi_remove_particle_local, -1, p_id);
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

  particle_node[p_id] = -1;
  mpi_remove_particle(-1, p_id);

  particle_node.erase(p_id);

  return ES_OK;
}

/** Locally rescale all particles on current node.
 *  @param dir   direction to scale (0/1/2 = x/y/z, 3 = x+y+z isotropically)
 *  @param scale factor by which to rescale (>1: stretch, <1: contract)
 */
void local_rescale_particles(int dir, double scale) {
  for (auto &p : cell_structure.local_particles()) {
    if (dir < 3)
      p.r.p[dir] *= scale;
    else {
      p.r.p *= scale;
    }
  }
}

void mpi_rescale_particles_local(int, int dir) {
  double scale = 0.0;
  MPI_Recv(&scale, 1, MPI_DOUBLE, 0, SOME_TAG, comm_cart, MPI_STATUS_IGNORE);
  local_rescale_particles(dir, scale);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_rescale_particles_local)

void mpi_rescale_particles(int dir, double scale) {
  mpi_call(mpi_rescale_particles_local, -1, dir);
  for (int pnode = 0; pnode < n_nodes; pnode++) {
    if (pnode == this_node) {
      local_rescale_particles(dir, scale);
    } else {
      MPI_Send(&scale, 1, MPI_DOUBLE, pnode, SOME_TAG, comm_cart);
    }
  }
  on_particle_change();
}

#ifdef EXCLUSIONS
/** Locally add an exclusion to a particle.
 *  @param part1 the identity of the first exclusion partner
 *  @param part2 the identity of the second exclusion partner
 *  @param _delete if true, delete the exclusion instead of add
 */
void local_change_exclusion(int part1, int part2, int _delete) {
  /* part1, if here */
  auto part = cell_structure.get_local_particle(part1);
  if (part) {
    if (_delete)
      delete_exclusion(part, part2);
    else
      add_exclusion(part, part2);
  }

  /* part2, if here */
  part = cell_structure.get_local_particle(part2);
  if (part) {
    if (_delete)
      delete_exclusion(part, part1);
    else
      add_exclusion(part, part1);
  }
}

namespace {
/* keep a unique list for particle i. Particle j is only added if it is not i
   and not already in the list. */
void add_partner(std::vector<int> &il, int i, int j, int distance) {
  if (j == i)
    return;
  for (int k = 0; k < il.size(); k += 2)
    if (il[k] == j)
      return;

  il.push_back(j);
  il.push_back(distance);
}
} // namespace

void mpi_send_exclusion_local(int part1, int part2, int _delete) {
  local_change_exclusion(part1, part2, _delete);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_send_exclusion_local)

/** Send exclusions.
 *  Also calls \ref on_particle_change.
 *  \param part1    identity of first particle of the exclusion.
 *  \param part2    identity of second particle of the exclusion.
 *  \param _delete  if true, do not add the exclusion, rather delete it if found
 */
void mpi_send_exclusion(int part1, int part2, int _delete) {
  mpi_call_all(mpi_send_exclusion_local, part1, part2, _delete);
}

int change_exclusion(int part1, int part2, int _delete) {
  if (particle_exists(part1) && particle_exists(part2)) {
    mpi_send_exclusion(part1, part2, _delete);
    return ES_OK;
  }
  return ES_ERROR;
}

void auto_exclusions(int distance) {
  /* partners is a list containing the currently found excluded particles for
     each particle, and their distance, as an interleaved list */
  std::unordered_map<int, std::vector<int>> partners;

  /* determine initial connectivity */
  for (auto const &part1 : partCfg()) {
    auto const p1 = part1.p.identity;
    for (auto const bond : part1.bonds()) {
      if ((bond.partner_ids().size() == 1) and (bond.partner_ids()[0] != p1)) {
        auto const p2 = bond.partner_ids()[0];
        add_partner(partners[p1], p1, p2, 1);
        add_partner(partners[p2], p2, p1, 1);
      }
    }
  }

  /* calculate transient connectivity. For each of the current neighbors,
     also exclude their close enough neighbors.
  */
  for (int count = 1; count < distance; count++) {
    for (auto const &p : partCfg()) {
      auto const p1 = p.identity();
      for (int i = 0; i < partners[p1].size(); i += 2) {
        auto const p2 = partners[p1][i];
        auto const dist1 = partners[p1][i + 1];
        if (dist1 > distance)
          continue;
        /* loop over all partners of the partner */
        for (int j = 0; j < partners[p2].size(); j += 2) {
          auto const p3 = partners[p2][j];
          auto const dist2 = dist1 + partners[p2][j + 1];
          if (dist2 > distance)
            continue;
          add_partner(partners[p1], p1, p3, dist2);
          add_partner(partners[p3], p3, p1, dist2);
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
  for (auto &kv : partners) {
    auto const id = kv.first;
    auto const partner_list = kv.second;
    for (int j : partner_list)
      if (id < j)
        change_exclusion(id, j, 0);
  }
}
#endif // EXCLUSIONS

void init_type_map(int type) {
  type_list_enable = true;
  if (type < 0)
    throw std::runtime_error("Types may not be negative");

  auto &map_for_type = particle_type_map[type];
  map_for_type.clear();
  for (auto const &p : partCfg()) {
    if (p.p.type == type)
      map_for_type.insert(p.p.identity);
  }
}

void remove_id_from_map(int part_id, int type) {
  auto it = particle_type_map.find(type);
  if (it != particle_type_map.end())
    it->second.erase(part_id);
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

void add_id_to_type_map(int part_id, int type) {
  auto it = particle_type_map.find(type);
  if (it != particle_type_map.end())
    it->second.insert(part_id);
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
void pointer_to_virtual(Particle const *p, bool const *&res) {
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
void pointer_to_ext_force(Particle const *p, double const *&res2) {
  res2 = p->p.ext_force.data();
}
#ifdef ROTATION
void pointer_to_ext_torque(Particle const *p, double const *&res2) {
  res2 = p->p.ext_torque.data();
}
#endif
void pointer_to_fix(Particle const *p, const uint8_t *&res) {
  res = &(p->p.ext_flag);
}
#endif // EXTERNAL_FORCES

#ifdef THERMOSTAT_PER_PARTICLE
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
#endif // PARTICLE_ANISOTROPY
}
#endif // ROTATION

#endif // THERMOSTAT_PER_PARTICLE

#ifdef ENGINE
void pointer_to_swimming(Particle const *p,
                         ParticleParametersSwimming const *&swim) {
  swim = &(p->p.swim);
}
#endif

#ifdef ROTATIONAL_INERTIA
void pointer_to_rotational_inertia(Particle const *p, double const *&res) {
  res = p->p.rinertia.data();
}
#endif

bool particle_exists(int part_id) {
  if (particle_node.empty())
    build_particle_node();
  return particle_node.count(part_id);
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

  return boost::accumulate(particle_node, -1,
                           [](int max, const std::pair<int, int> &kv) {
                             return std::max(max, kv.first);
                           });
}

int get_n_part() {
  if (particle_node.empty())
    build_particle_node();

  return static_cast<int>(particle_node.size());
}

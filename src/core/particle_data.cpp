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

#include "particle_data.hpp"

#include "Particle.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "event.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "particle_node.hpp"

#include "config/config.hpp"

#include <utils/Vector.hpp>

#include <boost/serialization/variant.hpp>
#include <boost/variant.hpp>

#include <cassert>

constexpr auto some_tag = 42;

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

template <typename T, T ParticleLocal::*m>
using UpdateLocalProperty = UpdateParticle<ParticleLocal, &Particle::l, T, m>;
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
#ifdef ELECTROSTATICS
        , UpdateProperty<double, &Prop::q>
#endif
#ifdef EXTERNAL_FORCES
        , UpdateProperty<Utils::Vector3d, &Prop::ext_force>
#endif
        >;

using UpdatePositionMessage = boost::variant
        < UpdatePosition<Utils::Vector3d, &ParticlePosition::p> >;

using UpdateMomentumMessage = boost::variant
      < UpdateMomentum<Utils::Vector3d, &ParticleMomentum::v> >;

using UpdateForceMessage = boost::variant
      < UpdateForce<Utils::Vector3d, &ParticleForce::f> >;

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

static void mpi_send_update_message_local(int node, int id) {
  if (node == comm_cart.rank()) {
    UpdateMessage msg{};
    comm_cart.recv(0, some_tag, msg);
    boost::apply_visitor(UpdateVisitor(id), msg);
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
static void mpi_send_update_message(int id, const UpdateMessage &msg) {
  auto const pnode = get_particle_node(id);

  mpi_call(mpi_send_update_message_local, pnode, id);

  /* If the particle is remote, send the
   * message to the target, otherwise we
   * can just apply the update directly. */
  if (pnode != comm_cart.rank()) {
    comm_cart.send(pnode, some_tag, msg);
  } else {
    boost::apply_visitor(UpdateVisitor(id), msg);
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

void set_particle_v(int part, Utils::Vector3d const &v) {
  mpi_update_particle<ParticleMomentum, &Particle::m, Utils::Vector3d,
                      &ParticleMomentum::v>(part, v);
}

#ifdef ELECTROSTATICS
void set_particle_q(int part, double q) {
  mpi_update_particle_property<double, &ParticleProperties::q>(part, q);
}
#else
const constexpr double ParticleProperties::q;
#endif

void set_particle_type(int p_id, int type) {
  make_particle_type_exist(type);
  on_particle_type_change(p_id, type);
  mpi_update_particle_property<int, &ParticleProperties::type>(p_id, type);
}

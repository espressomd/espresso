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
#include "exclusions.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"
#include "partCfg_global.hpp"
#include "particle_node.hpp"
#include "rotation.hpp"

#include "config/config.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/quaternion.hpp>

#include <boost/range/algorithm/copy.hpp>
#include <boost/range/algorithm/find.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/serialization/vector.hpp>
#include <boost/variant.hpp>

#include <iterator>
#include <stdexcept>
#include <unordered_map>
#include <vector>

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

using UpdateLocalPropertyMessage = boost::variant
        < UpdateLocalProperty<double, &ParticleLocal::lees_edwards_offset>>;

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
// clang-format on

/**
 * @brief Delete specific bond.
 */
struct RemoveBond {
  std::vector<int> bond;

  void operator()(Particle &p) const {
    assert(not bond.empty());
    auto const view =
        BondView(bond.front(), {bond.data() + 1, bond.size() - 1});
    auto it = boost::find(p.bonds(), view);

    if (it != p.bonds().end()) {
      p.bonds().erase(it);
    }
  }

  template <class Archive> void serialize(Archive &ar, long int) { ar &bond; }
};

/**
 * @brief Delete pair bonds to a specific partner
 */
struct RemovePairBondsTo {
  int other_pid;

  void operator()(Particle &p) const {
    using Bond = std::vector<int>;
    std::vector<Bond> to_delete;
    for (auto b : p.bonds()) {
      if (b.partner_ids().size() == 1 and b.partner_ids()[0] == other_pid)
        to_delete.push_back(Bond{b.bond_id(), other_pid});
    }
    for (auto b : to_delete) {
      RemoveBond{b}(p);
    }
  }
  template <class Archive> void serialize(Archive &ar, long int) {
    ar &other_pid;
  }
};

/**
 * @brief Delete all bonds.
 */
struct RemoveBonds {
  void operator()(Particle &p) const { p.bonds().clear(); }

  template <class Archive> void serialize(Archive &, long int) {}
};

struct AddBond {
  std::vector<int> bond;

  void operator()(Particle &p) const {
    auto const view = BondView(bond.at(0), {bond.data() + 1, bond.size() - 1});

    p.bonds().insert(view);
  }

  template <class Archive> void serialize(Archive &ar, long int) { ar &bond; }
};

// clang-format off
using UpdateBondMessage = boost::variant
        < RemoveBond
        , RemoveBonds
        , AddBond
        >;
// clang-format on

#ifdef ROTATION
struct UpdateOrientation {
  Utils::Vector3d axis;
  double angle;

  void operator()(Particle &p) const { local_rotate_particle(p, axis, angle); }

  template <class Archive> void serialize(Archive &ar, long int) {
    ar &axis &angle;
  }
};
#endif

// clang-format off
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
        < UpdateLocalPropertyMessage
        , UpdatePropertyMessage
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

template <> struct message_type<ParticleLocal, &Particle::l> {
  using type = UpdateLocalPropertyMessage;
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

void local_remove_bond(Particle &p, std::vector<int> const &bond) {
  RemoveBond{bond}(p);
}

void local_remove_pair_bonds_to(Particle &p, int other_pid) {
  RemovePairBondsTo{other_pid}(p);
}

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

void set_particle_lees_edwards_offset(int part, const double v) {
  mpi_update_particle<ParticleLocal, &Particle::l, double,
                      &ParticleLocal::lees_edwards_offset>(part, v);
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

#ifdef MASS
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
void set_particle_rotation(int part, Utils::Vector3i const &flag) {
  auto rot_flag = static_cast<uint8_t>(0u);
  if (flag[0])
    rot_flag |= static_cast<uint8_t>(1u);
  if (flag[1])
    rot_flag |= static_cast<uint8_t>(2u);
  if (flag[2])
    rot_flag |= static_cast<uint8_t>(4u);
  mpi_update_particle_property<uint8_t, &ParticleProperties::rotation>(
      part, rot_flag);
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
  auto const [quat, dipm] = convert_dip_to_quat(dip);
  set_particle_dipm(part, dipm);
  set_particle_quat(part, quat);
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
  auto vs_relative = get_particle_data(part).vs_relative();
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

#ifdef ELECTROSTATICS
void set_particle_q(int part, double q) {
  mpi_update_particle_property<double, &ParticleProperties::q>(part, q);
}
#else
const constexpr double ParticleProperties::q;
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
void set_particle_mu_E(int part, Utils::Vector3d const &mu_E) {
  mpi_update_particle_property<Utils::Vector3d, &ParticleProperties::mu_E>(
      part, mu_E);
}
#endif

void set_particle_type(int p_id, int type) {
  make_particle_type_exist(type);
  on_particle_type_change(p_id, type);
  mpi_update_particle_property<int, &ParticleProperties::type>(p_id, type);
}

void set_particle_mol_id(int part, int mid) {
  mpi_update_particle_property<int, &ParticleProperties::mol_id>(part, mid);
}

#ifdef ROTATION
void set_particle_quat(int part, Utils::Quaternion<double> const &quat) {
  mpi_update_particle<ParticlePosition, &Particle::r, Utils::Quaternion<double>,
                      &ParticlePosition::quat>(part, quat);
}

void set_particle_director(int part, const Utils::Vector3d &director) {
  Utils::Quaternion<double> quat =
      convert_director_to_quaternion(director.normalized());
  set_particle_quat(part, quat);
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

void set_particle_fix(int part, Utils::Vector3i const &flag) {
  auto ext_flag = static_cast<uint8_t>(0u);
  if (flag[0])
    ext_flag |= static_cast<uint8_t>(1u);
  if (flag[1])
    ext_flag |= static_cast<uint8_t>(2u);
  if (flag[2])
    ext_flag |= static_cast<uint8_t>(4u);
  mpi_update_particle_property<uint8_t, &ParticleProperties::ext_flag>(
      part, ext_flag);
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

void rescale_particles(int dir, double scale) {
  for (auto &p : cell_structure.local_particles()) {
    if (dir < 3)
      p.pos()[dir] *= scale;
    else {
      p.pos() *= scale;
    }
  }
  on_particle_change();
}

#ifdef EXCLUSIONS
/**
 * @brief Locally remove an exclusion to a particle.
 * @param part1 the identity of the first exclusion partner
 * @param part2 the identity of the second exclusion partner
 */
static void local_remove_exclusion(int part1, int part2) {
  auto *p1 = cell_structure.get_local_particle(part1);
  if (p1) {
    delete_exclusion(*p1, part2);
  }
  auto *p2 = cell_structure.get_local_particle(part2);
  if (p2) {
    delete_exclusion(*p2, part1);
  }
}

/**
 * @brief Locally add an exclusion to a particle.
 * @param part1 the identity of the first exclusion partner
 * @param part2 the identity of the second exclusion partner
 */
static void local_add_exclusion(int part1, int part2) {
  auto *p1 = cell_structure.get_local_particle(part1);
  if (p1) {
    add_exclusion(*p1, part2);
  }
  auto *p2 = cell_structure.get_local_particle(part2);
  if (p2) {
    add_exclusion(*p2, part1);
  }
}

/* keep a unique list for particle i. Particle j is only added if it is not i
   and not already in the list. */
static void add_partner(std::vector<int> &il, int i, int j, int distance) {
  if (j == i)
    return;
  for (int k = 0; k < il.size(); k += 2)
    if (il[k] == j)
      return;

  il.push_back(j);
  il.push_back(distance);
}

static void mpi_remove_exclusion_local(int part1, int part2) {
  local_remove_exclusion(part1, part2);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_remove_exclusion_local)

static void mpi_add_exclusion_local(int part1, int part2) {
  local_add_exclusion(part1, part2);
  on_particle_change();
}

REGISTER_CALLBACK(mpi_add_exclusion_local)

static void check_particle_exists(int p_id) {
  if (not particle_exists(p_id)) {
    throw std::runtime_error("Particle with id " + std::to_string(p_id) +
                             " not found");
  }
}

static void particle_exclusion_sanity_checks(int part1, int part2) {
  if (part1 == part2) {
    throw std::runtime_error("Particles cannot exclude themselves (id " +
                             std::to_string(part1) + ")");
  }
  check_particle_exists(part1);
  check_particle_exists(part2);
}

void remove_particle_exclusion(int part1, int part2) {
  particle_exclusion_sanity_checks(part1, part2);
  mpi_call_all(mpi_remove_exclusion_local, part1, part2);
}

void add_particle_exclusion(int part1, int part2) {
  particle_exclusion_sanity_checks(part1, part2);
  mpi_call_all(mpi_add_exclusion_local, part1, part2);
}

void auto_exclusions(int distance) {
  /* partners is a list containing the currently found excluded particles for
     each particle, and their distance, as an interleaved list */
  std::unordered_map<int, std::vector<int>> partners;

  /* determine initial connectivity */
  for (auto const &part1 : partCfg()) {
    auto const p1 = part1.id();
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
      auto const p1 = p.id();
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
        add_particle_exclusion(id, j);
  }
}
#endif // EXCLUSIONS

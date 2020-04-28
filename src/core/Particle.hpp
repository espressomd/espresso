/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#ifndef ESPRESSO_CORE_PARTICLE_HPP
#define ESPRESSO_CORE_PARTICLE_HPP

#include "config.hpp"

#include "BondList.hpp"

#include <utils/Vector.hpp>
#include <utils/math/quaternion.hpp>

#include <boost/serialization/vector.hpp>

#include <cstdint>

enum : uint8_t {
  ROTATION_FIXED = 0u,
  ROTATION_X = 1u,
  ROTATION_Y = 2u,
  ROTATION_Z = 4u
};

struct ParticleParametersSwimming {
  bool swimming = false;
  double f_swim = 0.;
  double v_swim = 0.;
  int push_pull = 0;
  double dipole_length = 0.;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &swimming &f_swim &v_swim &push_pull &dipole_length;
  }
};

/** Properties of a particle which are not supposed to
 *  change during the integration, but have to be known
 *  for all ghosts. Ghosts are particles which are
 *  needed in the interaction calculation, but are just copies of
 *  particles stored on different nodes.
 */
struct ParticleProperties {
  /** unique identifier for the particle. */
  int identity = -1;
  /** Molecule identifier. */
  int mol_id = 0;
  /** particle type, used for non bonded interactions. */
  int type = 0;

  /** particle mass */
#ifdef MASS
  double mass = 1.0;
#else
  constexpr static double mass{1.0};
#endif /* MASS */

  /** rotational inertia */
#ifdef ROTATIONAL_INERTIA
  Utils::Vector3d rinertia = {1., 1., 1.};
#else
  static constexpr Utils::Vector3d rinertia = {1., 1., 1.};
#endif

  /** bitfield for the particle axes of rotation */
#ifdef ROTATION
  uint8_t rotation = ROTATION_FIXED;
#else
  static constexpr uint8_t rotation = ROTATION_FIXED;
#endif

  /** charge. */
#ifdef ELECTROSTATICS
  double q = 0.0;
#else
  constexpr static double q{0.0};
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
  /** electrophoretic mobility times E-field: mu_0 * E */
  Utils::Vector3d mu_E = {0., 0., 0.};
#endif

#ifdef DIPOLES
  /** dipole moment (absolute value) */
  double dipm = 0.;
#endif

#ifdef VIRTUAL_SITES
  /** is particle virtual */
  bool is_virtual = false;
#ifdef VIRTUAL_SITES_RELATIVE
  /** The following properties define, with respect to which real particle a
   *  virtual site is placed and at what distance. The relative orientation of
   *  the vector pointing from real particle to virtual site with respect to the
   *  orientation of the real particle is stored in the virtual site's
   *  quaternion attribute.
   */
  struct VirtualSitesRelativeParameters {
    int to_particle_id = 0;
    double distance = 0;
    /** Relative position of the virtual site. */
    Utils::Vector4d rel_orientation = {0., 0., 0., 0.};
    /** Orientation of the virtual particle in the body fixed frame. */
    Utils::Vector4d quat = {0., 0., 0., 0.};

    template <class Archive> void serialize(Archive &ar, long int) {
      ar &to_particle_id;
      ar &distance;
      ar &rel_orientation;
      ar &quat;
    }
  } vs_relative;
#endif
#else  /* VIRTUAL_SITES */
  static constexpr bool is_virtual = false;
#endif /* VIRTUAL_SITES */

#ifdef LANGEVIN_PER_PARTICLE
  double T = -1.;
#ifndef PARTICLE_ANISOTROPY
  double gamma = -1.;
#else
  Utils::Vector3d gamma = {-1., -1., -1.};
#endif // PARTICLE_ANISOTROPY
/** Friction coefficient gamma for rotation */
#ifdef ROTATION
#ifndef PARTICLE_ANISOTROPY
  double gamma_rot = -1.;
#else
  Utils::Vector3d gamma_rot = {-1., -1., -1.};
#endif // ROTATIONAL_INERTIA
#endif // ROTATION
#endif // LANGEVIN_PER_PARTICLE

#ifdef EXTERNAL_FORCES
  /** flag whether to fix a particle in space.
      Values:
      <ul> <li> 0 no external influence
           <li> 1 apply external force \ref ParticleProperties::ext_force
           <li> 2,3,4 fix particle coordinate 0,1,2
           <li> 5 apply external torque \ref ParticleProperties::ext_torque
      </ul>
  */
  uint8_t ext_flag = 0;
  /** External force, apply if \ref ParticleProperties::ext_flag == 1. */
  Utils::Vector3d ext_force = {0, 0, 0};

#ifdef ROTATION
  /** External torque, apply if \ref ParticleProperties::ext_flag == 16. */
  Utils::Vector3d ext_torque = {0, 0, 0};
#endif
#else
  static constexpr const uint8_t ext_flag =
      0; // no external forces and fixed coordinates
#endif

#ifdef ENGINE
  ParticleParametersSwimming swim;
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &identity;
    ar &mol_id;
    ar &type;
#ifdef MASS
    ar &mass;
#endif /* MASS */
#ifdef ROTATIONAL_INERTIA
    ar &rinertia;
#endif
#ifdef ROTATION
    ar &rotation;
#endif
#ifdef ELECTROSTATICS
    ar &q;
#endif

#ifdef LB_ELECTROHYDRODYNAMICS
    ar &mu_E;
#endif
#ifdef DIPOLES
    ar &dipm;
#endif

#ifdef VIRTUAL_SITES
    ar &is_virtual;
#ifdef VIRTUAL_SITES_RELATIVE
    ar &vs_relative;
#endif
#endif /* VIRTUAL_SITES */

#ifdef LANGEVIN_PER_PARTICLE
    ar &T;
    ar &gamma;
#ifdef ROTATION
    ar &gamma_rot;
#endif
#endif // LANGEVIN_PER_PARTICLE
#ifdef EXTERNAL_FORCES
    ar &ext_flag;
    ar &ext_force;
#ifdef ROTATION
    ar &ext_torque;
#endif
#endif

#ifdef ENGINE
    ar &swim;
#endif
  }
};

/** Positional information on a particle. Information that is
 *  communicated to calculate interactions with ghost particles.
 */
struct ParticlePosition {
  /** periodically folded position. */
  Utils::Vector3d p = {0, 0, 0};

#ifdef ROTATION
  /** quaternion to define particle orientation */
  Utils::Vector4d quat = {1., 0., 0., 0.};
  /** unit director calculated from the quaternion */
  Utils::Vector3d calc_director() const {
    return Utils::convert_quaternion_to_director(quat);
  };
#endif

#ifdef BOND_CONSTRAINT
  /** particle position at the previous time step */
  Utils::Vector3d p_old = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &p;
#ifdef ROTATION
    ar &quat;
#endif
#ifdef BOND_CONSTRAINT
    ar &p_old;
#endif
  }
};

/** Force information on a particle. Forces of ghost particles are
 *  collected and added up to the force of the original particle.
 */
struct ParticleForce {
  ParticleForce() = default;
  ParticleForce(ParticleForce const &) = default;
  ParticleForce(const Utils::Vector3d &f) : f(f) {}
#ifdef ROTATION
  ParticleForce(const Utils::Vector3d &f, const Utils::Vector3d &torque)
      : f(f), torque(torque) {}
#endif

  friend ParticleForce operator+(ParticleForce const &lhs,
                                 ParticleForce const &rhs) {
#ifdef ROTATION
    return {lhs.f + rhs.f, lhs.torque + rhs.torque};
#else
    return lhs.f + rhs.f;
#endif
  }

  ParticleForce &operator+=(ParticleForce const &rhs) {
    return *this = *this + rhs;
  }

  /** force. */
  Utils::Vector3d f = {0., 0., 0.};

#ifdef ROTATION
  /** torque */
  Utils::Vector3d torque = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &f;
#ifdef ROTATION
    ar &torque;
#endif
  }
};

/** Momentum information on a particle. Information not contained in
    communication of ghost particles so far, but a communication would
    be necessary for velocity dependent potentials. */
struct ParticleMomentum {
  /** velocity. */
  Utils::Vector3d v = {0., 0., 0.};

#ifdef ROTATION
  /** angular velocity
      ALWAYS IN PARTICLE FIXED, I.E., CO-ROTATING COORDINATE SYSTEM */
  Utils::Vector3d omega = {0., 0., 0.};
#endif

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &v;
#ifdef ROTATION
    ar &omega;
#endif
  }
};

/** Information on a particle that is needed only on the
 *  node the particle belongs to
 */
struct ParticleLocal {
  /** check whether a particle is a ghost or not */
  bool ghost = false;
  /** position in the last time step before last Verlet list update. */
  Utils::Vector3d p_old = {0, 0, 0};
  /** index of the simulation box image where the particle really sits. */
  Utils::Vector3i i = {0, 0, 0};

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &ghost;
    ar &p_old;
    ar &i;
  }
};

/** Struct holding all information for one particle. */
struct Particle { // NOLINT(bugprone-exception-escape)
  int &identity() { return p.identity; }
  int const &identity() const { return p.identity; }

  bool operator==(Particle const &rhs) const {
    return identity() == rhs.identity();
  }

  bool operator!=(Particle const &rhs) const {
    return identity() != rhs.identity();
  }

  ///
  ParticleProperties p;
  ///
  ParticlePosition r;
#ifdef DIPOLES
  Utils::Vector3d calc_dip() const { return r.calc_director() * p.dipm; }
#endif
  ///
  ParticleMomentum m;
  ///
  ParticleForce f;
  ///
  ParticleLocal l;

private:
  BondList bl;

public:
  auto &bonds() { return bl; }
  auto const &bonds() const { return bl; }

private:
#ifdef EXCLUSIONS
  /** list of particles, with which this particle has no nonbonded
   *  interactions
   */

  std::vector<int> el;
#endif

public:
  std::vector<int> &exclusions() {
#ifdef EXCLUSIONS
    return el;
#else
    throw std::runtime_error{"Exclusions not enabled."};
#endif
  }

  std::vector<int> const &exclusions() const {
#ifdef EXCLUSIONS
    return el;
#else
    throw std::runtime_error{"Exclusions not enabled."};
#endif
  }

private:
  friend boost::serialization::access;
  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar &p;
    ar &r;
    ar &m;
    ar &f;
    ar &l;
    ar &bl;
#ifdef EXCLUSIONS
    ar &el;
#endif
  }
};

#endif

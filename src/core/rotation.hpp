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
#ifndef ROTATION_H
#define ROTATION_H
/** \file
 *  This file contains all subroutines required to process rotational motion.
 */

#include "config.hpp"

#ifdef ROTATION

#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/math/rotation_matrix.hpp>

/** Propagate angular velocities and update quaternions on a particle */
void propagate_omega_quat_particle(Particle &p);

/** Convert torques to the body-fixed frame and propagate
    angular velocities */
void convert_torques_propagate_omega(const ParticleRange &particles);

/** Convert torques to the body-fixed frame to start
    the integration loop */
void convert_initial_torques(const ParticleRange &particles);

// Frame conversion routines
inline Utils::Vector3d
convert_vector_body_to_space(const Particle &p, const Utils::Vector3d &vec) {
  return rotation_matrix(p.r.quat) * vec;
}

inline Utils::Vector3d convert_vector_space_to_body(const Particle &p,
                                                    const Utils::Vector3d &v) {
  return transpose(rotation_matrix(p.r.quat)) * v;
}

/**
 * @brief Transform matrix from body- to space-fixed frame.
 *
 * Given a linear map represented by \f$ A \in \mathbb{R}^{3 \times 3}\f$
 * in the body-fixed frame, this returns the matrix \f$ A \in \mathbb{R}^{3
 * \times 3}\f$ representing the map in the space-fixed frame. They are related
 * by the map between the space-fixed and body-fixed frame \f$O\f$ like
 *
 * \f[
 *     A' = O^T A O.
 * \f]
 *
 * @tparam T Scalar type
 * @param p Particle transforming from.
 * @param A Matrix to transform
 * @return Matrix representation in space-fixed coordinates.
 */
template <class T>
auto convert_body_to_space(const Particle &p, const Utils::Matrix<T, 3, 3> &A) {
  auto const O = rotation_matrix(p.r.quat);
  return transpose(O) * A * O;
}

#ifdef DIPOLES

/** convert a dipole moment to quaternions and dipolar strength  */
inline std::pair<Utils::Vector4d, double>
convert_dip_to_quat(const Utils::Vector3d &dip) {
  auto quat = Utils::convert_director_to_quaternion(dip);
  return {quat, dip.norm()};
}

#endif

/** Rotate the particle p around the NORMALIZED axis a by amount phi */
void local_rotate_particle(Particle &p, const Utils::Vector3d &a, double phi);

#endif // ROTATION
#endif

/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

/** \file
 *  Ghost particles and particle exchange.
 *
 *  In this file you find the payload-selector flags used by the
 *  ghost-communication engine (@ref GhostComm::halo_exchange).
 */

#include <config/config.hpp>

/** Transfer data classes, for the ghost-communication engine */
enum : unsigned {
  GHOSTTRANS_NONE = 0u,
  /// transfer \ref ParticleProperties
  GHOSTTRANS_PROPRTS = 1u,
  /// transfer \ref ParticlePosition
  GHOSTTRANS_POSITION = 2u,
  /// transfer \ref ParticleMomentum
  GHOSTTRANS_MOMENTUM = 8u,
  /// transfer \ref ParticleForce
  GHOSTTRANS_FORCE = 16u,
#ifdef ESPRESSO_BOND_CONSTRAINT
  /// transfer \ref ParticleRattle
  GHOSTTRANS_RATTLE = 32u,
#endif
  /// resize the receiver particle arrays to the size of the senders
  GHOSTTRANS_PARTNUM = 64u,
  GHOSTTRANS_BONDS = 128u,
#ifdef ESPRESSO_ROTATION
  /// transfer orientation quaternion (pushed with position;
  /// runtime-conditional)
  GHOSTTRANS_QUAT = 256u,
  /// transfer torque (reduced with force; runtime-conditional)
  GHOSTTRANS_TORQUE = 512u,
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
  /// transfer dipole field tracking data
  GHOSTTRANS_DIPFLD = 1024u,
#endif
#ifdef ESPRESSO_BOND_CONSTRAINT
  /// transfer the position at the start of the time step (pushed with
  /// position; only needed by RATTLE)
  GHOSTTRANS_POS_LAST = 2048u,
#endif
};

/// Union of all transfer data classes of this build (add new flags here)
inline constexpr unsigned GHOSTTRANS_ALL =
    GHOSTTRANS_PROPRTS | GHOSTTRANS_POSITION | GHOSTTRANS_MOMENTUM |
    GHOSTTRANS_FORCE | GHOSTTRANS_PARTNUM | GHOSTTRANS_BONDS
#ifdef ESPRESSO_BOND_CONSTRAINT
    | GHOSTTRANS_RATTLE | GHOSTTRANS_POS_LAST
#endif
#ifdef ESPRESSO_ROTATION
    | GHOSTTRANS_QUAT | GHOSTTRANS_TORQUE
#endif
#ifdef ESPRESSO_DIPOLE_FIELD_TRACKING
    | GHOSTTRANS_DIPFLD
#endif
    ;

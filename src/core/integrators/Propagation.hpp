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

#include "PropagationMode.hpp"

class Propagation {
public:
  int integ_switch = INTEG_METHOD_NVT;
  int used_propagations = PropagationMode::NONE;
  int default_propagation = PropagationMode::NONE;
  int lb_skipped_md_steps = 0;
  int ek_skipped_md_steps = 0;
  /** If true, forces will be recalculated before the next integration. */
  bool recalc_forces = true;
  bool recalc_used_propagations = true;

  void update_default_propagation(int thermo_switch);

  /** @brief Decide whether a particle with propagation bits @p particle_prop
   *  participates in propagation mode @p mode.
   *
   *  Overload taking the propagation value directly: hot per-particle loops
   * read
   *  @c p.propagation() ONCE into a local and pass it here for all ~10 mode
   *  queries, instead of re-reading the ParticleStore propagation column twice
   *  per query.
   */
  bool should_propagate_with(int particle_prop, int mode) const {
    return (particle_prop & mode) or
           ((default_propagation & mode) and
            (particle_prop & PropagationMode::SYSTEM_DEFAULT));
  }

  template <typename Particle>
  bool should_propagate_with(Particle const &p, int mode) const {
    return should_propagate_with(p.propagation(), mode);
  }

  void set_integ_switch(int value) {
    integ_switch = value;
    recalc_forces = true;
    recalc_used_propagations = true;
  }

  /** True for all integrators that use inertial dynamics (VV, SE, NPT).
   *  False for non-inertial ones (BD, SD, steepest descent), which cannot
   *  use RATTLE rigid bonds. */
  bool is_inertial() const {
    return integ_switch != INTEG_METHOD_STEEPEST_DESCENT &&
           integ_switch != INTEG_METHOD_BD && integ_switch != INTEG_METHOD_SD;
  }
};

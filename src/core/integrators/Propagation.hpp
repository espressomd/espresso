/*
 * Copyright (C) 2010-2023 The ESPResSo project
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

  void update_default_propagation();

  template <typename Particle>
  bool should_propagate_with(Particle const &p, int mode) const {
    return (p.propagation() & mode) or
           ((default_propagation & mode) and
            (p.propagation() & PropagationMode::SYSTEM_DEFAULT));
  }

  void set_integ_switch(int value) {
    integ_switch = value;
    recalc_forces = true;
  }
};

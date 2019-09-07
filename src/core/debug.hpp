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
 *  This file controls debug facilities.
 *
 *  Implementation in debug.cpp.
 *
 *  For every define there exists a macro that can be used to encapsulate
 *  short lines (like printf("...",...);) of code that should be executed
 *  iff the respective *_DEBUG macro is defined.
 */

#include "config.hpp"

/** this performs a lot of tests which will very likely detect corruptions of
 *  \ref local_particles and the cell structure.
 */
void check_particle_consistency();

/** check the consistency of the cells and particle_node.
 */
void check_particles();

void check_particle_sorting();

/** Print all particle positions contained in \ref cells::cells array. */
void print_particle_positions();
/** Print all particle forces contained in \ref cells::cells array. */
void print_particle_forces();

extern int this_node;

/** by setting this variable to 1, a regular exit is
 *  indicated. In that case, no core dump is generated.
 */
extern int regular_exit;

#ifdef EVENT_DEBUG
#define EVENT_TRACE(cmd)                                                       \
  { cmd; }
#else
/** Equals { cmd } iff EVENT_DEBUG is set. */
#define EVENT_TRACE(cmd)
#endif

#ifdef HALO_DEBUG
#define HALO_TRACE(cmd)                                                        \
  { cmd; }
#else
/** Equals { cmd  } iff HALO_DEBUG is set. */
#define HALO_TRACE(cmd)
#endif

#ifdef P3M_DEBUG
#define P3M_TRACE(cmd)                                                         \
  { cmd; }
#else
/** Equals { cmd } iff P3M_DEBUG is set. */
#define P3M_TRACE(cmd)
#endif

#ifdef THERMO_DEBUG
#define THERMO_TRACE(cmd)                                                      \
  { cmd; }
#else
#define THERMO_TRACE(cmd)
#endif

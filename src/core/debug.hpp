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

/** check the consistency of the cells and particle_node. Called from
 *  mpi_bcast_event(CHECK_PARTICLES)
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

/** Identity of the particle to check extensively if ONEPART_DEBUG is defined.
 */
extern int check_id;

#ifdef ESIF_DEBUG
#define ESIF_TRACE(cmd)                                                        \
  { cmd; }
#else
#define ESIF_TRACE(cmd)
#endif

#ifdef COMM_DEBUG
#define COMM_TRACE(cmd)                                                        \
  { cmd; }
#else
/** Equals { cmd } iff COMM_DEBUG is set. */
#define COMM_TRACE(cmd)
#define COMM_DBG_MSG(msg)
#endif

#ifdef EVENT_DEBUG
#define EVENT_TRACE(cmd)                                                       \
  { cmd; }
#else
/** Equals { cmd } iff EVENT_DEBUG is set. */
#define EVENT_TRACE(cmd)
#endif

#ifdef PARTICLE_DEBUG
#define PART_TRACE(cmd)                                                        \
  { cmd; }
#else
/** Equals { cmd } iff PARTICLE_DEBUG is set. */
#define PART_TRACE(cmd)
#endif

#ifdef INTEG_DEBUG
#define INTEG_TRACE(cmd)                                                       \
  { cmd; }
#else
/** Equals { cmd } iff INTEG_DEBUG is set. */
#define INTEG_TRACE(cmd)
#endif

#ifdef CELL_DEBUG
#define CELL_TRACE(cmd)                                                        \
  { cmd; }
#else
/** Equals { cmd } iff CELL_DEBUG is set. */
#define CELL_TRACE(cmd)
#endif

#ifdef GHOST_DEBUG
#define GHOST_TRACE(cmd)                                                       \
  { cmd; }
#else
/** Equals { cmd } iff GHOST_DEBUG is set. */
#define GHOST_TRACE(cmd)
#endif

#ifdef HALO_DEBUG
#define HALO_TRACE(cmd)                                                        \
  { cmd; }
#else
/** Equals { cmd  } iff HALO_DEBUG is set. */
#define HALO_TRACE(cmd)
#endif

#ifdef GRID_DEBUG
#define GRID_TRACE(cmd)                                                        \
  { cmd; }
#else
/** Equals { cmd } iff GRID_DEBUG is set. */
#define GRID_TRACE(cmd)
#endif

#ifdef LATTICE_DEBUG
#define LATTICE_TRACE(cmd)                                                     \
  { cmd; }
#else
/** Equals { cmd } iff LATTICE_DEBUG is set. */
#define LATTICE_TRACE(cmd)
#endif

#ifdef FORCE_DEBUG
#define FORCE_TRACE(cmd)                                                       \
  { cmd; }
#else
/** Equals { cmd } iff FORCE_DEBUG is set. */
#define FORCE_TRACE(cmd)
#endif

#ifdef VERLET_DEBUG
#define VERLET_TRACE(cmd)                                                      \
  { cmd; }
#else
/** Equals { cmd } iff VERLET_DEBUG is set. */
#define VERLET_TRACE(cmd)
#endif

#ifdef P3M_DEBUG
#define P3M_TRACE(cmd)                                                         \
  { cmd; }
#else
/** Equals { cmd } iff P3M_DEBUG is set. */
#define P3M_TRACE(cmd)
#endif

#ifdef MDLC_DEBUG
#define MDLC_TRACE(cmd)                                                        \
  { cmd; }
#else
#define MDLC_TRACE(cmd)
#endif

#ifdef FFT_DEBUG
#define FFT_TRACE(cmd)                                                         \
  { cmd; }
#else
/** Equals { cmd } iff FFT_DEBUG is set. */
#define FFT_TRACE(cmd)
#endif

#ifdef RANDOM_DEBUG
#define RANDOM_TRACE(cmd)                                                      \
  { cmd; }
#else
#define RANDOM_TRACE(cmd)
#endif

#ifdef THERMO_DEBUG
#define THERMO_TRACE(cmd)                                                      \
  { cmd; }
#else
#define THERMO_TRACE(cmd)
#endif

#ifdef LJ_DEBUG
#define LJ_TRACE(cmd)                                                          \
  { cmd; }
#else
#define LJ_TRACE(cmd)
#endif

#ifdef MORSE_DEBUG
#define MORSE_TRACE(cmd)                                                       \
  { cmd; }
#else
#define MORSE_TRACE(cmd)
#endif

#ifdef BUCK_DEBUG
#define BUCK_TRACE(cmd)                                                        \
  { cmd; }
#else
#define BUCK_TRACE(cmd)
#endif

#ifdef ESR_DEBUG
#define ESR_TRACE(cmd)                                                         \
  { cmd; }
#else
#define ESR_TRACE(cmd)
#endif

#ifdef ESK_DEBUG
#define ESK_TRACE(cmd)                                                         \
  { cmd; }
#else
#define ESK_TRACE(cmd)
#endif

#ifdef FENE_DEBUG
#define FENE_TRACE(cmd)                                                        \
  { cmd; }
#else
#define FENE_TRACE(cmd)
#endif

#ifdef GHOST_FORCE_DEBUG
#define GHOST_FORCE_TRACE(cmd)                                                 \
  { cmd; }
#else
#define GHOST_FORCE_TRACE(cmd)
#endif

#ifdef ONEPART_DEBUG
#define ONEPART_TRACE(cmd)                                                     \
  { cmd; }
#else
#define ONEPART_TRACE(cmd)
#endif

#ifdef STAT_DEBUG
#define STAT_TRACE(cmd)                                                        \
  { cmd; }
#else
/** Equals { cmd } iff STAT_DEBUG is set. */
#define STAT_TRACE(cmd)
#endif

#ifdef POLY_DEBUG
#define POLY_TRACE(cmd)                                                        \
  { cmd; }
#else
/** Equals { cmd } iff POLY_DEBUG is set. */
#define POLY_TRACE(cmd)
#endif

#ifdef PTENSOR_DEBUG
#define PTENSOR_TRACE(cmd)                                                     \
  { cmd; }
#else
#define PTENSOR_TRACE(cmd)
#endif

#ifdef LB_DEBUG
#define LB_TRACE(cmd)                                                          \
  { cmd; }
#else
/** Equals { cmd } iff LB_DEBUG is set. */
#define LB_TRACE(cmd)
#endif

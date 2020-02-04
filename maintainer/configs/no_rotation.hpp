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
/* This is the default myconfig.hpp-file. If no other myconfig-file is
   found, this file is used.

   DO NOT MODIFY THIS FILE! It should be modified *only* by the
   maintainers of ESPResSo, as it has a profound impact on many users,
   in particular newbies.
*/

// Geometry, equation of motion, thermostat/barostat
#define MASS
#define EXTERNAL_FORCES
#define LANGEVIN_PER_PARTICLE
#define BROWNIAN_PER_PARTICLE
#define BOND_CONSTRAINT
#define NPT
#define DPD

// Charges and dipoles
#define ELECTROSTATICS
#ifdef CUDA
#define MMM1D_GPU
#endif

// Hydrodynamics
#define LB_BOUNDARIES
#ifdef CUDA
#define LB_BOUNDARIES_GPU
#endif

// Electrokinetics
#ifdef CUDA
#define ELECTROKINETICS
#define EK_BOUNDARIES
#endif

// Force/energy calculation
#define EXCLUSIONS

#define TABULATED
#define LENNARD_JONES
#define LENNARD_JONES_GENERIC
#define LJGEN_SOFTCORE
#define LJCOS
#define LJCOS2
#define GAUSSIAN
#define HAT
#define SMOOTH_STEP
#define HERTZIAN
#define SOFT_SPHERE
#define WCA
#define THOLE

// Further features
#define VIRTUAL_SITES_INERTIALESS_TRACERS
#define COLLISION_DETECTION

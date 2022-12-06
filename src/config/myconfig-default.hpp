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

/*
 * This is the default myconfig.hpp file.
 * When users don't supply a myconfig.hpp file, this file is used.
 *
 * DO NOT MODIFY THIS FILE! It should be modified *only* by the
 * maintainers of ESPResSo, as it has a profound impact on many users.
 */

// Geometry, equation of motion, thermostat/barostat
#define ROTATION
#define ROTATIONAL_INERTIA
#define MASS
#define PARTICLE_ANISOTROPY
#define EXTERNAL_FORCES
#define THERMOSTAT_PER_PARTICLE
#define BOND_CONSTRAINT
#define NPT
#define DPD

// Charges and dipoles
#define ELECTROSTATICS
#ifdef CUDA
#define MMM1D_GPU
#endif
#define DIPOLES

// Active matter
#define ENGINE

// Force/energy calculation
#define EXCLUSIONS

// Long-range interactions
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

#ifdef FFTW
#define THOLE
#endif

// Further features
#define VIRTUAL_SITES_RELATIVE
#define VIRTUAL_SITES_INERTIALESS_TRACERS
#define COLLISION_DETECTION

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
/* maximal set of features usable at the same time */
#define ELECTROSTATICS
#define DIPOLES
#ifdef SCAFACOS
#define SCAFACOS_DIPOLES
#endif
#define ROTATION
#define ROTATIONAL_INERTIA
#define PARTICLE_ANISOTROPY
#define EXTERNAL_FORCES

#define MASS
#define EXCLUSIONS

#define BOND_CONSTRAINT
#define COLLISION_DETECTION
#define LANGEVIN_PER_PARTICLE
#define BROWNIAN_PER_PARTICLE
#define SWIMMER_REACTIONS

#define NPT

#define DPD

#define LB_BOUNDARIES
#define LB_ELECTROHYDRODYNAMICS

#define ENGINE

#ifdef CUDA
#define LB_BOUNDARIES_GPU
#define ELECTROKINETICS
#define EK_BOUNDARIES
#define MMM1D_GPU
#endif

#define TABULATED
#define LENNARD_JONES
#define LENNARD_JONES_GENERIC
#define LJGEN_SOFTCORE
#define LJCOS
#define LJCOS2
#define GAUSSIAN
#define HAT
#define GAY_BERNE
#define SMOOTH_STEP
#define HERTZIAN
#define BMHTF_NACL
#define MORSE
#define BUCKINGHAM
#define SOFT_SPHERE
#define WCA
#define THOLE

#define EXPERIMENTAL_FEATURES

#define VIRTUAL_SITES_RELATIVE
#define VIRTUAL_SITES_INERTIALESS_TRACERS

#define ADDITIONAL_CHECKS

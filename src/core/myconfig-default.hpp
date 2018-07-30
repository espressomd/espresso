/*
  Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project
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
/* This is the default myconfig.hpp-file. If no other myconfig-file is
   found, this file is used.

   DO NOT MODIFY THIS FILE! It should be modified *only* by the
   maintainers of ESPResSo, as it has a profound impact on many users,
   in particular newbies.
*/
/* global features */
#define PARTIAL_PERIODIC
#define ELECTROSTATICS
#define EXTERNAL_FORCES
#define MASS
#define EXCLUSIONS
#define NPT
#define COLLISION_DETECTION
#define LANGEVIN_PER_PARTICLE

/* potentials */
#define TABULATED
#define LENNARD_JONES
#define LENNARD_JONES_GENERIC
#define MORSE
#define LJCOS
#define LJCOS2
#define BUCKINGHAM
#define SOFT_SPHERE
#define BOND_ANGLE
#define GAUSSIAN
#define HERTZIAN

// Lattice Boltzmann
#define LB
#define LB_BOUNDARIES
#ifdef CUDA
  #define LB_GPU
  #define LB_BOUNDARIES_GPU
#endif  

// Electrokinetics
#ifdef CUDA
  #define ELECTROKINETICS
  #define EK_BOUNDARIES
  #define EK_ELECTROSTATIC_COUPLING
#endif


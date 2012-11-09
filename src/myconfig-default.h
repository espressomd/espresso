/*
  Copyright (C) 2010,2011,2012 The ESPResSo project
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
/* This is the default myconfig.h-file. If no other myconfig-file is
   found, this file is used.
*/
/* global features */
#define PARTIAL_PERIODIC
#define ELECTROSTATICS
#define EXTERNAL_FORCES
#define CONSTRAINTS
#define MASS
#define EXCLUSIONS
#define COMFORCE
#define COMFIXED
#define NPT

/* potentials */
#define TABULATED
#define LENNARD_JONES
#define LENNARD_JONES_GENERIC
#define MORSE
#define LJCOS
#define LJCOS2
#define BUCKINGHAM
#define SOFT_SPHERE
#define BOND_ANGLE_COSINE

/* NEW Interactions */
#define AREA_FORCE_LOCAL
#define AREA_FORCE_GLOBAL   
#define STRETCHING_FORCE
#define BENDING_FORCE
#define VOLUME_FORCE	

/* Bond angle */
/* Note: Activate ONLY ONE bonded angle potential out of the following! */
//#define BOND_ANGLE_HARMONIC
//#define BOND_ANGLE_COSINE
//#define BOND_ANGLE_COSSQUARE

//#define BOND_ANGLEDIST
//#define BOND_ANGLEDIST_HARMONIC

//#define BOND_ENDANGLEDIST
//#define BOND_ENDANGLEDIST_HARMONIC

/* Strange features. Use only if you know what you are doing! */
/* activate the old dihedral form */
//#define OLD_DIHEDRAL
/* turn off nonbonded interactions within molecules */
//#define NO_INTRA_NB

/* Debugging */
//#define ADDITIONAL_CHECKS
//#define ASYNC_BARRIER

//#define COMM_DEBUG
//#define EVENT_DEBUG
//#define INTEG_DEBUG
//#define CELL_DEBUG
//#define GHOST_DEBUG
//#define LATTICE_DEBUG
//#define HALO_DEBUG
//#define GRID_DEBUG
//#define VERLET_DEBUG
//#define PARTICLE_DEBUG
//#define P3M_DEBUG
//#define EWALD_DEBUG
//#define FFT_DEBUG
//#define RANDOM_DEBUG
//#define FORCE_DEBUG
//#define THERMO_DEBUG
//#define LJ_DEBUG
//#define MORSE_DEBUG
//#define ESR_DEBUG
//#define ESK_DEBUG
//#define FENE_DEBUG
//#define GHOST_FORCE_DEBUG
//#define STAT_DEBUG
//#define POLY_DEBUG
//#define MOLFORCES_DEBUG
//#define PTENSOR_DEBUG
//#define MEM_DEBUG
//#define MAGGS_DEBUG
//#define LB_DEBUG
//#define VIRTUAL_SITES_DEBUG

//#define MPI_CORE
//#define FORCE_CORE

/* Single particle debugging */
//#define ONEPART_DEBUG
// which particle id to debug
//#define ONEPART_DEBUG_ID 3149


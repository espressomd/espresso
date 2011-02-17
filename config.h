/*
  Copyright (C) 2010,2011 The ESPResSo project
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
#ifndef CONFIG_H
#define CONFIG_H

/** \file config.h

    This file contains the defaults for Espresso. To modify them, add
    an appropriate line in myconfig.h. To find a list of features that
    can be compiled into Espresso, refer to myconfig-sample.h or to
    \ref tcl_features "the documentation of the features".
 */

/* Include the defines created by configure. */
#include <acconfig.h>

/************************************************/
/** \name Default Parameter Settings            */
/** These values can be changed from the Tcl    */
/** */
/************************************************/
/*@{*/

/** CELLS: Default value for the maximal number of cells per node. */
#define CELLS_MAX_NUM_CELLS 32768

/** P3M: Default for number of interpolation points of the charge
    assignment function. */
#define P3M_N_INTERPOL 32768

/** P3M: Default for boundary condition: Epsilon of the surrounding
    medium. */
#define P3M_EPSILON 0.0

/** P3M: Default for boundary condition in magnetic calculations */
#define P3M_EPSILON_MAGNETIC 1.0

/** P3M: Default for offset of first mesh point from the origin (left
    down corner of the simulation box. */
#define P3M_MESHOFF 0.5

/*@}*/

#ifdef MYCONFIG_H
#include MYCONFIG_H
#endif

/*********************************************************/
/** \name Parameters from myconfig.h that need to be set */
/*********************************************************/
/*@{*/

/** P3M: Number of Brillouin zones taken into account
    in the calculation of the optimal influence function (aliasing
    sums). */
#ifndef P3M_BRILLOUIN
#define P3M_BRILLOUIN 0
#endif
/** P3M: Maximal mesh size that will be checked. The current setting
         limits the memory consumption to below 1GB, which is probably
	 reasonable for a while. */
#ifndef P3M_MAX_MESH
#define P3M_MAX_MESH 128
#endif

/** Whether to use the approximation of Abramowitz/Stegun
    AS_erfc_part() for \f$\exp(d^2) erfc(d)\f$, or the C function erfc
    in P3M and Ewald summation. */
#ifndef USE_ERFC_APPROXIMATION
#define USE_ERFC_APPROXIMATION 1
#endif

/** Precision for capture of round off errors. */
#ifndef ROUND_ERROR_PREC
#define ROUND_ERROR_PREC 1.0e-14
#endif

/** Tiny angle cutoff for sinus calculations */
#ifndef TINY_SIN_VALUE
#define TINY_SIN_VALUE 1e-10
#endif
/** Tiny angle cutoff for cosine calculations */
#ifndef TINY_COS_VALUE
#define TINY_COS_VALUE 0.9999999999
#endif
/** Tiny length cutoff */
#ifndef TINY_LENGTH_VALUE
#define TINY_LENGTH_VALUE 0.0001
#endif

/** maximal number of iterations in the RATTLE algorithm before it bails out. */
#ifndef SHAKE_MAX_ITERATIONS
#define SHAKE_MAX_ITERATIONS 1000
#endif

/*@}*/

/*********************************************************/
/** \name  */
/*********************************************************/
/*@{*/

//inter_rf needs ELECTROSTATICS
#ifdef INTER_RF
#ifndef ELECTROSTATICS
#define ELECTROSTATICS
#endif
#endif

#ifdef GAY_BERNE
#ifndef ROTATION
#define ROTATION
#endif
#endif

/* activate P3M only with FFTW */
#if defined(ELECTROSTATICS) && defined(FFTW)
#define ELP3M
#endif


/* activate dipolar P3M only with FFTW */
#if defined(MAGNETOSTATICS) && defined(FFTW)
#define ELP3M
#ifndef DIPOLES
#define DIPOLES
#endif
#endif


/* MAGNETOSTATICS implies the use of DIPOLES */
#if defined(MAGNETOSTATICS)
#ifndef DIPOLES
#define DIPOLES
#endif
#endif

/* LB_ELECTROHYDRODYNAMICS needs LB, obviously... */
#ifdef LB_ELECTROHYDRODYNAMICS
#ifndef LB
#define LB
#endif
#endif

/* Lattice Boltzmann needs lattice structures and temporary particle data */
#ifdef LB
#define USE_TEMPORARY
#define LATTICE
//#define ALTERNATIVE_INTEGRATOR
#endif

//adress needs mol_cut
#ifdef ADRESS
#ifndef MOL_CUT
#define MOL_CUT
#endif
#endif

//mol_cut needs virtual sites
#ifdef MOL_CUT
#ifndef VIRTUAL_SITES_COM
#define VIRTUAL_SITES_COM
#endif
#endif

#if defined(DPD_MASS_RED) || defined(DPD_MASS_LIN)
#ifndef DPD_MASS
#define DPD_MASS
#endif
#endif

/*DPD with mass needs MASS and DPD */
#ifdef DPD_MASS
#ifndef MASS
#define MASS
#endif
#endif

/*Transversal DPD -> needs normal DPD*/
#ifdef TRANS_DPD
#ifndef DPD
#define DPD
#endif
#endif

/* If any bond angle potential is activated, actiavte the whole bond angle code */
#if defined(BOND_ANGLE_HARMONIC) || defined(BOND_ANGLE_COSINE) || defined(BOND_ANGLE_COSSQUARE)
#define BOND_ANGLE
#endif

/* If any bond angledist potential is activated, activate the whole bond angle code and constraints */
#if defined(BOND_ANGLEDIST_HARMONIC)
#define BOND_ANGLEDIST
#define CONSTRAINTS
#endif

#if defined(BOND_ENDANGLEDIST_HARMONIC)
#define BOND_ENDANGLEDIST
#define CONSTRAINTS
#endif

#if defined(VIRTUAL_SITES_COM) || defined(VIRTUAL_SITES_RELATIVE)
#define VIRTUAL_SITES
#endif

#ifdef VIRTUAL_SITES_RELATIVE
#ifndef ROTATION
#define ROTATION
#endif
#endif

/*@}*/

/********************************************/
/* \name exported functions of config.c     */
/********************************************/
/*@{*/
#include <tcl.h>

/** callback for version status. */
int tclcallback_version(Tcl_Interp *interp);

/** callback for compilation status. */
int tclcallback_compilation(Tcl_Interp *interp);
/*@}*/

#endif

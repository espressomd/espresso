/*
  Copyright (C) 2010,2011,2012,2013 The ESPResSo project
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
#ifndef _CONFIG_HPP
#define _CONFIG_HPP

/** \file config.hpp

    This file contains the defaults for Espresso. To modify them, add
    an appropriate line in myconfig.hpp. To find a list of features that
    can be compiled into Espresso, refer to myconfig-sample.hpp or to
    the documentation of the features.
 */

/* Include the defines created by configure. */
#include <acconfig.hpp>

/* Prevent C++ bindings in MPI (there is a DataType called LB in there) */
#define OMPI_SKIP_MPICXX
#define MPICH_SKIP_MPICXX

#include "myconfig-final.hpp"
#include "config-features.hpp"

extern const char* ESPRESSO_VERSION;

/*********************************************************/
/** \name Parameters from myconfig.hpp that need to be set */
/*********************************************************/
/*@{*/

#ifndef ONEPART_DEBUG_ID
#define ONEPART_DEBUG_ID 13
#endif

/** CELLS: Default value for the maximal number of cells per node. */
#ifndef CELLS_MAX_NUM_CELLS
#define CELLS_MAX_NUM_CELLS 32768
#endif

/** P3M: Default for number of interpolation points of the charge
    assignment function. */
#ifndef P3M_N_INTERPOL
#define P3M_N_INTERPOL 32768
#endif

/** P3M: Default for boundary condition: Epsilon of the surrounding
    medium. */
#ifndef P3M_EPSILON
#define P3M_EPSILON 0.0
#endif

/** P3M: Default for boundary condition in magnetic calculations */
#ifndef P3M_EPSILON_MAGNETIC
#define P3M_EPSILON_MAGNETIC 0.0
#endif

/** P3M: Default for offset of first mesh point from the origin (left
    down corner of the simulation box. */
#ifndef P3M_MESHOFF
#define P3M_MESHOFF 0.5
#endif

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

/** maximal number of objects in the object-in-fluid framework. */
#ifndef MAX_OBJECTS_IN_FLUID
#define MAX_OBJECTS_IN_FLUID 10000
#endif

/** number of fluid components for lattice boltzmann  */
#ifndef LB_COMPONENTS
#ifdef SHANCHEN
#define LB_COMPONENTS 2
#else 
#define LB_COMPONENTS 1
#endif
#endif


/* Mathematical constants, from gcc's math.hpp */
#ifndef M_PI
#define M_E		2.7182818284590452353602874713526625L  /* e */
#define M_LOG2E		1.4426950408889634073599246810018921L  /* log_2 e */
#define M_LOG10E	0.4342944819032518276511289189166051L  /* log_10 e */
#define M_LN2		0.6931471805599453094172321214581766L  /* log_e 2 */
#define M_LN10	2.3025850929940456840179914546843642L  /* log_e 10 */
#define M_PI		3.1415926535897932384626433832795029L  /* pi */
#define M_PI_2	1.5707963267948966192313216916397514L  /* pi/2 */
#define M_PI_4	0.7853981633974483096156608458198757L  /* pi/4 */
#define M_1_PI	0.3183098861837906715377675267450287L  /* 1/pi */
#define M_2_PI	0.6366197723675813430755350534900574L  /* 2/pi */
#define M_2_SQRTPI	1.1283791670955125738961589031215452L  /* 2/sqrt(pi) */
#define M_SQRT2	       	1.4142135623730950488016887242096981L  /* sqrt(2) */
#define M_SQRT1_2	0.7071067811865475244008443621048490L  /* 1/sqrt(2) */
#endif

#endif

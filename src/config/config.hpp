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
#ifndef ESPRESSO_CONFIG_HPP
#define ESPRESSO_CONFIG_HPP

/** \file
 *
 *  This file contains the defaults for ESPResSo. To modify them, add
 *  an appropriate line in myconfig.hpp. To find a list of features that
 *  can be compiled into ESPResSo, refer to myconfig-sample.hpp or to
 *  the documentation of the features.
 */

/* Prevent C++ bindings in MPI (there is a DataType called LB in there) */
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
#ifndef MPICH_SKIP_MPICXX
#define MPICH_SKIP_MPICXX
#endif

#include "config-features.hpp"

/** P3M: Default for offset of first mesh point from the origin (left
 *  down corner of the simulation box).
 */
#ifndef P3M_MESHOFF
#define P3M_MESHOFF 0.5
#endif

/** P3M: Number of Brillouin zones taken into account
 *  in the calculation of the optimal influence function (aliasing sums).
 */
#ifndef P3M_BRILLOUIN
#define P3M_BRILLOUIN 0
#endif
/** P3M: Maximal mesh size that will be checked. The current setting
 *  limits the memory consumption to below 1GB, which is probably
 *  reasonable for a while.
 */
#ifndef P3M_MAX_MESH
#define P3M_MAX_MESH 128
#endif

/** Whether to use the approximation of Abramowitz/Stegun @cite abramowitz65a
 *  @ref AS_erfc_part() for \f$\exp(d^2) \mathrm{erfc}(d)\f$,
 *  or the C function <tt>std::erfc()</tt> in P3M and Ewald summation.
 */
#ifndef USE_ERFC_APPROXIMATION
#define USE_ERFC_APPROXIMATION 1
#endif

/** Precision for capture of round off errors. */
#ifndef ROUND_ERROR_PREC
#define ROUND_ERROR_PREC 1.0e-14
#endif

/** Tiny angle cutoff for sinus calculations. */
#ifndef TINY_SIN_VALUE
#define TINY_SIN_VALUE 1e-10
#endif
/** Tiny angle cutoff for cosine calculations. */
#ifndef TINY_COS_VALUE
#define TINY_COS_VALUE 0.9999999999
#endif
/** Tiny length cutoff. */
#ifndef TINY_LENGTH_VALUE
#define TINY_LENGTH_VALUE 0.0001
#endif
/** Tiny oif elasticity cutoff. */
#ifndef TINY_OIF_ELASTICITY_COEFFICIENT
#define TINY_OIF_ELASTICITY_COEFFICIENT 1e-10
#endif

/** Maximal number of iterations in the RATTLE algorithm before it bails out. */
#ifndef SHAKE_MAX_ITERATIONS
#define SHAKE_MAX_ITERATIONS 1000
#endif

/** Maximal number of objects in the object-in-fluid framework. */
#ifndef MAX_OBJECTS_IN_FLUID
#define MAX_OBJECTS_IN_FLUID 10000
#endif

/** Maximal number of soft particles per immersed boundary */
#ifndef IBM_MAX_NUM
#define IBM_MAX_NUM 1000
#endif

#endif

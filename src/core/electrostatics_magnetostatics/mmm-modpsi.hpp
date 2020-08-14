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
/** \file
 *  Common parts of the MMM family of methods for the electrostatic
 *  interaction, MMM1D, MMM2D and ELC. This file contains the code for the
 *  CPU and GPU implementations.
 */

#ifndef MMM_BASE_H
#define MMM_BASE_H

#include <vector>

/** \name Math constants */
/*@{*/
#define C_2PI (2 * M_PI)
#define C_GAMMA (0.57721566490153286060651209008)
/*@}*/

/** Table of the Taylor expansions of the modified polygamma functions */
extern std::vector<std::vector<double>> modPsi;

/** Create both the even and odd polygamma functions up to order 2*n */
void create_mod_psi_up_to(int n);

#endif

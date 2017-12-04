/*
  Copyright (C) 2010,2012,2013,2014,2015,2016 The ESPResSo project
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

/** \file ljangle.hpp
 *  Routines to calculate the lennard-jones 12-10 with angular dependance.
 *  The potential is a product of a 12-10 LJ potential with two cos^2.
 *  The potential actually relies on 6 particles: the 2 primary beads
 *  and each bead needs two other particles to define an orientation.
 *  We calculate the interaction explicitly *without* the use of the ROTATION
 * feature.
 *
 *  Optional: simulate two different environments in which the interaction
 *  strengths differ. For example: simulate hydrogen-bonds both in water and
 *  inside a bilayer. The two environments are distinguished by their
 *  z-position. Input: the midplane of the 2nd environment, its total
 *  thickness, the thickness of the interface, and the value of the
 *  interaction strength in this medium. The interaction strengh of the second
 *  environment must be *stronger* than of the first one.
 *
 *  \ref forces.cpp
 */

#include "utils.hpp"

#ifdef LJ_ANGLE
#include "communication.hpp"
#include "interaction_data.hpp"
#include "ljangle.hpp"
#include <cmath>

int ljangle_set_params(int part_type_a, int part_type_b, double eps, double sig,
                       double cut, int b1p, int b1n, int b2p, int b2n,
                       double z0, double dz, double kappa,
                       double epsprime) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data)
    return ES_ERROR;

  data->LJANGLE_eps = eps;
  data->LJANGLE_sig = sig;
  data->LJANGLE_cut = cut;
  data->LJANGLE_bonded1pos = b1p;
  data->LJANGLE_bonded1neg = b1n;
  data->LJANGLE_bonded2pos = b2p;
  data->LJANGLE_bonded2neg = b2n;

  data->LJANGLE_bonded1type = part_type_a;

  if (dz > 0.) {
    data->LJANGLE_z0 = z0;
    data->LJANGLE_dz = dz;
    data->LJANGLE_kappa = kappa;
    data->LJANGLE_epsprime = epsprime;
  }

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif

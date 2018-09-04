/*
  Copyright (C) 2010-2018 The ESPResSo project
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

/** \file cos2.hpp
 *  Routines to calculate a flat potential with cosine tail energy and/or  force
 *  for a particle pair.  Cosine tail is different from that in ljcos.hpp
 *  Used for attractive tail/tail interactions in lipid bilayer calculations.
 *  Same potential as ljcos2 without Lennard-Jones part.
 *  \ref forces.cpp
 */
#include "cos2.hpp"

#ifdef COS2
#include <cmath>

#include "communication.hpp"
#include "lj.hpp"

int cos2_set_params(int part_type_a, int part_type_b, double eps, double offset,
                    double w) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data)
    return ES_ERROR;

  data->COS2_eps = eps;
  data->COS2_offset = offset;
  data->COS2_w = w;

  /* calculate dependent parameters */
  data->COS2_cut = w + offset;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif /* ifdef COS2 */

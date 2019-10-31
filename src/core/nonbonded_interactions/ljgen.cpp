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
 *
 *  Implementation of \ref ljgen.hpp
 */
#include "ljgen.hpp"

#ifdef LENNARD_JONES_GENERIC
#include "communication.hpp"
#include "lj.hpp"

#include <utils/constants.hpp>

int ljgen_set_params(int part_type_a, int part_type_b, double eps, double sig,
                     double cut, double shift, double offset, double a1,
                     double a2, double b1, double b2

#ifdef LJGEN_SOFTCORE
                     ,
                     double lambda, double softrad
#endif
) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data)
    return ES_ERROR;

  data->ljgen.eps = eps;
  data->ljgen.sig = sig;
  data->ljgen.cut = cut;
  data->ljgen.shift = shift;
  data->ljgen.offset = offset;
  data->ljgen.a1 = a1;
  data->ljgen.a2 = a2;
  data->ljgen.b1 = b1;
  data->ljgen.b2 = b2;
#ifdef LJGEN_SOFTCORE
  data->ljgen.lambda1 = lambda;
  data->ljgen.softrad = softrad;
#endif

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif

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
 *  Implementation of \ref ljcos.hpp
 */
#include "ljcos.hpp"

#ifdef LJCOS
#include "communication.hpp"

#include <utils/constants.hpp>

int ljcos_set_params(int part_type_a, int part_type_b, double eps, double sig,
                     double cut, double offset) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  if (!data)
    return ES_ERROR;

  data->ljcos.eps = eps;
  data->ljcos.sig = sig;
  data->ljcos.cut = cut;
  data->ljcos.offset = offset;

  /* Calculate dependent parameters */
  auto const driwu2 = 1.25992104989487316476721060728;
  auto const facsq = driwu2 * Utils::sqr(sig);
  data->ljcos.rmin = sqrt(driwu2) * sig;
  data->ljcos.alfa = Utils::pi() / (Utils::sqr(data->ljcos.cut) - facsq);
  data->ljcos.beta =
      Utils::pi() * (1. - (1. / (Utils::sqr(data->ljcos.cut) / facsq - 1.)));

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}
#endif

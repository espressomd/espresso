/*
  Copyright (C) 2018 The ESPResSo project

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

#include "wca.hpp"

#ifdef WCA
#include "communication.hpp"

#include <utils/constants.hpp>

int wca_set_params(int part_type_a, int part_type_b, double eps, double sig) {
  IA_parameters *data = get_ia_param_safe(part_type_a, part_type_b);

  data->WCA_eps = eps;
  data->WCA_sig = sig;
  data->WCA_cut = sig * std::pow(2., 1. / 6.);

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(part_type_a, part_type_b);

  return ES_OK;
}

#endif /* ifdef WCA */

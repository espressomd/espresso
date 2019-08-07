/*
Copyright (C) 2010-2018 The ESPResSo project

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
#include "npt.hpp"
#include "communication.hpp"
#include "config.hpp"

void synchronize_npt_state(int n_steps) {
  nptiso.invalidate_p_vel = 0;
  MPI_Bcast(&nptiso.p_inst, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&nptiso.p_diff, 1, MPI_DOUBLE, 0, comm_cart);
  MPI_Bcast(&nptiso.volume, 1, MPI_DOUBLE, 0, comm_cart);
  if (this_node == 0)
    nptiso.p_inst_av /= 1.0 * n_steps;
  MPI_Bcast(&nptiso.p_inst_av, 1, MPI_DOUBLE, 0, comm_cart);
}

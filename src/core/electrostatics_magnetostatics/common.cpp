/*
 * Copyright (C) 2010-2020 The ESPResSo project
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

#include "electrostatics_magnetostatics/common.hpp"

#include "electrostatics_magnetostatics/coulomb.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#include "communication.hpp"
#include "config.hpp"
#include "event.hpp"

#include <mpi.h>

static void mpi_bcast_coulomb_params_local() {
#ifdef ELECTROSTATICS
  MPI_Bcast(&coulomb, sizeof(Coulomb_parameters), MPI_BYTE, 0, comm_cart);
  Coulomb::bcast_coulomb_params();
#endif

#ifdef DIPOLES
  MPI_Bcast(&dipole, sizeof(Dipole_parameters), MPI_BYTE, 0, comm_cart);
  Dipole::set_method_local(dipole.method);
  Dipole::bcast_params(comm_cart);
#endif

#if defined(ELECTROSTATICS) || defined(DIPOLES)
  on_coulomb_change();
#endif
}

REGISTER_CALLBACK(mpi_bcast_coulomb_params_local)

void mpi_bcast_coulomb_params() {
#if defined(ELECTROSTATICS) || defined(DIPOLES)
  mpi_call_all(mpi_bcast_coulomb_params_local);
#endif
}

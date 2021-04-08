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
 *  Implementation of \ref debye_hueckel.hpp
 */

#include "debye_hueckel.hpp"

#ifdef ELECTROSTATICS
#include "electrostatics_magnetostatics/common.hpp"

#include <stdexcept>

Debye_hueckel_params dh_params{};

void dh_set_params(double kappa, double r_cut) {
  if (kappa < 0.0)
    throw std::domain_error("kappa should be a non-negative number");
  if (r_cut < 0.0)
    throw std::domain_error("r_cut should be a non-negative number");

  dh_params = {r_cut, kappa};

  mpi_bcast_coulomb_params();
}

#endif

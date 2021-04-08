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
 *  Implementation of \ref reaction_field.hpp
 */
#include "reaction_field.hpp"

#ifdef ELECTROSTATICS
#include "common.hpp"

Reaction_field_params rf_params{};

void rf_set_params(double kappa, double epsilon1, double epsilon2,
                   double r_cut) {
  if (kappa < 0.0)
    throw std::domain_error("kappa should be a non-negative number");
  if (epsilon1 < 0.0)
    throw std::domain_error("epsilon1 should be a non-negative number");
  if (epsilon2 < 0.0)
    throw std::domain_error("epsilon2 should be a non-negative number");
  if (r_cut < 0.0)
    throw std::domain_error("r_cut should be a non-negative number");

  auto const B = (2 * (epsilon1 - epsilon2) * (1 + kappa * r_cut) -
                  epsilon2 * kappa * kappa * r_cut * r_cut) /
                 ((epsilon1 + 2 * epsilon2) * (1 + kappa * r_cut) +
                  epsilon2 * kappa * kappa * r_cut * r_cut);
  rf_params = {kappa, epsilon1, epsilon2, r_cut, B};

  mpi_bcast_coulomb_params();
}
#endif

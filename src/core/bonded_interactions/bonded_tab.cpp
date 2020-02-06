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
#include "bonded_interactions/bonded_tab.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"

int tabulated_bonded_set_params(int bond_type,
                                TabulatedBondedInteraction tab_type, double min,
                                double max, std::vector<double> const &energy,
                                std::vector<double> const &force) {
  if (bond_type < 0)
    return ES_ERROR;

  assert(max >= min);
  assert((max == min) || force.size() > 1);
  assert(force.size() == energy.size());

  make_bond_type_exist(bond_type);

  /* set types */
  auto tab_pot = bonded_ia_params[bond_type].p.tab.pot = new TabulatedPotential;
  bonded_ia_params[bond_type].p.tab.type = tab_type;

  /* set number of interaction partners */
  switch (tab_type) {
  case TAB_BOND_LENGTH:
    tab_pot->minval = min;
    tab_pot->maxval = max;
    bonded_ia_params[bond_type].num = 1;
    bonded_ia_params[bond_type].type = BONDED_IA_TABULATED_DISTANCE;
    break;
  case TAB_BOND_ANGLE:
    tab_pot->minval = 0.0;
    tab_pot->maxval = Utils::pi() + ROUND_ERROR_PREC;
    bonded_ia_params[bond_type].num = 2;
    bonded_ia_params[bond_type].type = BONDED_IA_TABULATED_ANGLE;
    break;
  case TAB_BOND_DIHEDRAL:
    tab_pot->minval = 0.0;
    tab_pot->maxval = 2.0 * Utils::pi() + ROUND_ERROR_PREC;
    bonded_ia_params[bond_type].num = 3;
    bonded_ia_params[bond_type].type = BONDED_IA_TABULATED_DIHEDRAL;
    break;
  default:
    runtimeErrorMsg() << "Unsupported tabulated bond type.";
    return ES_ERROR;
  }

  tab_pot->invstepsize = static_cast<double>(force.size() - 1) / (max - min);

  tab_pot->force_tab = force;
  tab_pot->energy_tab = energy;

  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

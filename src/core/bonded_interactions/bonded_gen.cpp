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
#include "bonded_interactions/bonded_gen.hpp"

#include "communication.hpp"
#include "errorhandling.hpp"

int generic_bonded_set_params(int bond_type, GenericBondedInteraction type,
                              double max, std::string const &energy,
                              std::string const &force) {
  if (bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  /* set types */
  auto gen_pot = bonded_ia_params[bond_type].p.gen.pot = new GenericPotential;

  /* set number of interaction partners */
  switch (type) {
  case GEN_BOND_LENGTH:
    gen_pot->maxval = max;
    bonded_ia_params[bond_type].num = 1;
    bonded_ia_params[bond_type].type = BONDED_IA_GENERIC_DISTANCE;
    break;
  case GEN_BOND_ANGLE:
    gen_pot->maxval = Utils::pi() + ROUND_ERROR_PREC;
    bonded_ia_params[bond_type].num = 2;
    bonded_ia_params[bond_type].type = BONDED_IA_GENERIC_ANGLE;
    break;
  case GEN_BOND_DIHEDRAL:
    gen_pot->maxval = 2.0 * Utils::pi() + ROUND_ERROR_PREC;
    bonded_ia_params[bond_type].num = 3;
    bonded_ia_params[bond_type].type = BONDED_IA_GENERIC_DIHEDRAL;
    break;
  default:
    runtimeErrorMsg() << "Unsupported generic bond type.";
    return ES_ERROR;
  }

  gen_pot->force_expr = force;
  gen_pot->energy_expr = energy;

  gen_pot->force_parser = std::make_shared<matheval::Parser>();
  gen_pot->energy_parser = std::make_shared<matheval::Parser>();

  gen_pot->parse();

  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

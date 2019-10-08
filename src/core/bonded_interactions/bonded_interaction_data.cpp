/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
#include "bonded_interaction_data.hpp"
#include "communication.hpp"

#include <utils/constants.hpp>

std::vector<Bonded_ia_parameters> bonded_ia_params;

double recalc_maximal_cutoff_bonded() {
  auto max_cut_bonded = -1.;

  for (auto const &bonded_ia_param : bonded_ia_params) {
    switch (bonded_ia_param.type) {
    case BONDED_IA_FENE:
      max_cut_bonded =
          std::max(max_cut_bonded,
                   bonded_ia_param.p.fene.r0 + bonded_ia_param.p.fene.drmax);
      break;
    case BONDED_IA_HARMONIC:
      max_cut_bonded =
          std::max(max_cut_bonded, bonded_ia_param.p.harmonic.r_cut);
      break;
    case BONDED_IA_THERMALIZED_DIST:
      max_cut_bonded =
          std::max(max_cut_bonded, bonded_ia_param.p.thermalized_bond.r_cut);
      break;
    case BONDED_IA_RIGID_BOND:
      max_cut_bonded =
          std::max(max_cut_bonded, std::sqrt(bonded_ia_param.p.rigid_bond.d2));
      break;
    case BONDED_IA_TABULATED_DISTANCE:
      max_cut_bonded =
          std::max(max_cut_bonded, bonded_ia_param.p.tab.pot->cutoff());
      break;
    case BONDED_IA_IBM_TRIEL:
      max_cut_bonded =
          std::max(max_cut_bonded, bonded_ia_param.p.ibm_triel.maxDist);
      break;
    default:
      break;
    }
  }

  /* dihedrals: the central particle is indirectly connected to the fourth
   * particle via the third particle, so we have to double the cutoff */
  for (auto const &bonded_ia_param : bonded_ia_params) {
    switch (bonded_ia_param.type) {
    case BONDED_IA_DIHEDRAL:
    case BONDED_IA_TABULATED_DIHEDRAL:
      max_cut_bonded *= 2;
      break;
    default:
      break;
    }
  }

  return max_cut_bonded;
}

void make_bond_type_exist(int type) {
  int i, ns = type + 1;
  const auto old_size = bonded_ia_params.size();
  if (ns <= bonded_ia_params.size()) {
    return;
  }
  /* else allocate new memory */
  bonded_ia_params.resize(ns);
  /* set bond types not used as undefined */
  for (i = old_size; i < ns; i++)
    bonded_ia_params[i].type = BONDED_IA_NONE;
}

int virtual_set_params(int bond_type) {
  if (bond_type < 0)
    return ES_ERROR;

  make_bond_type_exist(bond_type);

  bonded_ia_params[bond_type].type = BONDED_IA_VIRTUAL_BOND;
  bonded_ia_params[bond_type].num = 1;

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

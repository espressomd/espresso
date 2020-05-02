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

#include <boost/algorithm/cxx11/any_of.hpp>
#include <boost/range/algorithm/copy.hpp>
#include <boost/range/numeric.hpp>

#include <boost/range/algorithm/find.hpp>
#include <errorhandling.hpp>

#include <utils/constants.hpp>

std::vector<Bonded_ia_parameters> bonded_ia_params;

auto cutoff(int type, Bond_parameters const &bp) {
  switch (type) {
  case BONDED_IA_NONE:
    return -1.;
  case BONDED_IA_FENE:
    return bp.fene.cutoff();
  case BONDED_IA_HARMONIC:
    return bp.harmonic.cutoff();
  case BONDED_IA_HARMONIC_DUMBBELL:
    return bp.harmonic_dumbbell.cutoff();
  case BONDED_IA_QUARTIC:
    return bp.quartic.cutoff();
  case BONDED_IA_BONDED_COULOMB:
    return bp.bonded_coulomb.cutoff();
  case BONDED_IA_BONDED_COULOMB_SR:
    return bp.bonded_coulomb_sr.cutoff();
  case BONDED_IA_DIHEDRAL:
    return bp.dihedral.cutoff();
  case BONDED_IA_TABULATED_DISTANCE:
  case BONDED_IA_TABULATED_ANGLE:
  case BONDED_IA_TABULATED_DIHEDRAL:
    return bp.tab.cutoff();
  case BONDED_IA_RIGID_BOND:
    return bp.rigid_bond.cutoff();
  case BONDED_IA_VIRTUAL_BOND:
    return bp.virt.cutoff();
  case BONDED_IA_ANGLE_HARMONIC:
    return bp.angle_harmonic.cutoff();
  case BONDED_IA_ANGLE_COSINE:
    return bp.angle_cosine.cutoff();
  case BONDED_IA_ANGLE_COSSQUARE:
    return bp.angle_cossquare.cutoff();
  case BONDED_IA_OIF_LOCAL_FORCES:
    return bp.oif_local_forces.cutoff();
  case BONDED_IA_OIF_GLOBAL_FORCES:
    return bp.oif_global_forces.cutoff();
  case BONDED_IA_IBM_TRIEL:
    return bp.ibm_triel.cutoff();
  case BONDED_IA_IBM_VOLUME_CONSERVATION:
    return bp.ibmVolConsParameters.cutoff();
  case BONDED_IA_IBM_TRIBEND:
    return bp.ibm_tribend.cutoff();
  case BONDED_IA_UMBRELLA:
    return bp.umbrella.cutoff();
  case BONDED_IA_THERMALIZED_DIST:
    return bp.thermalized_bond.cutoff();
  default:
    throw std::runtime_error("Unknown bond type.");
  }
}

double maximal_cutoff_bonded() {
  auto const max_cut_bonded =
      boost::accumulate(bonded_ia_params, -1.,
                        [](auto max_cut, Bonded_ia_parameters const &bond) {
                          return std::max(max_cut, cutoff(bond.type, bond.p));
                        });

  /* Check if there are dihedrals */
  auto const any_dihedrals = std::any_of(
      bonded_ia_params.begin(), bonded_ia_params.end(), [](auto const &bond) {
        switch (bond.type) {
        case BONDED_IA_DIHEDRAL:
        case BONDED_IA_TABULATED_DIHEDRAL:
          return true;
        default:
          return false;
        }
      });

  /* dihedrals: the central particle is indirectly connected to the fourth
   * particle via the third particle, so we have to double the cutoff */
  return (any_dihedrals) ? 2 * max_cut_bonded : max_cut_bonded;
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
  bonded_ia_params[bond_type].p.virt = VirtualBond_Parameters{};

  /* broadcast interaction parameters */
  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

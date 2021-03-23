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
#ifndef _BONDED_INTERACTION_DATA_HPP
#define _BONDED_INTERACTION_DATA_HPP
/** @file
 *  Data structures for bonded interactions.
 *  For more information on how to add new interactions, see @ref bondedIA_new.
 */

#include "angle_common.hpp"
#include "angle_cosine.hpp"
#include "angle_cossquare.hpp"
#include "angle_harmonic.hpp"
#include "bonded_coulomb.hpp"
#include "bonded_coulomb_sr.hpp"
#include "bonded_tab.hpp"
#include "dihedral.hpp"
#include "fene.hpp"
#include "harmonic.hpp"
#include "immersed_boundary/ibm_tribend.hpp"
#include "immersed_boundary/ibm_triel.hpp"
#include "immersed_boundary/ibm_volcons.hpp"
#include "object-in-fluid/oif_global_forces_params.hpp"
#include "object-in-fluid/oif_local_forces.hpp"
#include "quartic.hpp"
#include "rattle.hpp"
#include "rigid_bond.hpp"
#include "thermalized_bond.hpp"

#include "TabulatedPotential.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/variant.hpp>

#include <cassert>
#include <cmath>
#include <stdexcept>
#include <vector>

/** Interaction type for unused bonded interaction slots */
struct None_bond_parameters {
  static constexpr int num = 0;
  double cutoff() const { return -1.; }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {}
};

/** Interaction type for virtual bonds */
struct VirtualBond_Parameters {
  static constexpr int num = 1;
  double cutoff() const { return -1.; }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {}
};

/** Visitor to get the number of bound partners from the bond parameter
 *  variant.
 */
class BondNumPartners : public boost::static_visitor<int> {
public:
  template <typename T> int operator()(T const &) const { return T::num; }
};

/** Variant in which to store the parameters of an individual bonded
 *  interaction
 */
using Bonded_ia_parameters = boost::variant<
    None_bond_parameters, Fene_bond_parameters, Harmonic_bond_parameters,
    Quartic_bond_parameters, Bonded_coulomb_bond_parameters,
    Bonded_coulomb_sr_bond_parameters, Angle_harmonic_bond_parameters,
    Angle_cosine_bond_parameters, Angle_cossquare_bond_parameters,
    Dihedral_bond_parameters, Tabulated_distance_bond_parameters,
    Tabulated_angle_bond_parameters, Tabulated_dihedral_bond_parameters,
    Thermalized_bond_parameters, Rigid_bond_parameters, IBM_Triel_Parameters,
    IBM_VolCons_Parameters, IBM_Tribend_Parameters,
    Oif_global_forces_bond_parameters, Oif_local_forces_bond_parameters,
    VirtualBond_Parameters>;

/** Field containing the parameters of the bonded ia types */
extern std::vector<Bonded_ia_parameters> bonded_ia_params;

/** Makes sure that \ref bonded_ia_params is large enough to cover the
 *  parameters for the bonded interaction type.
 *  Attention: 1: There is no initialization done here.
 *  2: Use only in connection with creating new or overwriting old bond types
 */
void make_bond_type_exist(int type);

/** Calculate the maximal cutoff of bonded interactions, required to
 *  determine the cell size for communication.
 *
 *  Bond angle and dihedral potentials do not contain a cutoff intrinsically.
 *  The cutoff for these potentials depends on the bond length potentials
 *  (it is assumed that particles participating in a bond angle or dihedral
 *  potential are bound to each other by some bond length potential). For bond
 *  angle potentials nothing has to be done. For dihedral potentials the cutoff
 *  is set to twice the maximal cutoff because the particle in which the bond
 *  is stored is only bonded to the first two partners, one of which has an
 *  additional bond to the third partner.
 */
double maximal_cutoff_bonded();

/** Return the number of bonded partners for the specified bond */
inline int number_of_partners(Bonded_ia_parameters const &iaparams) {
  return boost::apply_visitor(BondNumPartners(), iaparams);
}

#endif

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
#include "rigid_bond.hpp"
#include "thermalized_bond.hpp"

#include "TabulatedPotential.hpp"

#include <boost/serialization/access.hpp>
#include <boost/serialization/variant.hpp>
#include <boost/variant.hpp>

#include <algorithm>
#include <cassert>
#include <cmath>
#include <stdexcept>
#include <unordered_map>
#include <vector>

/* Special cutoff value for a disabled bond.
 * Bonds that have this cutoff are not visited during bond evaluation.
 */
static constexpr double BONDED_INACTIVE_CUTOFF = -1.;

/** Interaction type for unused bonded interaction slots */
struct NoneBond {
  static constexpr int num = 0;
  double cutoff() const { return BONDED_INACTIVE_CUTOFF; }

private:
  friend boost::serialization::access;
  template <typename Archive>
  void serialize(Archive &ar, long int /* version */) {}
};

/** Interaction type for virtual bonds */
struct VirtualBond {
  static constexpr int num = 1;
  double cutoff() const { return BONDED_INACTIVE_CUTOFF; }

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
using Bonded_IA_Parameters =
    boost::variant<NoneBond, FeneBond, HarmonicBond, QuarticBond, BondedCoulomb,
                   BondedCoulombSR, AngleHarmonicBond, AngleCosineBond,
                   AngleCossquareBond, DihedralBond, TabulatedDistanceBond,
                   TabulatedAngleBond, TabulatedDihedralBond, ThermalizedBond,
                   RigidBond, IBMTriel, IBMVolCons, IBMTribend,
                   OifGlobalForcesBond, OifLocalForcesBond, VirtualBond>;

class BondedInteractionsMap {
  using container_type =
      std::unordered_map<int, std::shared_ptr<Bonded_IA_Parameters>>;

public:
  using key_type = typename container_type::key_type;
  using mapped_type = typename container_type::mapped_type;
  using value_type = typename container_type::value_type;
  using iterator = typename container_type::iterator;
  using const_iterator = typename container_type::const_iterator;

  explicit BondedInteractionsMap() = default;

  iterator begin() { return m_params.begin(); }
  iterator end() { return m_params.end(); }
  const_iterator begin() const { return m_params.begin(); }
  const_iterator end() const { return m_params.end(); }

  void insert(key_type const &key, mapped_type const &ptr) {
    next_key = std::max(next_key, key + 1);
    m_params[key] = ptr;
  }
  key_type insert(mapped_type const &ptr) {
    auto const key = next_key++;
    m_params[key] = ptr;
    return key;
  }
  auto erase(key_type const &key) { return m_params.erase(key); }
  mapped_type at(key_type const &key) const { return m_params.at(key); }
  auto count(key_type const &key) const { return m_params.count(key); }
  bool contains(key_type const &key) const { return m_params.count(key); }
  bool empty() const { return m_params.empty(); }
  auto size() const { return m_params.size(); }
  auto get_next_key() const { return next_key; }

private:
  container_type m_params = {};
  key_type next_key = static_cast<key_type>(0);
};

/** Notify the cell system about changes to the maximal interaction range. */
void mpi_update_cell_system_ia_range_local();

/** Field containing the parameters of the bonded ia types */
extern BondedInteractionsMap bonded_ia_params;

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
inline int number_of_partners(Bonded_IA_Parameters const &iaparams) {
  return boost::apply_visitor(BondNumPartners(), iaparams);
}

#endif

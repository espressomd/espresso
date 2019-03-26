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
#ifndef PAIR_CRITERIA_HPP
#define PAIR_CRITERIA_HPP

#include "energy_inline.hpp"
#include "particle_data.hpp"
#include <stdexcept>

namespace PairCriteria {

/** @brief Criterion which provides a true/false for a pair of particles */
class PairCriterion {
public:
  /** @brief Make a decision based on two Particle objects */
  virtual bool decide(const Particle &p1, const Particle &p2) const = 0;
  /** @brief Make a decision based on particle ids.
   * This can only run on the master node outside the integration loop */
  bool decide(int id1, int id2) const {
    // Retrieve particle data
    auto const &p1 = get_particle_data(id1);
    auto const &p2 = get_particle_data(id2);
    const bool res = decide(p1, p2);
    return res;
  }
  virtual ~PairCriterion() = default;
};

/** @brief True if two particles are closer than a cut off distance, respecting
 * minimum image convention */
class DistanceCriterion : public PairCriterion {
public:
  bool decide(const Particle &p1, const Particle &p2) const override {
    double vec21[3];
    get_mi_vector(vec21, p1.r.p, p2.r.p);
    return sqrt(sqrlen(vec21)) <= m_cut_off;
  };
  double get_cut_off() { return m_cut_off; }
  void set_cut_off(double c) { m_cut_off = c; }

private:
  double m_cut_off;
};

/** True if the short range energy is larger than a cut_off */
class EnergyCriterion : public PairCriterion {
public:
  bool decide(const Particle &p1, const Particle &p2) const override {
    // Distance between particles
    double vec21[3];
    get_mi_vector(vec21, p1.r.p, p2.r.p);
    const double dist_betw_part = sqrt(sqrlen(vec21));

    // Interaction parameters for particle types
    IA_parameters *ia_params = get_ia_param(p1.p.type, p2.p.type);

    return (calc_non_bonded_pair_energy(
               &p1, &p2, ia_params, vec21, dist_betw_part,
               dist_betw_part * dist_betw_part)) >= m_cut_off;
  };
  double get_cut_off() { return m_cut_off; }
  void set_cut_off(double c) { m_cut_off = c; }

private:
  double m_cut_off;
};

/** True if a bond of given type exists between the two particles */
class BondCriterion : public PairCriterion {
public:
  bool decide(const Particle &p1, const Particle &p2) const override {
    return pair_bond_exists_on(&p1, &p2, m_bond_type) ||
           pair_bond_exists_on(&p2, &p1, m_bond_type);
  };
  int get_bond_type() { return m_bond_type; };
  void set_bond_type(int t) { m_bond_type = t; }

private:
  int m_bond_type;
};
} // namespace PairCriteria

#endif

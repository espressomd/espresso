/*
 * Copyright (C) 2010-2022 The ESPResSo project
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

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_SCAFACOS_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_SCAFACOS_HPP

#include "config/config.hpp"

#ifdef SCAFACOS

#include "electrostatics/actor.hpp"

#include "scafacos/ScafacosContextBase.hpp"

#include <utils/Vector.hpp>

#include <memory>
#include <string>

struct CoulombScafacos : virtual public ScafacosContextBase,
                         public Coulomb::Actor<CoulombScafacos> {
  ~CoulombScafacos() override = default;

  void on_activation() {
    update_system_params();
    sanity_checks();
    tune();
  }
  /** @brief Recalculate all box-length-dependent parameters. */
  void on_boxl_change() { update_system_params(); }
  void on_node_grid_change() const {}
  void on_periodicity_change() { update_system_params(); }
  void on_cell_structure_change() const {}
  void init() const {}

  void sanity_checks() const override { sanity_checks_charge_neutrality(); }

  bool is_tuned() const { return m_is_tuned; }
  void tune() {
    if (not is_tuned()) {
      tune_impl();
    }
  }
  virtual double get_r_cut() const = 0;
  virtual double get_pair_force(double dist) const = 0;
  virtual double get_pair_energy(double dist) const = 0;
  virtual void set_near_field_delegation(bool delegate) = 0;
  virtual bool get_near_field_delegation() const = 0;

  /** @brief Calculate near-field pair force. */
  Utils::Vector3d pair_force(double q1q2, Utils::Vector3d const &d,
                             double dist) const {
    if (dist > get_r_cut())
      return {};

    return d * (-q1q2 * prefactor * get_pair_force(dist) / dist);
  }

  /** @brief Calculate near-field pair energy. */
  double pair_energy(double q1q2, double dist) const {
    if (dist > get_r_cut())
      return 0.;

    return q1q2 * prefactor * get_pair_energy(dist);
  }

protected:
  virtual void tune_impl() = 0;

private:
  bool m_is_tuned = false;
};

std::shared_ptr<CoulombScafacos>
make_coulomb_scafacos(std::string const &method, std::string const &parameters);

#endif // SCAFACOS
#endif

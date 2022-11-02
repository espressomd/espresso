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

#ifndef ESPRESSO_SRC_CORE_MAGNETOSTATICS_DLC_HPP
#define ESPRESSO_SRC_CORE_MAGNETOSTATICS_DLC_HPP

#include "config/config.hpp"

#ifdef DIPOLES

#include "actor/traits.hpp"

#include "magnetostatics/dds.hpp"
#include "magnetostatics/dp3m.hpp"

#include <ParticleRange.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <memory>

struct DipolarLayerCorrection;

namespace traits {
template <>
struct is_layer_correction<DipolarLayerCorrection> : std::true_type {};
} // namespace traits

/** @brief Parameters for the DLC method */
struct dlc_data {
  dlc_data(double maxPWerror, double gap_size, double far_cut);

  /** maximal pairwise error of the potential and force */
  double maxPWerror;

  /** Size of the empty gap. Note that MDLC relies on the user to make sure
   *  that this condition is fulfilled.
   */
  double gap_size;

  /** Up to where particles can be found */
  double box_h;

  /** Cutoff of the exponential sum. Since in all other MMM methods this is
   *  the far formula, we call it here the same, although in the ELC context
   *  it does not make much sense.
   */
  double far_cut;

  /** Flag whether #far_cut was set by the user, or calculated by ESPResSo.
   *  In the latter case, the cutoff will be adapted if important parameters,
   *  such as the box dimensions, change.
   */
  bool far_calculated;
};

/**
 * @brief Adapt a magnetostatics solver to remove contributions from the
 * z-direction. For details see @cite brodka04a.
 */
struct DipolarLayerCorrection {
  using BaseSolver = boost::variant<
#ifdef DP3M
      std::shared_ptr<DipolarP3M>,
#endif
      std::shared_ptr<DipolarDirectSum>>;

  /** @name Variables from the adapted solver. */
  /**@{*/
  double prefactor;
  double epsilon;
  double epsilon_correction;
  /**@}*/
  dlc_data dlc;

  /** Magnetostatics solver that is adapted. */
  BaseSolver base_solver;

  DipolarLayerCorrection(dlc_data &&parameters, BaseSolver &&solver);

  void on_activation() {
    sanity_checks_node_grid();
    /* None of the DLC parameters depend on the DP3M parameters,
     * but the DP3M parameters depend on the DLC parameters during tuning,
     * therefore DLC needs to be tuned before DP3M. */
    recalc_box_h();
    recalc_far_cut();
    visit_base_solver([](auto &solver) { solver->on_activation(); });
  }
  /** @brief Recalculate all box-length-dependent parameters. */
  void on_boxl_change() {
    recalc_box_h();
    recalc_far_cut();
    visit_base_solver([](auto &actor) { actor->on_boxl_change(); });
  }
  void on_node_grid_change() const {
    sanity_checks_node_grid();
    visit_base_solver([](auto &solver) { solver->on_node_grid_change(); });
  }
  void on_periodicity_change() const {
    visit_base_solver([](auto &solver) { solver->on_periodicity_change(); });
  }
  void on_cell_structure_change() const {
    visit_base_solver([](auto &solver) { solver->on_cell_structure_change(); });
  }
  void init() {
    recalc_box_h();
    recalc_far_cut();
    visit_base_solver([](auto &solver) { solver->init(); });
  }

  void sanity_checks() const {
    sanity_checks_node_grid();
    visit_base_solver([](auto &actor) { actor->sanity_checks(); });
  }

  void recalc_box_h();
  void recalc_far_cut() {
    if (dlc.far_calculated) {
      dlc.far_cut = tune_far_cut();
    }
  }

  /** @brief Calculate the dipolar energy correction. */
  double energy_correction(ParticleRange const &particles) const;
  /** @brief Add the dipolar force and torque corrections. */
  void add_force_corrections(ParticleRange const &particles) const;

  void adapt_solver();

  void release_solver() {
    prefactor = -1.;
    epsilon = -1.;
    epsilon_correction = 0.;
  }

private:
  /** Check if a magnetic particle is in the gap region. */
  void check_gap(Particle const &p) const;
  double tune_far_cut() const;

  void sanity_checks_node_grid() const;

  template <class Visitor> void visit_base_solver(Visitor &&visitor) const {
    boost::apply_visitor(visitor, base_solver);
  }
};

#endif // DIPOLES

#endif

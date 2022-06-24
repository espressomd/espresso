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

/** @file
 *  @brief ELC algorithm for long-range Coulomb interactions.
 *
 *  Implementation of the ELC method for the calculation of the electrostatic
 *  interaction in two dimensional periodic systems. For details on the method
 *  see MMM in general. The ELC method works together with any three-dimensional
 *  method, for example @ref p3m.hpp "P3M", with metallic boundary conditions.
 */

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_ELC_HPP
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_ELC_HPP

#include "config.hpp"

#ifdef P3M

#include "actor/traits.hpp"

#include "electrostatics/p3m.hpp"
#include "electrostatics/p3m_gpu.hpp"

#include "Particle.hpp"
#include "ParticleRange.hpp"

#include <utils/Vector.hpp>

#include <boost/optional.hpp>
#include <boost/variant.hpp>

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>

struct ElectrostaticLayerCorrection;

namespace traits {
template <>
struct is_layer_correction<ElectrostaticLayerCorrection> : std::true_type {};
} // namespace traits

/** @brief Parameters for the ELC method */
struct elc_data {
  elc_data(double maxPWerror, double gap_size, double far_cut, bool neutralize,
           double delta_top, double delta_bot, bool const_pot, double pot_diff);

  /**
   * @brief Maximal allowed pairwise error for the potential and force.
   * Note that this counts for the plain 1/r contribution
   * alone, without the prefactor and the charge prefactor.
   */
  double maxPWerror;
  /** Size of the empty gap. Note that ELC relies on the user to make sure
   *  that this condition is fulfilled.
   */
  double gap_size;
  /** Up to where particles can be found. */
  double box_h;
  /**
   * @brief Cutoff of the exponential sum.
   * Since in all other MMM methods this is the far formula,
   * it is given the same name here.
   */
  double far_cut;
  /** Squared value of #far_cut. */
  double far_cut2;
  /** Flag whether #far_cut was set by the user, or calculated by ESPResSo.
   *  In the latter case, the cutoff will be adapted if important parameters,
   *  such as the box dimensions, change.
   */
  bool far_calculated;

  /// @brief Flag whether there is any dielectric contrast in the system.
  bool dielectric_contrast_on;
  /// @brief Flag whether a constant potential difference is applied.
  bool const_pot;
  /**
   * @brief Flag whether the box is neutralized by a homogeneous background.
   * If true, use a homogeneous neutralizing background for non-neutral
   * systems. Unlike the 3D case, this background adds an additional
   * force pointing towards the system center, and the gap size
   * enters into the value of the forces, so be careful with this.
   */
  bool neutralize;

  /// dielectric contrast in the upper part of the simulation cell.
  double delta_mid_top;
  /// dielectric contrast in the lower part of the simulation cell.
  double delta_mid_bot;
  /// @brief Constant potential difference.
  double pot_diff;

  /** Layer around the dielectric contrast in which we trick around. */
  double space_layer;
  /** The space that is finally left. */
  double space_box;

  /// pairwise contributions from the lowest and top layers
  template <typename Kernel>
  void dielectric_layers_contribution(CoulombP3M const &p3m,
                                      Utils::Vector3d const &pos1,
                                      Utils::Vector3d const &pos2, double q1q2,
                                      Kernel &&kernel) const {
    if (pos1[2] < space_layer) {
      auto const q_eff = delta_mid_bot * q1q2;
      auto const d = get_mi_vector(pos2, {pos1[0], pos1[1], -pos1[2]});
      kernel(q_eff, d);
    }
    if (pos1[2] > (box_h - space_layer)) {
      auto const q_eff = delta_mid_top * q1q2;
      auto const l = 2. * box_h;
      auto const d = get_mi_vector(pos2, {pos1[0], pos1[1], l - pos1[2]});
      kernel(q_eff, d);
    }
  }

  /// self energies of top and bottom layers with their virtual images
  double dielectric_layers_self_energy(CoulombP3M const &p3m,
                                       ParticleRange const &particles) const {
    auto energy = 0.;
    for (auto const &p : particles) {
      dielectric_layers_contribution(
          p3m, p.pos(), p.pos(), Utils::sqr(p.q()),
          [&](double q1q2, Utils::Vector3d const &d) {
            energy += p3m.pair_energy(q1q2, d.norm());
          });
    }
    return energy;
  }

  /// forces of particles in border layers with themselves
  void dielectric_layers_self_forces(CoulombP3M const &p3m,
                                     ParticleRange const &particles) const {
    for (auto &p : particles) {
      dielectric_layers_contribution(
          p3m, p.pos(), p.pos(), Utils::sqr(p.q()),
          [&](double q1q2, Utils::Vector3d const &d) {
            p.force() += p3m.pair_force(q1q2, d, d.norm());
          });
    }
  }

private:
  Utils::Vector3d get_mi_vector(Utils::Vector3d const &a,
                                Utils::Vector3d const &b) const;
};

struct ElectrostaticLayerCorrection
    : public Coulomb::Actor<ElectrostaticLayerCorrection> {
  using BaseSolver = boost::variant<
#ifdef CUDA
      std::shared_ptr<CoulombP3MGPU>,
#endif // CUDA
      std::shared_ptr<CoulombP3M>>;

  elc_data elc;

  /** Electrostatics solver that is adapted. */
  BaseSolver base_solver;

  ElectrostaticLayerCorrection(elc_data &&parameters, BaseSolver &&solver);

  void on_activation() {
    sanity_checks_periodicity();
    sanity_checks_cell_structure();
    sanity_checks_charge_neutrality();
    sanity_checks_dielectric_contrasts();
    /* Most ELC parameters do not depend on the P3M parameters,
     * but the P3M parameters depend on the ELC parameters during tuning,
     * therefore ELC needs to be tuned before P3M. */
    recalc_box_h();
    recalc_far_cut();
    visit_base_solver([](auto &solver) { solver->on_activation(); });
    /* With dielectric contrasts, the ELC space layer depends
     * on the P3M real-space cutoff, and the P3M FFT parameters
     * depend on the ELC space layer */
    if (elc.dielectric_contrast_on) {
      recalc_space_layer();
      visit_base_solver([](auto &actor) { actor->init(); });
    }
  }
  /** @brief Recalculate all box-length-dependent parameters. */
  void on_boxl_change() {
    visit_base_solver([](auto &actor) { actor->on_boxl_change(); });
    recalc_box_h();
    recalc_far_cut();
    recalc_space_layer();
  }
  void on_node_grid_change() const {
    visit_base_solver([](auto &solver) { solver->on_node_grid_change(); });
  }
  void on_periodicity_change() const {
    sanity_checks_periodicity();
    visit_base_solver([](auto &solver) { solver->on_periodicity_change(); });
  }
  void on_cell_structure_change() {
    sanity_checks_cell_structure();
    visit_base_solver([](auto &solver) { solver->on_cell_structure_change(); });
    recalc_box_h();
    recalc_far_cut();
    if (elc.dielectric_contrast_on) {
      recalc_space_layer();
      visit_base_solver([](auto &actor) { actor->init(); });
    }
  }
  /** @brief Recalculate all derived parameters. */
  void init() {
    visit_base_solver([](auto &actor) { actor->init(); });
    recalc_box_h();
    recalc_far_cut();
    recalc_space_layer();
  }

  void sanity_checks() const {
    sanity_checks_periodicity();
    sanity_checks_cell_structure();
    sanity_checks_charge_neutrality();
    sanity_checks_dielectric_contrasts();
    visit_base_solver([](auto &actor) { actor->sanity_checks(); });
  }

  /**
   * @brief Veto real-space cutoff values that are incompatible with ELC.
   * When ELC is used with dielectric contrasts, the short-range cutoff needs
   * to be smaller than the gap size to allow placement of the image charges.
   */
  boost::optional<std::string> veto_r_cut(double r_cut) const {
    if (elc.dielectric_contrast_on and r_cut >= elc.gap_size) {
      return {std::string("conflict with ELC w/ dielectric contrasts")};
    }
    return {};
  }

  /** @brief Calculate short-range pair energy correction. */
  double pair_energy_correction(double q1q2, Particle const &p1,
                                Particle const &p2) const {
    double energy = 0.;
    if (elc.dielectric_contrast_on) {
      energy = boost::apply_visitor(
          [this, &p1, &p2, q1q2](auto &p3m_ptr) {
            auto const &pos1 = p1.pos();
            auto const &pos2 = p2.pos();
            auto const &p3m = *p3m_ptr;
            auto energy = 0.;
            elc.dielectric_layers_contribution(
                p3m, pos1, pos2, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  energy += p3m.pair_energy(q_eff, d.norm());
                });
            elc.dielectric_layers_contribution(
                p3m, pos2, pos1, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  energy += p3m.pair_energy(q_eff, d.norm());
                });
            return energy / 2.;
          },
          base_solver);
    }
    return energy;
  }

  /** @brief Add short-range pair force corrections. */
  void add_pair_force_corrections(Particle &p1, Particle &p2,
                                  double q1q2) const {
    if (elc.dielectric_contrast_on) {
      boost::apply_visitor(
          [this, &p1, &p2, q1q2](auto &p3m_ptr) {
            auto const &pos1 = p1.pos();
            auto const &pos2 = p2.pos();
            auto const &p3m = *p3m_ptr;
            elc.dielectric_layers_contribution(
                p3m, pos1, pos2, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  p1.force() += p3m.pair_force(q_eff, d, d.norm());
                });
            elc.dielectric_layers_contribution(
                p3m, pos2, pos1, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  p2.force() += p3m.pair_force(q_eff, d, d.norm());
                });
          },
          base_solver);
    }
  }

  double long_range_energy(ParticleRange const &particles) const;
  void add_long_range_forces(ParticleRange const &particles) const;

private:
  /** Check if a charged particle is in the gap region. */
  void check_gap(Particle const &p) const;
  double tune_far_cut() const;
  void adapt_solver();
  /** pairwise contributions from the lowest and top layers to the energy */
  double dipole_energy(ParticleRange const &particles) const;
  void add_dipole_force(ParticleRange const &particles) const;
  double z_energy(ParticleRange const &particles) const;
  void add_z_force(ParticleRange const &particles) const;

  void recalc_box_h();
  void recalc_far_cut() {
    if (elc.far_calculated) {
      elc.far_cut = tune_far_cut();
    }
    elc.far_cut2 = Utils::sqr(elc.far_cut);
  }
  void recalc_space_layer();

  void sanity_checks_cell_structure() const {}
  void sanity_checks_periodicity() const;
  void sanity_checks_dielectric_contrasts() const;

  /// the force calculation
  void add_force(ParticleRange const &particles) const;
  /// the energy calculation
  double calc_energy(ParticleRange const &particles) const;

  template <class Visitor> void visit_base_solver(Visitor &&visitor) const {
    boost::apply_visitor(visitor, base_solver);
  }
};

#endif // P3M

#endif

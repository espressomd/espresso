/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

#include <config/config.hpp>

#ifdef ESPRESSO_P3M

#include "actor/traits.hpp"

#include "electrostatics/p3m.hpp"

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "ParticleRange.hpp"
#include "aosoa_pack.hpp"

#include <utils/Vector.hpp>
#include <utils/math/sqr.hpp>

#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <type_traits>
#include <utility>
#include <variant>

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

  /// pairwise contributions from lower and upper layers
  void dielectric_layers_contribution(BoxGeometry const &box_geo,
                                      std::size_t p1, std::size_t p2,
                                      auto &aosoa, double q1q2,
                                      auto &&kernel) const {
    if (aosoa.position(p1, 2) < space_layer) {
      auto const q_eff = delta_mid_bot * q1q2;
      auto pos2 = aosoa.get_vector_at(aosoa.position, p2);
      auto pos1 = aosoa.get_vector_at(aosoa.position, p1);
      pos1[2] *= -1.;
      auto const d = box_geo.get_mi_vector(pos2, pos1);
      kernel(q_eff, d);
    }
    if (aosoa.position(p1, 2) > (box_h - space_layer)) {
      auto const q_eff = delta_mid_top * q1q2;
      auto const z = 2. * box_h - aosoa.position(p1, 2);
      auto pos2 = aosoa.get_vector_at(aosoa.position, p2);
      auto pos1 = aosoa.get_vector_at(aosoa.position, p1);
      pos1[2] = 2. * box_h - pos1[2];
      auto const d = box_geo.get_mi_vector(pos2, pos1);
      kernel(q_eff, d);
    }
  }

  /// pairwise contributions from lower and upper layers
  void dielectric_layers_contribution(BoxGeometry const &box_geo,
                                      Utils::Vector3d const &pos1,
                                      Utils::Vector3d const &pos2, double q1q2,
                                      auto &&kernel) const {
    if (pos1[2] < space_layer) {
      auto const q_eff = delta_mid_bot * q1q2;
      auto const d = box_geo.get_mi_vector(pos2, {pos1[0], pos1[1], -pos1[2]});
      kernel(q_eff, d);
    }
    if (pos1[2] > (box_h - space_layer)) {
      auto const q_eff = delta_mid_top * q1q2;
      auto const z = 2. * box_h - pos1[2];
      auto const d = box_geo.get_mi_vector(pos2, {pos1[0], pos1[1], z});
      kernel(q_eff, d);
    }
  }

  /// self energies of top and bottom layers with their virtual images
  double dielectric_layers_self_energy(CoulombP3M const &p3m,
                                       BoxGeometry const &box_geo,
                                       ParticleRange const &particles) const {
    auto energy = 0.;
    for (auto const &p : particles) {
      dielectric_layers_contribution(
          box_geo, p.pos(), p.pos(), Utils::sqr(p.q()),
          [&](double q1q2, Utils::Vector3d const &d) {
            energy += p3m.pair_energy(q1q2, d.norm());
          });
    }
    return energy;
  }

  /// forces of particles in border layers with themselves
  void dielectric_layers_self_forces(CoulombP3M const &p3m,
                                     BoxGeometry const &box_geo,
                                     ParticleRange const &particles) const {
    for (auto &p : particles) {
      dielectric_layers_contribution(
          box_geo, p.pos(), p.pos(), Utils::sqr(p.q()),
          [&](double q1q2, Utils::Vector3d const &d) {
            p.force() += p3m.pair_force(q1q2, d, d.norm());
          });
    }
  }
};

struct ElectrostaticLayerCorrection
    : public Coulomb::Actor<ElectrostaticLayerCorrection> {
  using BaseSolver = std::variant<std::shared_ptr<CoulombP3M>>;

  elc_data elc;
  BoxGeometry *m_box_geo;

  /** Electrostatics solver that is adapted. */
  BaseSolver base_solver;

  ElectrostaticLayerCorrection(elc_data &&parameters, BaseSolver &&solver);

  void on_activation() {
    visit_base_solver(
        [this](auto &solver) { solver->bind_system(m_system.lock()); });
    sanity_checks_periodicity();
    sanity_checks_cell_structure();
    sanity_checks_charge_neutrality();
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
    visit_base_solver([](auto &actor) { actor->sanity_checks(); });
  }

  /**
   * @brief Veto real-space cutoff values that are incompatible with ELC.
   * When ELC is used with dielectric contrasts, the short-range cutoff needs
   * to be smaller than the gap size to allow placement of the image charges.
   */
  std::optional<std::string> veto_r_cut(double r_cut) const {
    if (elc.dielectric_contrast_on and r_cut >= elc.gap_size) {
      return {std::string("conflict with ELC w/ dielectric contrasts")};
    }
    return {};
  }

  Utils::Vector3d pair_force(double q1q2, Utils::Vector3d const &d,
                             double dist) const {
    return std::visit(
        [&](auto &solver) { return solver->pair_force(q1q2, d, dist); },
        base_solver);
  }

  /** @brief Calculate short-range pair energy correction. */
  double pair_energy_correction(std::size_t p1, std::size_t p2, auto &aosoa,
                                double q1q2) const {
    double energy = 0.;
    if (elc.dielectric_contrast_on) {
      energy = std::visit(
          [this, &aosoa, p1, p2, q1q2](auto &p3m_ptr) {
            auto const &p3m = *p3m_ptr;
            auto energy = 0.;
            elc.dielectric_layers_contribution(
                *m_box_geo, p1, p2, aosoa, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  energy += p3m.pair_energy(q_eff, d.norm());
                });
            elc.dielectric_layers_contribution(
                *m_box_geo, p2, p1, aosoa, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  energy += p3m.pair_energy(q_eff, d.norm());
                });
            return energy / 2.;
          },
          base_solver);
    }
    return energy;
  }

  /** @brief Calculate short-range pair energy correction. */
  double pair_energy_correction(Utils::Vector3d const &pos1,
                                Utils::Vector3d const &pos2,
                                double q1q2) const {
    double energy = 0.;
    if (elc.dielectric_contrast_on) {
      energy = std::visit(
          [this, &pos1, &pos2, q1q2](auto &p3m_ptr) {
            auto const &p3m = *p3m_ptr;
            auto energy = 0.;
            elc.dielectric_layers_contribution(
                *m_box_geo, pos1, pos2, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  energy += p3m.pair_energy(q_eff, d.norm());
                });
            elc.dielectric_layers_contribution(
                *m_box_geo, pos2, pos1, q1q2,
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
  void add_pair_force_corrections(Utils::Vector3d const &pos1,
                                  Utils::Vector3d const &pos2,
                                  Utils::Vector3d &p1f_asym,
                                  Utils::Vector3d &p2f_asym,
                                  double q1q2) const {
    if (elc.dielectric_contrast_on) {
      std::visit(
          [this, &pos1, &pos2, &p1f_asym, &p2f_asym, q1q2](auto &p3m_ptr) {
            auto const &p3m = *p3m_ptr;
            elc.dielectric_layers_contribution(
                *m_box_geo, pos1, pos2, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  p1f_asym += p3m.pair_force(q_eff, d, d.norm());
                });
            elc.dielectric_layers_contribution(
                *m_box_geo, pos2, pos1, q1q2,
                [&](double q_eff, Utils::Vector3d const &d) {
                  p2f_asym += p3m.pair_force(q_eff, d, d.norm());
                });
          },
          base_solver);
    }
  }

  /** @brief Calculate long-range electrostatic energy with corrections. */
  double long_range_energy() const;
  /** @brief Accumulate long-range electrostatic forces with corrections. */
  void add_long_range_forces() const;

private:
  /** Check if a charged particle is in the gap region. */
  void check_gap(Particle const &p) const;
  double tune_far_cut() const;
  void adapt_solver();
  /** pairwise contributions from the lowest and top layers to the energy */
  double dipole_energy() const;
  void add_dipole_force() const;
  double z_energy() const;
  void add_z_force() const;

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

  /// the force calculation
  void add_force() const;
  /// the energy calculation
  double calc_energy() const;

  void visit_base_solver(auto &&visitor) const {
    std::visit(visitor, base_solver);
  }
};

#endif // ESPRESSO_P3M

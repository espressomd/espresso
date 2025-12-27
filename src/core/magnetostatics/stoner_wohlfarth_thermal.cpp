/*
 * Copyright (C) 2025 The ESPResSo project
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

#include <config/config.hpp>

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH

#include "Particle.hpp"
#include "cell_system/CellStructure.hpp"
#include "cells.hpp"
#include "constraints/Constraints.hpp"
#include "constraints/HomogeneousMagneticField.hpp"
#include "errorhandling.hpp"
#include "random.hpp"
#include "rotation.hpp"
#include "system/System.hpp"
#include "thermostat.hpp"
#include "virtual_sites/relative.hpp"

#include <utils/Vector.hpp>
#include <utils/uniform.hpp>

#include <nlopt.hpp>

#include <cassert>
#include <cmath>
#include <numbers>
#include <tuple>
#include <utility>
#include <vector>

// small perturbation to avoid starting exactly at a stationary point
constexpr static double eps_phi = 1e-3;
// absolute error precision required for the optimiser
constexpr static double eps_abs = 1e-15;
// relative error precision required for the optimiser
constexpr static double eps_rel = 1e-15;

/**
 * @brief Objective (energy) function for the Stoner-Wohlfarth phi minimisation.
 *
 * Evaluate the magnetic energy (normalized by the anisotropy field) for a
 * given in-plane angle phi according to Eq. 5 in @cite mostarac25a.
 * Assumes minima lie in the plane phi = zeta and uses trig identities
 * to reduce the expression.
 *
 * @param n Number of optimization variables (should be 1: phi).
 * @param x Pointer to variables; x[0] is the angle phi.
 * @param grad If non-null, gradient is written to grad[0].
 * @param my_func_data Pointer to a double[2] array with {theta, h}.
 * @return Energy value for the given phi.
 */
static double phi_objective(unsigned n, const double *x, double *grad,
                            void *my_func_data) {
  auto const phi = x[0];
  auto const *params = reinterpret_cast<double const *>(my_func_data);
  auto const theta = params[0];
  auto const h = params[1];
  if (grad) {
    grad[0] = std::sin(2. * (phi - theta)) + 2. * h * std::sin(phi);
  }
  return -0.5 - 0.5 * std::cos(2. * (phi - theta)) - 2. * h * std::cos(phi);
}

/**
 * @brief Find the in-plane angle phi corresponding to the correct
 * energy minimum for the thermal Stoner-Wohlfarth particles.
 *
 * @param theta Angle between anisotropy director and external field (rad).
 * @param h Reduced field (external + dipolar) normalised by H_k.
 * @param phi0 Initial in-plane angle guess (rad).
 * @param ani_param Inverse thermal energy factor (@f$1/(k_B T V)@f$ scaled).
 * @param tau0_inv Attempt frequency inverse (1/tau0).
 * @param dt Time increment for switching probability.
 * @param noise Uniform random number in (0,1) used for the kinetic MC step.
 * @return In-plane angle phi in range @f$ [0,2\pi) @f$.
 */
static double get_phi_at_energy_min(double theta, double h, double phi0,
                                    double ani_param, double tau0_inv,
                                    double dt, double const &noise) {

  auto constexpr pi = std::numbers::pi_v<double>;
  auto constexpr two_pi = 2. * pi;

  // critical field, above which there is only one minimum (no need to do the
  // thermal step); Eq. 6 in @cite mostarac25a.
  auto const h_crit = std::pow(std::pow(std::sin(theta), 2. / 3.) +
                                   std::pow(std::cos(theta), 2. / 3.),
                               -3. / 2.);
  nlopt::opt opt(nlopt::LD_MMA, 1);
  double params[] = {theta, h};

  opt.set_min_objective(phi_objective, &params);
  opt.set_ftol_rel(eps_rel);
  opt.set_ftol_abs(eps_abs);
  std::vector<double> phi(1);

  // make initial guess from previous position plus an arbitrary perturbation
  phi[0] = phi0 + eps_phi;
  auto min1 = 0.;
  opt.optimize(phi, min1);
  auto const phi_min1 = std::fmod(phi[0], two_pi);
  auto solution = phi_min1;
  if (std::fabs(h) < h_crit) {
    opt.set_max_objective(phi_objective, &params);
    phi[0] = phi0 + eps_phi;
    auto max1 = 0.;
    opt.optimize(phi, max1);
    auto const phi_max1 = std::fmod(phi[0], two_pi);
    phi[0] = std::fmod(phi_max1 + pi, two_pi);
    auto max2 = 0.;
    opt.optimize(phi, max2);
    // Eqs. 12 in @cite mostarac25a.
    auto const b1 = std::abs(max1 - min1) * ani_param;
    auto const b2 = std::abs(max2 - min1) * ani_param;
    auto const b_min = (b1 < b2) ? b1 : b2;
    // Eq. 13 in @cite mostarac25a.
    auto const tau_inv = tau0_inv * exp(-b_min);
    // switching probability (without backflip)
    auto const p12 = 1. - exp(-dt * tau_inv);
    // if MC move accepted, find the location of the other minimum
    if (noise < p12) {
      opt.set_min_objective(phi_objective, &params);
      phi[0] = std::fmod(phi_min1 + pi + eps_phi, two_pi);
      /*try to find another minimimum from the other side*/
      auto min2 = 0.;
      opt.optimize(phi, min2);
      auto const phi_min2 = std::fmod(phi[0], two_pi);
      solution = phi_min2;
    }
  }
  return std::fmod(solution + two_pi, two_pi);
}

/**
 * @brief Collect external homogeneous magnetic field from active constraints.
 *
 * Iterate over constraints and sum the homogeneous magnetic field vectors
 * provided by @ref Constraints::HomogeneousMagneticField objects.
 *
 * @return The total external homogeneous magnetic field.
 */
static auto get_external_field(Constraints::Constraints const &constraints) {
  using HomogeneousMagneticField = ::Constraints::HomogeneousMagneticField;
  Utils::Vector3d ext_fld = {0., 0., 0.};
  for (auto const &constraint : constraints) {
    auto ptr = std::dynamic_pointer_cast<HomogeneousMagneticField>(constraint);
    if (ptr) {
      ext_fld += ptr->H();
    }
  }
  return ext_fld;
}

/**
 * @brief Simplified Stoner-Wohlfarth update in field-free case.
 *
 * @param[in,out] p Virtual particle to update.
 * @param e_k Anisotropy director of the reference particle.
 * @param kT Thermal energy from thermostat.
 * @param noise Uniform random number in (0,1) used for the kinetic MC step.
 */
void stoner_wohlfarth_no_field(Particle &p, Utils::Vector3d const &e_k,
                               double const kT, double const noise) {
  auto constexpr pi = std::numbers::pi_v<double>;
  auto const ani_param = p.magnetic_anisotropy_energy() / kT;
  auto const tau_inv = p.stoner_wohlfarth_tau0_inv() * std::exp(-ani_param);
  auto const p12 = 1. - std::exp(-p.stoner_wohlfarth_dt_incr() * tau_inv);
  auto const kernel = [&](bool flip) {
    auto const sat_mag = (flip ? -1. : +1.) * p.saturation_magnetization();
    auto const [quat, dipm] = convert_dip_to_quat(sat_mag * e_k);
    p.stoner_wohlfarth_phi_0() = flip ? pi : 0.;
    p.dipm() = dipm;
    p.quat() = quat;
  };
  auto const flip = noise < p12;
  if (p.stoner_wohlfarth_phi_0() == 0.) {
    kernel(flip);
  } else if (p.stoner_wohlfarth_phi_0() == pi) {
    kernel(not flip);
  } else {
    auto const dist_0 = std::abs(p.stoner_wohlfarth_phi_0() - 0.);
    auto const dist_pi = std::abs(p.stoner_wohlfarth_phi_0() - pi);
    // Compare the differences and determine the nearest angle (simplifies
    // down to an XNOR operation, i.e. equality operator for boolean types)
    kernel(flip == (dist_0 < dist_pi));
  }
}

/**
 * @brief Update virtual site dipole moment according to the full in-field
 * (incl. dipole field) thermal Stoner-Wohlfarth model (incl. the kinetic MC
 * step)
 *
 * @param[in,out] p Virtual particle to update (modified).
 * @param e_k Anisotropy director of the reference particle.
 * @param ext_fld_dpl External homogeneous magnetic field + total dipolar field
 * acting on the particle.
 * @param kT Thermal energy from thermostat.
 * @param noise Uniform random number in (0,1) used for the kinetic MC step.
 */
static void stoner_wohlfarth_main(Particle &p, Utils::Vector3d const &e_k,
                                  Utils::Vector3d const &ext_fld_dpl,
                                  double const kT, double const noise) {

  auto constexpr pi = std::numbers::pi_v<double>;
  auto constexpr pi_half = pi / 2.;

  // reduced field; Eq. 4 in @cite mostarac25a.
  auto h = ext_fld_dpl.norm() * p.magnetic_anisotropy_field_inv();
  auto e_h = ext_fld_dpl.normalized();
  auto theta = std::acos(e_h * e_k);
  if (theta > pi_half) {
    theta = pi - theta;
    h = -h;
    e_h = -e_h;
  }
  auto const rot_axis =
      vector_product(vector_product(e_h, e_k), e_h).normalized();
  auto const ani_param = p.magnetic_anisotropy_energy() / kT;
  auto const phi = get_phi_at_energy_min(
      theta, h, p.stoner_wohlfarth_phi_0(), ani_param,
      p.stoner_wohlfarth_tau0_inv(), p.stoner_wohlfarth_dt_incr(), noise);
  auto const mom = e_h * std::cos(phi) + rot_axis * std::sin(phi);
  p.stoner_wohlfarth_phi_0() = phi;
  auto const [quat, dipm] =
      convert_dip_to_quat(mom * p.saturation_magnetization());
  p.dipm() = dipm;
  p.quat() = quat;
}

/**
 * @brief Run magnetodynamics update for local virtual particles.
 *
 * Iterate over local particles and update the dipole moment of virtual
 * particles according to the thermal Stoner-Wohlfarth model.
 * Collect active homogeneous external magnetic fields from constraints and
 * add the per-particle dipolar contribution before performing either the
 * simplified no-field update or the full thermal Stoner-Wohlfarth update.
 */
void System::System::integrate_magnetodynamics() {
  // collect HomogeneousMagneticFields if active
  auto const ext_fld = get_external_field(*constraints);
  auto const kT = thermostat->kT;
  cell_structure->for_each_local_particle([&](Particle &p) {
    if (not p.is_virtual() or not p.stoner_wohlfarth_is_enabled()) {
      return;
    }
    auto *p_ref = get_reference_particle(*cell_structure, p);
    if (not p_ref) {
      return;
    }
    assert(thermostat->thermo_switch & THERMO_LANGEVIN);
    auto const &langevin = *thermostat->langevin;
    auto const e_k = p_ref->calc_director();
    auto const ext_fld_dpl = ext_fld + p.dip_fld();
    auto const random_ints =
        Random::philox_4_uint64s<RNGSalt::THERMAL_STONER_WOHLFARTH>(
            langevin.rng_counter(), langevin.rng_seed(), p.id());
    auto const noise = Utils::uniform(random_ints[0]);
    if (ext_fld_dpl.norm2() == 0.) {
      stoner_wohlfarth_no_field(p, e_k, kT, noise);
    } else {
      // full Stoner-Wohlfarth update with external + dipolar field
      stoner_wohlfarth_main(p, e_k, ext_fld_dpl, kT, noise);
    }
  });
}

#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH

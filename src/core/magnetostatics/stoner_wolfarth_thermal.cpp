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

#include "config/config.hpp"

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
#define TWO_M_PI 2 * M_PI

#include "cells.hpp"
#include "constraints/Constraints.hpp"
#include "constraints/HomogeneousMagneticField.hpp"
#include "errorhandling.hpp"

#include <nlopt.hpp>

#include "magnetostatics/stoner_wolfarth_thermal.hpp"
#include "random.hpp"
#include "rotation.hpp"
#include "thermostat.hpp"
#include <cmath>
#include <tuple>
#include <utility>
#include <vector>

namespace {
// small perturbation to avoid starting exactly at a stationary point
constexpr double eps_phi = 1e-3;
// absolute error precision required for the optimiser
constexpr double eps_abs = 1e-15;
// relative error precision required for the optimiser
constexpr double eps_rel = 1e-15;

/**
 * @brief Get real particle tracked by a virtual site.
 *
 * @param cell_structure Cell structure.
 * @param p Virtual site.
 * @return Pointer to real particle.
 */
static Particle *get_reference_particle(CellStructure &cell_structure,
                                        Particle const &p) {
  auto const &vs_rel = p.vs_relative();
  if (vs_rel.to_particle_id == -1) {
    runtimeErrorMsg() << "Particle with id " << p.id()
                      << " is a dangling virtual site";
    return nullptr;
  }
  auto p_ref_ptr = cell_structure.get_local_particle(vs_rel.to_particle_id);
  if (!p_ref_ptr) {
    runtimeErrorMsg() << "No real particle with id " << vs_rel.to_particle_id
                      << " for virtual site with id " << p.id();
  }
  return p_ref_ptr;
}

/**
 * @brief Objective (energy) function for the Stoner–Wohlfarth phi minimisation.
 *
 * Evaluates the magnetic energy (normalized by the anisotropy field) for a
 * given in-plane angle phi according to Eq. 5 in
 * https://doi.org/10.1103/PhysRevB.111.014438. Assumes minima lie in the
 * plane phi = zeta and uses trig identities to reduce the expression.
 *
 * @param n Number of optimization variables (should be 1: phi).
 * @param x Pointer to variables; x[0] is the angle phi.
 * @param grad If non-null, gradient is written to grad[0].
 * @param my_func_data Pointer to a double[2] array with {theta, h}.
 * @return Energy value for the given phi.
 */
double phi_objective(unsigned n, const double *x, double *grad,
                     void *my_func_data) {
  double const phi = x[0];
  double const *params = (double const *)my_func_data;
  double const theta = params[0];
  double const h = params[1];
  if (grad) {
    grad[0] = std::sin(2 * (phi - theta)) + 2 * h * std::sin(phi);
  }
  return -0.5 - 0.5 * std::cos(2 * (phi - theta)) - 2 * h * std::cos(phi);
}

/**
 * @brief Find the in-plane angle phi corresponding to the correct
 *        energy minimum for the thermal Stoner–Wohlfarth particles.
 *
 * @param theta Angle between anisotropy director and external field (rad).
 * @param h Reduced field (external + dipolar) normalised by H_k.
 * @param phi0 Initial in‑plane angle guess (rad).
 * @param ani_param Inverse thermal energy factor (1/(k_B T V) scaled).
 * @param tau0_inv Attempt frequency inverse (1/tau0).
 * @param dt Time increment for switching probability.
 * @param noise Uniform random number in (0,1) used for the kinetic Monte‑Carlo
 * step.
 * @return In‑plane angle phi in range [0,2π).
 */
double get_phi_at_energy_min(double theta, double h, double phi0,
                             double ani_param, double tau0_inv, double dt,
                             const double &noise) {

  // critical filed, above which there is only one minimum  (no need to do the
  // thermal step); Eq. 6 in https://doi.org/10.1103/PhysRevB.111.014438.
  double const h_crit = std::pow(std::pow(std::sin(theta), 2.0 / 3) +
                                     std::pow(std::cos(theta), 2.0 / 3),
                                 -3.0 / 2);
  nlopt::opt opt(nlopt::LD_MMA, 1);
  double params[] = {theta, h};

  opt.set_min_objective(phi_objective, &params);
  opt.set_ftol_rel(
      eps_rel); // Set the relative tolerance for the objective function value
  opt.set_ftol_abs(
      eps_abs); // Set the relative tolerance for the objective function value
  std::vector<double> phi(1);

  phi[0] = phi0 + eps_phi; /* make initial guess from previos position plus
                                   an arbitrary perturbation*/
  double min1; /* this is the actuall value of the energy from minimiser */
  opt.optimize(phi, min1);
  double const phi_min1 = fmod(phi[0], TWO_M_PI);
  double sol = phi_min1;
  if (fabs(h) < h_crit) {
    opt.set_max_objective(phi_objective, &params);
    phi[0] = phi0 + eps_phi;
    double max1;
    opt.optimize(phi, max1);
    double const phi_max1 = fmod(phi[0], TWO_M_PI);
    phi[0] = fmod(phi_max1 + M_PI, 2 * M_PI);
    double max2;
    opt.optimize(phi, max2);
    // Eqs. 12 in https://doi.org/10.1103/PhysRevB.111.014438.
    double const b1 = std::abs(max1 - min1) * ani_param;
    double const b2 = std::abs(max2 - min1) * ani_param;
    double const b_min = (b1 < b2) ? b1 : b2;
    // Eq. 13 in https://doi.org/10.1103/PhysRevB.111.014438.
    double const tau_inv = tau0_inv * exp(-b_min);
    // switching probability (without backflip)
    double const p12 = 1. - exp(-dt * tau_inv);
    // if MC move accepted, find the location of the other minimum
    if (noise < p12) {
      opt.set_min_objective(phi_objective, &params);
      phi[0] = fmod(phi_min1 + M_PI + eps_phi, 2 * M_PI);
      /*try to find another minimimum from the other side*/
      double min2;
      opt.optimize(phi, min2);
      double const phi_min2 = fmod(phi[0], TWO_M_PI);
      sol = phi_min2;
    }
  }
  return fmod(sol + TWO_M_PI, TWO_M_PI);
}

/**
 * @brief Collect external homogeneous magnetic field from active constraints.
 *
 * Iterates over System::get_system().constraints and sums the homogeneous
 * magnetic field vectors provided by Constraints::HomogeneousMagneticField
 * constraint objects.
 *
 * @return Utils::Vector3d The total external homogeneous magnetic field.
 */
const Utils::Vector3d get_external_field() {
  Utils::Vector3d ext_fld = {0., 0., 0.};
  auto &system = System::get_system();
  for (auto const &constraint : *system.constraints) {
    auto ptr =
        std::dynamic_pointer_cast<::Constraints::HomogeneousMagneticField>(
            constraint);
    if (ptr) {
      ext_fld += ptr->H();
    }
  }
  return ext_fld;
}
} // namespace
/**
 * @brief Simplified Stoner–Wohlfarth update in field free case.
 *
 * @param p Virtual particle to update (modified).
 * @param pi Reference particle providing the anisotropy director (read-only).
 * @param kT Thermal energy from thermostat.
 * @param noise Uniform random number in (0,1) used for the kinetic Monte‑Carlo
 * step.
 */
void stoner_wohlfarth_no_field(Particle &p, Particle &pi, double const kT,
                               const double &noise) {

  auto const e_k = pi.calc_director();
  double const ani_param = p.magnetic_anisotropy_energy() / kT;
  double const tau_inv = p.stoner_wolfarth_tau0_inv() * exp(-ani_param);
  double const p12 = 1. - exp(-p.stoner_wolfarth_dt_incr() * tau_inv);
  if (noise < p12) {
    if (p.stoner_wolfarth_phi_0() == 0) {
      auto const [quat, dipm] =
          convert_dip_to_quat(p.saturation_magnetization() * -e_k);
      p.stoner_wolfarth_phi_0() = M_PI;
      p.dipm() = dipm;
      p.quat() = quat;
    } else if (p.stoner_wolfarth_phi_0() == M_PI) {
      auto const [quat, dipm] =
          convert_dip_to_quat(p.saturation_magnetization() * e_k);
      p.stoner_wolfarth_phi_0() = 0;
      p.dipm() = dipm;
      p.quat() = quat;
    } else {
      double const diff_0 = std::abs(p.stoner_wolfarth_phi_0() - 0);
      double const diff_PI = std::abs(p.stoner_wolfarth_phi_0() - M_PI);
      // Compare the differences and determine the closer angle
      if (diff_0 < diff_PI) {
        auto const [quat, dipm] =
            convert_dip_to_quat(p.saturation_magnetization() * -e_k);
        p.stoner_wolfarth_phi_0() = M_PI;
        p.dipm() = dipm;
        p.quat() = quat;
      } else {
        auto const [quat, dipm] =
            convert_dip_to_quat(p.saturation_magnetization() * e_k);
        p.stoner_wolfarth_phi_0() = 0;
        p.dipm() = dipm;
        p.quat() = quat;
      }
    }
  } else {
    if (p.stoner_wolfarth_phi_0() == 0) {
      auto const [quat, dipm] =
          convert_dip_to_quat(p.saturation_magnetization() * e_k);
      p.dipm() = dipm;
      p.quat() = quat;
    } else if (p.stoner_wolfarth_phi_0() == M_PI) {
      auto const [quat, dipm] =
          convert_dip_to_quat(p.saturation_magnetization() * -e_k);
      p.dipm() = dipm;
      p.quat() = quat;
    } else {
      double const diff_0 = std::abs(p.stoner_wolfarth_phi_0() - 0);
      double const diff_PI = std::abs(p.stoner_wolfarth_phi_0() - M_PI);
      // Compare the differences and determine the closer angle
      if (diff_0 < diff_PI) {
        auto const [quat, dipm] =
            convert_dip_to_quat(p.saturation_magnetization() * e_k);
        p.stoner_wolfarth_phi_0() = 0;
        p.dipm() = dipm;
        p.quat() = quat;
      } else {
        auto const [quat, dipm] =
            convert_dip_to_quat(p.saturation_magnetization() * -e_k);
        p.stoner_wolfarth_phi_0() = M_PI;
        p.dipm() = dipm;
        p.quat() = quat;
      }
    }
  }
}
/**
 * @brief Update virtual site dipole moment accodring to the full in-field
 * (incl. dipole field) thermal Stoner-Wohlfarth model (incl. the kinetic MC
 * step)
 *
 * @param p Virtual particle to update (modified).
 * @param pi Reference particle providing the anisotropy director (read-only).
 * @param ext_fld_dpl External homogeneous magnetic field + total dipolar field
 * acting on the particle.
 * @param kT Thermal energy from thermostat.
 * @param noise Uniform random number in (0,1) used for the kinetic Monte‑Carlo
 * step.
 */
void stoner_wohlfarth_main(Particle &p, Particle &pi,
                           const Utils::Vector3d &ext_fld_dpl, double const kT,
                           const double &noise) {
  // reduced field; Eq. 4 in https://doi.org/10.1103/PhysRevB.111.014438.
  double h = ext_fld_dpl.norm() * p.magnetic_anisotropy_field_inv();
  auto e_h = ext_fld_dpl.normalized();
  // calc_director() result already normalised
  auto const e_k = pi.calc_director();
  double theta = std::acos(e_h * e_k);
  if (theta > M_PI_2) {
    theta = M_PI - theta;
    h = -h;
    e_h = -e_h;
  }
  auto const rot_axis =
      vector_product(vector_product(e_h, e_k), e_h).normalized();
  double const ani_param = p.magnetic_anisotropy_energy() / kT;
  auto const phi = get_phi_at_energy_min(
      theta, h, p.stoner_wolfarth_phi_0(), ani_param,
      p.stoner_wolfarth_tau0_inv(), p.stoner_wolfarth_dt_incr(), noise);
  auto const mom = e_h * std::cos(phi) + rot_axis * std::sin(phi);
  p.stoner_wolfarth_phi_0() = phi;
  auto const [quat, dipm] =
      convert_dip_to_quat(mom * p.saturation_magnetization());
  p.dipm() = dipm;
  p.quat() = quat;
}
/**
 * @brief Run magnetodynamics update for local virtual particles.
 *
 * Iterates over local particles and updates the dipole moment of virtual
 * particles according to the thermal Stoner–Wohlfarth model. Collects
 * active homogeneous external magnetic fields from constraints and adds the
 * per-particle dipolar contribution before performing either the simplified
 * no-field update or the full thermal Stoner–Wohlfarth update.
 *
 * @param cell_structure CellStructure providing access to local particles.
 * @param thermostat Const reference to Thermostat used to access Philox RNG
 * state and seeds.
 */
void run_magnetodynamics(CellStructure &cell_structure,
                         Thermostat::Thermostat const &thermostat) {
  /* collect HomogeneousMagneticFields if active */
  auto const ext_fld = get_external_field();
  auto const kT = thermostat.kT;
  cell_structure.for_each_local_particle([&](Particle &p) {
    /* collect particle data */
    if (!p.is_virtual() || !p.stoner_wolfarth_is_enabled()) {
      return;
    }

    auto *pref = get_reference_particle(cell_structure, p);
    if (!pref) {
      return;
    }
    auto &pi = *pref;
    auto const ext_fld_dpl = ext_fld + p.dip_fld();
    // if no external field and no dipolar field, do simplified Stoner-Wohlfarth
    // update
    auto const random_int =
        Random::philox_4_uint64s<RNGSalt::THERMAL_STONER_WOHLFARTH>(
            thermostat.get_philox_counter(), thermostat.get_philox_seed(),
            p.id());
    double const random_uniform_dist_cast = Utils::uniform(
        static_cast<std::size_t>(random_int[0])); // uniform (0,1)
    if (ext_fld_dpl.norm() == 0.) {
      stoner_wohlfarth_no_field(p, pi, kT, random_uniform_dist_cast);
      return;
    }
    // full Stoner-Wohlfarth update with external + dipolar field
    stoner_wohlfarth_main(p, pi, ext_fld_dpl, kT, random_uniform_dist_cast);
  });
}
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH

/*
 * Copyright (C) 2010-2026 The ESPResSo project
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

#pragma once

/**
 * @file
 * @brief Cold per-particle parameter PODs stored in @ref ParticleStore host
 * sidecars.
 *
 * These structs live in this shared header so both @ref ParticleStore (which
 * owns the sidecar @c std::vector storage) and @ref Particle can name the
 * complete type. Their @c serialize members define the checkpoint / inter-rank
 * wire format.
 */

#include <config/config.hpp>

#include <utils/quaternion.hpp>

#ifdef ESPRESSO_ENGINE
/** Properties of a self-propelled particle. */
struct ParticleParametersSwimming {
  /** Imposed constant force. */
  double f_swim = 0.;
  /** Is the particle a swimmer. */
  bool swimming = false;
  /** Whether f_swim is applied to the particle or to the fluid. */
  bool is_engine_force_on_fluid = false;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & f_swim & swimming & is_engine_force_on_fluid;
  }
};
#endif

#ifdef ESPRESSO_THERMAL_STONER_WOHLFARTH
/** Properties for thermal Stoner-Wohlfarth magnetodynamics. */
struct ThermalStonerWohlfarthParameters {
  /**
   * Flag to distinguish virtual particles carrying the dipole moment in
   * the thermal Stoner-Wohlfarth model from other types of virtual sites.
   */
  bool is_enabled = false;
  /** angle between the director and dipole moment of a Stoner-Wohlfarth
   * particle */
  double phi0 = 0.;
  /** saturation magnetisation of a polarizable particle */
  double sat_mag = 0.;
  /**
   * @brief Inverse anisotropy field in reduced units.
   * anisotropy field = 2.*K1/(mu0 * Ms) in [A / m] where K1 is the magnetic
   * anisotropy constant [kg / (m s^2)].
   */
  double ani_fld_inv = 0.;
  /**
   * @brief Magnetic anisotropy energy (K1 * V) in units of energy.
   * Related to ani_param from Eq.3 in @cite mostarac25a
   * by: ani_param = ani_energy / kT
   */
  double ani_energy = 0.;
  /** Browns attempt frequency. Prefactor from Eq.9 in @cite mostarac25a.  */
  double tau0_inv = 0.;
  /** time units parameter for the kinetic Monte Carlo step */
  double dt_incr = 0.;

  template <class Archive> void serialize(Archive &ar, long int /* version */) {
    ar & is_enabled & phi0 & sat_mag & ani_fld_inv & ani_energy & tau0_inv &
        dt_incr;
  }
};
#endif // ESPRESSO_THERMAL_STONER_WOHLFARTH

#ifdef ESPRESSO_VIRTUAL_SITES_RELATIVE
/** The following properties define, with respect to which real particle a
 *  virtual site is placed and at what distance. The relative orientation of
 *  the vector pointing from real particle to virtual site with respect to the
 *  orientation of the real particle is stored in the virtual site's
 *  quaternion attribute.
 */
struct VirtualSitesRelativeParameters {
  int to_particle_id = -1;
  double distance = 0.;
  /** Relative position of the virtual site. */
  Utils::Quaternion<double> rel_orientation =
      Utils::Quaternion<double>::identity();
  /** Orientation of the virtual particle in the body fixed frame. */
  Utils::Quaternion<double> quat = Utils::Quaternion<double>::identity();

  template <class Archive> void serialize(Archive &ar, long int) {
    ar & to_particle_id;
    ar & distance;
    ar & rel_orientation;
    ar & quat;
  }
};
#endif // ESPRESSO_VIRTUAL_SITES_RELATIVE

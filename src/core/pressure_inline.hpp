/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
/** \file
 *  Pressure calculation. Really similar to energy.hpp.
 */

#ifndef CORE_PRESSURE_INLINE_HPP
#define CORE_PRESSURE_INLINE_HPP

#include "Observable_stat.hpp"
#include "exclusions.hpp"
#include "forces_inline.hpp"
#include "integrate.hpp"
#include "npt.hpp"
#include "pressure.hpp"

#include <utils/math/tensor_product.hpp>

extern Observable_stat_wrapper obs_scalar_pressure;
extern Observable_stat_wrapper obs_pressure_tensor;
extern Observable_stat_non_bonded_wrapper obs_scalar_pressure_non_bonded;
extern Observable_stat_non_bonded_wrapper obs_pressure_tensor_non_bonded;

/** Calculate non bonded energies between a pair of particles.
 *  @param p1        pointer to particle 1.
 *  @param p2        pointer to particle 2.
 *  @param d         vector between p1 and p2.
 *  @param dist      distance between p1 and p2.
 */
inline void add_non_bonded_pair_virials(Particle const &p1, Particle const &p2,
                                        Utils::Vector3d const &d, double dist) {
#ifdef EXCLUSIONS
  if (do_nonbonded(p1, p2))
#endif
  {
    auto const force = calc_non_bonded_pair_force(p1, p2, d, dist);
    obs_scalar_pressure.local.non_bonded_contribution(
        p1.p.type, p2.p.type)[0] += d * force;

    /* pressure tensor part */
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        obs_pressure_tensor.local.non_bonded_contribution(
            p1.p.type, p2.p.type)[k * 3 + l] += force[k] * d[l];

    auto const p1molid = p1.p.mol_id;
    auto const p2molid = p2.p.mol_id;
    if (p1molid == p2molid) {
      obs_scalar_pressure_non_bonded.local.non_bonded_intra_contribution(
          p1.p.type, p2.p.type)[0] += d * force;

      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          obs_pressure_tensor_non_bonded.local.non_bonded_intra_contribution(
              p1.p.type, p2.p.type)[k * 3 + l] += force[k] * d[l];
    } else {
      obs_scalar_pressure_non_bonded.local.non_bonded_inter_contribution(
          p1.p.type, p2.p.type)[0] += d * force;

      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          obs_pressure_tensor_non_bonded.local.non_bonded_inter_contribution(
              p1.p.type, p2.p.type)[k * 3 + l] += force[k] * d[l];
    }
  }

#ifdef ELECTROSTATICS
  if (!obs_scalar_pressure.local.coulomb.empty()) {
    /* real space Coulomb */
    auto const p_coulomb = Coulomb::pair_pressure(p1, p2, d, dist);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        obs_pressure_tensor.local.coulomb[i * 3 + j] += p_coulomb[i][j];
      }
    }

    obs_scalar_pressure.local.coulomb[0] += trace(p_coulomb);
  }
#endif /*ifdef ELECTROSTATICS */

#ifdef DIPOLES
  /* real space magnetic dipole-dipole */
  if (dipole.method != DIPOLAR_NONE) {
    fprintf(stderr, "calculating pressure for magnetostatics which doesn't "
                    "have it implemented\n");
  }
#endif /*ifdef DIPOLES */
}

boost::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_virial_pressure_tensor(Bonded_ia_parameters const &iaparams,
                                   Particle const &p1, Particle const &p2) {
  auto const dx = get_mi_vector(p1.r.p, p2.r.p, box_geo);
  Utils::Vector3d torque{};
  auto const result = calc_bond_pair_force(p1, p2, iaparams, dx, torque);
  if (result) {
    auto const &force = result.get();

    return tensor_product(force, dx);
  }

  return {};
}

boost::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_three_body_pressure_tensor(Bonded_ia_parameters const &iaparams,
                                       Particle const &p1, Particle const &p2,
                                       Particle const &p3) {
  switch (iaparams.type) {
  case BONDED_IA_ANGLE_HARMONIC:
  case BONDED_IA_ANGLE_COSINE:
  case BONDED_IA_ANGLE_COSSQUARE:
  case BONDED_IA_TABULATED_ANGLE: {
    auto const dx21 = -get_mi_vector(p1.r.p, p2.r.p, box_geo);
    auto const dx31 = get_mi_vector(p3.r.p, p1.r.p, box_geo);

    auto const result = calc_bonded_three_body_force(iaparams, p1, p2, p3);
    if (result) {
      Utils::Vector3d force2, force3;
      std::tie(std::ignore, force2, force3) = result.get();

      return tensor_product(force2, dx21) + tensor_product(force3, dx31);
    }
  }
  default:
    runtimeWarningMsg() << "Unsupported bond type " +
                               std::to_string(iaparams.type) +
                               " in pressure calculation.";
    return Utils::Matrix<double, 3, 3>{};
  }
}

inline boost::optional<Utils::Matrix<double, 3, 3>>
calc_bonded_pressure_tensor(Bonded_ia_parameters const &iaparams, Particle &p1,
                            Utils::Span<Particle *> partners) {
  switch (iaparams.num) {
  case 1:
    return calc_bonded_virial_pressure_tensor(iaparams, p1, *partners[0]);
  case 2:
    return calc_bonded_three_body_pressure_tensor(iaparams, p1, *partners[0],
                                                  *partners[1]);
  default:
    runtimeWarningMsg() << "Unsupported bond type " +
                               std::to_string(iaparams.type) +
                               " in pressure calculation.";
    return Utils::Matrix<double, 3, 3>{};
  }
}

inline bool add_bonded_pressure_tensor(Particle &p1, int bond_id,
                                       Utils::Span<Particle *> partners) {
  auto const &iaparams = bonded_ia_params[bond_id];

  auto const result = calc_bonded_pressure_tensor(iaparams, p1, partners);
  if (result) {
    auto const &tensor = result.get();

    obs_scalar_pressure.local.bonded_contribution(bond_id)[0] += trace(tensor);

    /* pressure tensor part */
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        obs_pressure_tensor.local.bonded_contribution(bond_id)[k * 3 + l] +=
            tensor[k][l];

    return false;
  }

  return true;
}

/** Calculate kinetic pressure (aka energy) for one particle.
 *  @param p1 particle for which to calculate pressure
 *  @param v_comp flag which enables compensation of the velocities required
 *     for deriving a pressure reflecting \ref nptiso_struct::p_inst
 *     (hence it only works with domain decomposition); naturally it
 *     therefore doesn't make sense to use it without NpT.
 */
inline void add_kinetic_virials(Particle const &p1, bool v_comp) {
  if (p1.p.is_virtual)
    return;

  /* kinetic energy */
  if (v_comp) {
    // velocity at half the time step: v(t + dt / 2) =
    // v(t + dt) - a(t) * dt / 2 = v(t + dt) - F(t) * dt / m / 2
    auto const v = p1.m.v - p1.f.f * (time_step / (2. * p1.p.mass));
    obs_scalar_pressure.local.kinetic[0] += v.norm2() * p1.p.mass;
  } else {
    obs_scalar_pressure.local.kinetic[0] += p1.m.v.norm2() * p1.p.mass;
  }

  for (int k = 0; k < 3; k++)
    for (int l = 0; l < 3; l++)
      obs_pressure_tensor.local.kinetic[k * 3 + l] +=
          (p1.m.v[k]) * (p1.m.v[l]) * p1.p.mass;
}

#endif

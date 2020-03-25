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

#include "exclusions.hpp"
#include "forces_inline.hpp"
#include "integrate.hpp"
#include "npt.hpp"
#include "pressure.hpp"

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
    *obsstat_nonbonded(&virials, p1.p.type, p2.p.type) += d * force;

    /* stress tensor part */
    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        obsstat_nonbonded(&p_tensor, p1.p.type, p2.p.type)[k * 3 + l] +=
            force[k] * d[l];

    auto const p1molid = p1.p.mol_id;
    auto const p2molid = p2.p.mol_id;
    if (p1molid == p2molid) {
      *obsstat_nonbonded_intra(&virials_non_bonded, p1.p.type, p2.p.type) +=
          d * force;

      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          obsstat_nonbonded_intra(&p_tensor_non_bonded, p1.p.type,
                                  p2.p.type)[k * 3 + l] += force[k] * d[l];
    } else {
      *obsstat_nonbonded_inter(&virials_non_bonded, p1.p.type, p2.p.type) +=
          d * force;

      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          obsstat_nonbonded_inter(&p_tensor_non_bonded, p1.p.type,
                                  p2.p.type)[k * 3 + l] += force[k] * d[l];
    }
  }

#ifdef ELECTROSTATICS
  {
    /* real space Coulomb */
    auto const p_coulomb = Coulomb::pair_pressure(p1, p2, d, dist);

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        p_tensor.coulomb[i * 3 + j] += p_coulomb[i][j];
      }
    }

    virials.coulomb[0] += trace(p_coulomb);
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

/** Calculate bonded virials for one particle.
 *  @param p1 particle for which to calculate virials
 */
inline void add_bonded_virials(Particle const &p1) {

  int i = 0;
  while (i < p1.bl.n) {
    auto const type_num = p1.bl.e[i++];
    Bonded_ia_parameters const &iaparams = bonded_ia_params[type_num];
    if (iaparams.num != 1) {
      i += iaparams.num;
      continue;
    }

    /* fetch particle 2 */
    Particle const *const p2 = cell_structure.get_local_particle(p1.bl.e[i++]);
    if (!p2) {
      // for harmonic spring:
      // if cutoff was defined and p2 is not there it is anyway outside the
      // cutoff, see calc_maximal_cutoff()
      if ((type_num == BONDED_IA_HARMONIC) && (iaparams.p.harmonic.r_cut > 0))
        return;
      runtimeErrorMsg() << "bond broken between particles " << p1.p.identity
                        << " and " << p1.bl.e[i - 1]
                        << " (particles not stored on the same node)";
      return;
    }

    auto const dx = get_mi_vector(p1.r.p, p2->r.p, box_geo);
    Utils::Vector3d torque{};
    auto result = calc_bond_pair_force(p1, *p2, iaparams, dx, torque);
    if (result) {
      auto const &force = result.get();

      *obsstat_bonded(&virials, type_num) += dx * force;

      /* stress tensor part */
      for (int k = 0; k < 3; k++)
        for (int l = 0; l < 3; l++)
          obsstat_bonded(&p_tensor, type_num)[k * 3 + l] += force[k] * dx[l];
    }
  }
}

/** Calculate the contribution to the stress tensor from angular potentials.
 *  The central particle of the three-particle interaction is responsible
 *  for the contribution of the entire interaction - this is the coding
 *  not the physics.
 */
inline void add_three_body_bonded_stress(Particle const &p1) {
  int i = 0;
  while (i < p1.bl.n) {
    /* scan bond list for angular interactions */
    auto const type_num = p1.bl.e[i];
    Bonded_ia_parameters const &iaparams = bonded_ia_params[type_num];

    // Skip non-three-particle-bonds
    if (iaparams.num != 2) // number of partners
    {
      i += 1 + iaparams.num;
      continue;
    }

    switch (iaparams.type) {
    case BONDED_IA_ANGLE_HARMONIC:
    case BONDED_IA_ANGLE_COSINE:
    case BONDED_IA_ANGLE_COSSQUARE:
    case BONDED_IA_TABULATED_ANGLE: {
      Particle const &p2 = *cell_structure.get_local_particle(p1.bl.e[i + 1]);
      Particle const &p3 = *cell_structure.get_local_particle(p1.bl.e[i + 2]);

      auto const dx21 = -get_mi_vector(p1.r.p, p2.r.p, box_geo);
      auto const dx31 = get_mi_vector(p3.r.p, p1.r.p, box_geo);

      auto const result = calc_bonded_three_body_force(iaparams, p1, p2, p3);
      /* three-body bonded interactions contribute to the stress but not the
       * scalar pressure */
      if (result) {
        Utils::Vector3d force2, force3;
        std::tie(std::ignore, force2, force3) = result.get();

        for (int k = 0; k < 3; k++) {
          for (int l = 0; l < 3; l++) {
            obsstat_bonded(&p_tensor, type_num)[3 * k + l] +=
                force2[k] * dx21[l] + force3[k] * dx31[l];
          }
        }
      }
    }
    default:
      runtimeWarningMsg() << "Unsupported bond type " +
                                 std::to_string(iaparams.type) +
                                 " in pressure calculation.";
    }
    i += 3; // bond type and 2 partners
  }
}

/** Calculate kinetic pressure (aka energy) for one particle.
 *  @param p1 particle for which to calculate pressure
 *  @param v_comp flag which enables compensation of the velocities required
 *                for deriving a pressure reflecting \ref nptiso_struct::p_inst
 *                (hence it only works with domain decomposition); naturally it
 *                therefore doesn't make sense to use it without NpT.
 */
inline void add_kinetic_virials(Particle const &p1, int v_comp) {
  if (not p1.p.is_virtual) {
    /* kinetic energy */
    if (v_comp) {
      virials.data.e[0] +=
          ((p1.m.v * time_step) -
           (p1.f.f * (0.5 * Utils::sqr(time_step) / p1.p.mass)))
              .norm2() *
          p1.p.mass;
    } else {
      virials.data.e[0] += Utils::sqr(time_step) * p1.m.v.norm2() * p1.p.mass;
    }

    for (int k = 0; k < 3; k++)
      for (int l = 0; l < 3; l++)
        p_tensor.data.e[k * 3 + l] +=
            (p1.m.v[k] * time_step) * (p1.m.v[l] * time_step) * p1.p.mass;
  }
}

#endif

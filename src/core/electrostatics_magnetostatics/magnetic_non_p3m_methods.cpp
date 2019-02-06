/*
  Copyright (C) 2010-2018 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
    Max-Planck-Institute for Polymer Research, Theory Group

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
/** \file
 * All 3d non P3M methods to deal with the
 * magnetic dipoles
 *
 *  DAWAANR => DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA
 *   Handling of a system of dipoles where no replicas
 *   Assumes minimum image convention for those axis in which the
 *   system is periodic as defined by setmd.
 *
 *   MDDS => Calculates dipole-dipole interaction of a periodic system
 *   by explicitly summing the dipole-dipole interaction over several copies of
 * the system
 *   Uses spherical summation order
 *
 */

#include "electrostatics_magnetostatics/magnetic_non_p3m_methods.hpp"

#ifdef DIPOLES
#include "cells.hpp"
#include "grid.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include "utils/math/int_pow.hpp"

/************************************************************/

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS
   =============================================================================
*/

int Ncut_off_magnetic_dipolar_direct_sum = 0;

/************************************************************/

int magnetic_dipolar_direct_sum_sanity_checks() {
  /* left for the future , at this moment nothing to do */

  return 0;
}

/************************************************************/

double magnetic_dipolar_direct_sum_calculations(int force_flag,
                                                int energy_flag) {
  std::vector<Vector3d> m;
  m.reserve(n_part);

  double u;

  std::vector<Vector3d> pos;
  pos.reserve(n_part);

  std::vector<ParticleForce> f;
  if (force_flag) {
      f.reserve(n_part);
  }

  for (auto const &p : local_cells.particles()) {
    if (p.p.dipm != 0.0) {
      m.emplace_back(p.calc_dip());
      pos.emplace_back(folded_position(p.r.p));
    }
  }

  auto f_kernel = [](Vector3d const& d, Vector3d const& m1, Vector3d const& m2) -> ParticleForce {
    auto const pe2 = m1 * d;
    auto const pe3 = m2 * d;

    auto const r2 = d.norm2();
    auto const r = std::sqrt(r2);
    auto const r5 = r2 * r2 *r;
    auto const r7 = r5 * r2;

    auto const a = 3.0 * (m1 * m2) / r5;
    auto const b = -15.0 * pe2 * pe3 / r7;

    auto const f = (a + b) * d + 3.0 * (pe3 * m1 + pe2 * m2) / r5;
#ifdef ROTATION
    auto const r3 = r2*r;
    auto const t = - m1.cross(m2) / r3 + 3.0 * pe3 * m1.cross(d) / r5;

    return {f, t};
#else
    return f;
#endif
  };

  auto u_kernel = [] (Vector3d const& d, Vector3d const& m1, Vector3d const& m2) -> double {
      auto const r2 = d*d;
      auto const r = sqrt(r2);
      auto const r3 = r2 * r;
      auto const r5 = r3 * r2;

      auto const pe1 = m1 * m2;
      auto const pe2 = m1 * d;
      auto const pe3 = m2 * d;

      // Energy ............................
      return pe1 / r3 - 3.0 * pe2 * pe3 / r5;
  };

  /*now we do the calculations */

  { /* beginning of the area of calculation */
    const Vector3i ncut = Ncut_off_magnetic_dipolar_direct_sum * Vector3i{static_cast<int>(PERIODIC(0)),
                                                                          static_cast<int>(PERIODIC(1)),
                                                                          static_cast<int>(PERIODIC(2))};
    auto const ncut2 = ncut.norm2();
    u = 0;

    for (int i = 0; i < n_part; i++) {
        ParticleForce fi{};

      for (int j = 0; j < n_part; j++) {
          auto const d = pos[i] - pos[j];
          auto const rx = d[0];
          auto const ry = d[1];
          auto const rz = d[2];

          for (int nx = -ncut[0]; nx <= ncut[0]; nx++) {
              auto const rnx = rx + nx * box_l[0];
              for (int ny = -ncut[1]; ny <= ncut[1]; ny++) {
                  auto const rny = ry + ny * box_l[1];
                  for (int nz = -ncut[2]; nz <= ncut[2]; nz++) {
                      if (!(i == j && nx == 0 && ny == 0 && nz == 0)) {
                          if (nx * nx + ny * ny + nz * nz <= ncut2) {
                              auto const rnz = rz + nz * box_l[2];
                              if (energy_flag) {
                                  u += u_kernel({rnx, rny, rnz}, m[i], m[j]);
                              }

                              if (force_flag) {
                                  // force ............................
                                  fi += f_kernel({rnx, rny, rnz}, m[i], m[j]);
                              } /* of force_flag  */
                          }
                      } /* of nx*nx+ny*ny +nz*nz< NCUT*NCUT   and   !(i==j && nx==0 &&
                   ny==0 && nz==0) */
                  }   /* of  for nz */
              }     /* of  for ny  */
          }       /* of  for nx  */
      }
      if(force_flag) {
          f.emplace_back(fi);
      }
    } /* of  j and i  */
  }   /* end of the area of calculation */

  /* set the forces, and torques of the particles within Espresso */
  if (force_flag) {
    int dip_particles2 = 0;
    for (auto &p : local_cells.particles()) {
      if (p.p.dipm != 0.0) {
        p.f.f += coulomb.Dprefactor * f[dip_particles2].f;
        p.f.torque += coulomb.Dprefactor * f[dip_particles2].torque;
      }
        dip_particles2++;
    }
  } /*of if force_flag */

  return 0.5 * coulomb.Dprefactor * u;
}

void mdds_set_params(int n_cut) {
  Ncut_off_magnetic_dipolar_direct_sum = n_cut;

  if ((PERIODIC(0) || PERIODIC(1) || PERIODIC(2)) && Ncut_off_magnetic_dipolar_direct_sum == 0) {
    fprintf(stderr, "Careful:  the number of extra replicas to take into "
                    "account during the direct sum calculation is zero \n");
  }

  if (coulomb.Dmethod != DIPOLAR_DS && coulomb.Dmethod != DIPOLAR_MDLC_DS) {
    set_dipolar_method_local(DIPOLAR_DS);
  }

  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();
}

#endif

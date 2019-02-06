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
  std::vector<double> x, y, z;
  std::vector<double> mx, my, mz;
  std::vector<double> fx, fy, fz;
#ifdef ROTATION
  std::vector<double> tx, ty, tz;
#endif
  int dip_particles, dip_particles2;
  double u;

  if (n_nodes != 1) {
    fprintf(stderr, "error: magnetic Direct Sum is just for one cpu .... \n");
    errexit();
  }
  if (!(force_flag) && !(energy_flag)) {
    fprintf(stderr, " I don't know why you call dawaanr_caclulations with all "
                    "flags zero \n");
    return 0;
  }

  x.resize(n_part);
  y.resize(n_part);
  z.resize(n_part);

  mx.resize(n_part);
  my.resize(n_part);
  mz.resize(n_part);

  if (force_flag) {
    fx.resize(n_part);
    fy.resize(n_part);
    fz.resize(n_part);

#ifdef ROTATION
    tx.resize(n_part);
    ty.resize(n_part);
    tz.resize(n_part);
#endif
  }

  dip_particles = 0;
  for (auto const &p : local_cells.particles()) {
    if (p.p.dipm != 0.0) {
      const Vector3d dip = p.calc_dip();

      mx[dip_particles] = dip[0];
      my[dip_particles] = dip[1];
      mz[dip_particles] = dip[2];

      /* here we wish the coordinates to be folded into the primary box */
      auto const ppos = folded_position(p.r.p);
      x[dip_particles] = ppos[0];
      y[dip_particles] = ppos[1];
      z[dip_particles] = ppos[2];

      if (force_flag) {
        fx[dip_particles] = 0;
        fy[dip_particles] = 0;
        fz[dip_particles] = 0;

#ifdef ROTATION
        tx[dip_particles] = 0;
        ty[dip_particles] = 0;
        tz[dip_particles] = 0;
#endif
      }

      dip_particles++;
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

  /*now we do the calculations */

  { /* beginning of the area of calculation */
    const Vector3i ncut = Ncut_off_magnetic_dipolar_direct_sum * Vector3i{static_cast<int>(PERIODIC(0)),
                                                                          static_cast<int>(PERIODIC(1)),
                                                                          static_cast<int>(PERIODIC(2))};
    auto const ncut2 = ncut.norm2();
    u = 0;

    for (int i = 0; i < dip_particles; i++) {
      for (int j = 0; j < dip_particles; j++) {
        auto const pe1 = mx[i] * mx[j] + my[i] * my[j] + mz[i] * mz[j];
        auto const rx = x[i] - x[j];
        auto const ry = y[i] - y[j];
        auto const rz = z[i] - z[j];

        for (int nx = -ncut[0]; nx <= ncut[0]; nx++) {
          auto const rnx = rx + nx * box_l[0];
          auto const rnx2 = rnx * rnx;
          for (int ny = -ncut[1]; ny <= ncut[1]; ny++) {
            auto const rny = ry + ny * box_l[1];
            auto const rny2 = rny * rny;
            for (int nz = -ncut[2]; nz <= ncut[2]; nz++) {
              if (!(i == j && nx == 0 && ny == 0 && nz == 0)) {
                if (nx * nx + ny * ny + nz * nz <= ncut2) {
                  auto const rnz = rz + nz * box_l[2];
                  if(energy_flag)
                    {
                        auto const r2 = rnx2 + rny2 + rnz * rnz;
                        auto const r = sqrt(r2);
                        auto const r3 = r2 * r;
                        auto const r5 = r3 * r2;

                        auto const pe2 = mx[i] * rnx + my[i] * rny + mz[i] * rnz;
                        auto const pe3 = mx[j] * rnx + my[j] * rny + mz[j] * rnz;

                        // Energy ............................

                        u += pe1 / r3 - 3.0 * pe2 * pe3 / r5;
                    }

                  if (force_flag) {
                    // force ............................
                          auto const f = f_kernel({rnx, rny, rnz}, {mx[i], my[i], mz[i]}, {mx[j], my[j], mz[j]});

                          fx[i] += f.f[0];
                          fy[i] += f.f[1];
                          fz[i] += f.f[2];
#ifdef ROTATION
                    // torque ............................
                          tx[i] += f.torque[0];
                          ty[i] += f.torque[1];
                          tz[i] += f.torque[2];
#endif
                  } /* of force_flag  */
                }
              } /* of nx*nx+ny*ny +nz*nz< NCUT*NCUT   and   !(i==j && nx==0 &&
                   ny==0 && nz==0) */
            }   /* of  for nz */
          }     /* of  for ny  */
        }       /* of  for nx  */
      }

    } /* of  j and i  */
  }   /* end of the area of calculation */

  /* set the forces, and torques of the particles within Espresso */
  if (force_flag) {

    dip_particles2 = 0;

    for (auto &p : local_cells.particles()) {
      if (p.p.dipm != 0.0) {

        p.f.f[0] += coulomb.Dprefactor * fx[dip_particles2];
        p.f.f[1] += coulomb.Dprefactor * fy[dip_particles2];
        p.f.f[2] += coulomb.Dprefactor * fz[dip_particles2];

#ifdef ROTATION
        p.f.torque[0] += coulomb.Dprefactor * tx[dip_particles2];
        p.f.torque[1] += coulomb.Dprefactor * ty[dip_particles2];
        p.f.torque[2] += coulomb.Dprefactor * tz[dip_particles2];
#endif
        dip_particles2++;
      }
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

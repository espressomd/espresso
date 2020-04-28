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

#include "electrostatics_magnetostatics/magnetic_non_p3m_methods.hpp"

#ifdef DIPOLES
#include "cells.hpp"
#include "communication.hpp"
#include "dipole.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/constants.hpp>

double calc_dipole_dipole_ia(Particle &p1, Utils::Vector3d const &dip1,
                             Particle &p2, bool force_flag) {

  // Cache dipole moment
  auto const dip2 = p2.calc_dip();

  // Distance between particles
  auto const dr = get_mi_vector(p1.r.p, p2.r.p, box_geo);

  // Powers of distance
  auto const r2 = dr.norm2();
  auto const r = sqrt(r2);
  auto const r3 = r2 * r;
  auto const r5 = r3 * r2;
  auto const r7 = r5 * r2;

  // Dot products
  auto const pe1 = dip1 * dip2;
  auto const pe2 = dip1 * dr;
  auto const pe3 = dip2 * dr;
  auto const pe4 = 3.0 / r5;

  // Energy
  auto const u = dipole.prefactor * (pe1 / r3 - pe4 * pe2 * pe3);

  // Forces, if requested
  if (force_flag) {
    auto const a = pe4 * pe1;
    auto const b = -15.0 * pe2 * pe3 / r7;
    auto const ab = a + b;
    auto const cc = pe4 * pe3;
    auto const dd = pe4 * pe2;

    // Forces
    auto const ff = ab * dr + cc * dip1 + dd * dip2;
    p1.f.f += dipole.prefactor * ff;
    p2.f.f -= dipole.prefactor * ff;

#ifdef ROTATION
    // Torques
    auto const aa = vector_product(dip1, dip2);
    auto const b1 = vector_product(dip1, dr);
    auto const b2 = vector_product(dip2, dr);
    p1.f.torque += dipole.prefactor * (-aa / r3 + b1 * cc);
    p2.f.torque += dipole.prefactor * (aa / r3 + b2 * dd);
#endif
  }

  // Return energy
  return u;
}

/* =============================================================================
                  DAWAANR => DIPOLAR ALL WITH ALL AND NO REPLICA
   =============================================================================
*/

double dawaanr_calculations(bool force_flag, bool energy_flag,
                            const ParticleRange &particles) {

  if (n_nodes != 1) {
    fprintf(stderr, "error: DAWAANR is just for one cpu...\n");
    errexit();
  }
  if (!(force_flag) && !(energy_flag)) {
    fprintf(stderr, "I don't know why you call dawaanr_calculations() "
                    "with all flags zero.\n");
    return 0;
  }

  // Variable to sum up the energy
  double u = 0;

  auto parts = particles;

  // Iterate over all cells
  for (auto it = parts.begin(), end = parts.end(); it != end; ++it) {
    // If the particle has no dipole moment, ignore it
    if (it->p.dipm == 0.0)
      continue;

    const Utils::Vector3d dip1 = it->calc_dip();
    auto jt = it;
    /* Skip diagonal */
    ++jt;
    for (; jt != end; ++jt) {
      // If the particle has no dipole moment, ignore it
      if (jt->p.dipm == 0.0)
        continue;
      // Calculate energy and/or force between the particles
      u += calc_dipole_dipole_ia(*it, dip1, *jt, force_flag);
    }
  }

  // Return energy
  return u;
}

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS
   =============================================================================
*/

int Ncut_off_magnetic_dipolar_direct_sum = 0;

int magnetic_dipolar_direct_sum_sanity_checks() {
  /* left for the future, at this moment nothing to do */

  return 0;
}

double
magnetic_dipolar_direct_sum_calculations(bool force_flag, bool energy_flag,
                                         ParticleRange const &particles) {
  std::vector<double> x, y, z;
  std::vector<double> mx, my, mz;
  std::vector<double> fx, fy, fz;
#ifdef ROTATION
  std::vector<double> tx, ty, tz;
#endif
  double u;

  if (n_nodes != 1) {
    fprintf(stderr, "error: magnetic Direct Sum is just for one cpu...\n");
    errexit();
  }
  if (!(force_flag) && !(energy_flag)) {
    fprintf(stderr, "I don't know why you call magnetic_dipolar_direct_sum_"
                    "calculations() with all flags zero\n");
    return 0;
  }

  auto const n_part = particles.size();

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

  int dip_particles = 0;
  for (auto const &p : particles) {
    if (p.p.dipm != 0.0) {
      const Utils::Vector3d dip = p.calc_dip();

      mx[dip_particles] = dip[0];
      my[dip_particles] = dip[1];
      mz[dip_particles] = dip[2];

      /* here we wish the coordinates to be folded into the primary box */
      auto const ppos = folded_position(p.r.p, box_geo);
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

  /*now we do the calculations */

  { /* beginning of the area of calculation */
    int NCUT[3], NCUT2;

    for (int i = 0; i < 3; i++) {
      NCUT[i] = Ncut_off_magnetic_dipolar_direct_sum;
      if (box_geo.periodic(i) == 0) {
        NCUT[i] = 0;
      }
    }
    NCUT2 = Ncut_off_magnetic_dipolar_direct_sum *
            Ncut_off_magnetic_dipolar_direct_sum;

    u = 0;

    for (int i = 0; i < dip_particles; i++) {
      for (int j = 0; j < dip_particles; j++) {
        auto const pe1 = mx[i] * mx[j] + my[i] * my[j] + mz[i] * mz[j];
        auto const rx = x[i] - x[j];
        auto const ry = y[i] - y[j];
        auto const rz = z[i] - z[j];

        for (int nx = -NCUT[0]; nx <= NCUT[0]; nx++) {
          auto const rnx = rx + nx * box_geo.length()[0];
          auto const rnx2 = rnx * rnx;
          for (int ny = -NCUT[1]; ny <= NCUT[1]; ny++) {
            auto const rny = ry + ny * box_geo.length()[1];
            auto const rny2 = rny * rny;
            for (int nz = -NCUT[2]; nz <= NCUT[2]; nz++) {
              if (!(i == j && nx == 0 && ny == 0 && nz == 0)) {
                if (nx * nx + ny * ny + nz * nz <= NCUT2) {
                  auto const rnz = rz + nz * box_geo.length()[2];
                  auto const r2 = rnx2 + rny2 + rnz * rnz;
                  auto const r = sqrt(r2);
                  auto const r3 = r2 * r;
                  auto const r5 = r3 * r2;
                  auto const r7 = r5 * r2;

                  auto const pe2 = mx[i] * rnx + my[i] * rny + mz[i] * rnz;
                  auto const pe3 = mx[j] * rnx + my[j] * rny + mz[j] * rnz;

                  // Energy ............................

                  u += pe1 / r3 - 3.0 * pe2 * pe3 / r5;

                  if (force_flag) {
                    double a, b, c, d;
                    // force ............................
                    a = mx[i] * mx[j] + my[i] * my[j] + mz[i] * mz[j];
                    a = 3.0 * a / r5;
                    b = -15.0 * pe2 * pe3 / r7;
                    c = 3.0 * pe3 / r5;
                    d = 3.0 * pe2 / r5;

                    fx[i] += (a + b) * rnx + c * mx[i] + d * mx[j];
                    fy[i] += (a + b) * rny + c * my[i] + d * my[j];
                    fz[i] += (a + b) * rnz + c * mz[i] + d * mz[j];

#ifdef ROTATION
                    // torque ............................
                    c = 3.0 / r5 * pe3;
                    auto const ax = my[i] * mz[j] - my[j] * mz[i];
                    auto const ay = mx[j] * mz[i] - mx[i] * mz[j];
                    auto const az = mx[i] * my[j] - mx[j] * my[i];

                    auto const bx = my[i] * rnz - rny * mz[i];
                    auto const by = rnx * mz[i] - mx[i] * rnz;
                    auto const bz = mx[i] * rny - rnx * my[i];

                    tx[i] += -ax / r3 + bx * c;
                    ty[i] += -ay / r3 + by * c;
                    tz[i] += -az / r3 + bz * c;
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

  /* set the forces, and torques of the particles within ESPResSo */
  if (force_flag) {

    int dip_particles2 = 0;

    for (auto &p : particles) {
      if (p.p.dipm != 0.0) {

        p.f.f[0] += dipole.prefactor * fx[dip_particles2];
        p.f.f[1] += dipole.prefactor * fy[dip_particles2];
        p.f.f[2] += dipole.prefactor * fz[dip_particles2];

#ifdef ROTATION
        p.f.torque[0] += dipole.prefactor * tx[dip_particles2];
        p.f.torque[1] += dipole.prefactor * ty[dip_particles2];
        p.f.torque[2] += dipole.prefactor * tz[dip_particles2];
#endif
        dip_particles2++;
      }
    }
  } /*of if force_flag */

  return 0.5 * dipole.prefactor * u;
}

int dawaanr_set_params() {
  if (n_nodes > 1) {
    runtimeErrorMsg() << "MPI parallelization not supported by "
                      << "DipolarDirectSumCpu.";
    return ES_ERROR;
  }
  if (dipole.method != DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA) {
    Dipole::set_method_local(DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA);
  }
  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();

  return ES_OK;
}

int mdds_set_params(int n_cut) {
  if (n_nodes > 1) {
    runtimeErrorMsg() << "MPI parallelization not supported by "
                      << "DipolarDirectSumWithReplicaCpu.";
    return ES_ERROR;
  }

  Ncut_off_magnetic_dipolar_direct_sum = n_cut;

  if (Ncut_off_magnetic_dipolar_direct_sum == 0) {
    fprintf(stderr, "Careful: the number of extra replicas to take into "
                    "account during the direct sum calculation is zero\n");
  }

  if (dipole.method != DIPOLAR_DS && dipole.method != DIPOLAR_MDLC_DS) {
    Dipole::set_method_local(DIPOLAR_DS);
  }

  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();
  return ES_OK;
}

#endif

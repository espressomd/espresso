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

#include "config.hpp"

#ifdef DIPOLES

#include "electrostatics_magnetostatics/magnetic_non_p3m_methods.hpp"

#include "electrostatics_magnetostatics/common.hpp"
#include "electrostatics_magnetostatics/dipole.hpp"

#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <cassert>
#include <cmath>
#include <cstdio>
#include <stdexcept>
#include <vector>

/**
 * Calculate dipolar energy and optionally force between two particles.
 * @param[in,out] p1          First particle
 * @param[in]     dip1        Cached dipole moment of the first particle
 * @param[in,out] p2          Second particle
 * @param[in]     force_flag  If true, update the particle forces and torques
 */
static double calc_dipole_dipole_ia(Particle &p1, Utils::Vector3d const &dip1,
                                    Particle &p2, bool force_flag) {

  // Cache dipole moment
  auto const dip2 = p2.calc_dip();

  // Distance between particles
  auto const dr = box_geo.get_mi_vector(p1.r.p, p2.r.p);

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
  auto const energy = dipole.prefactor * (pe1 / r3 - pe4 * pe2 * pe3);

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

    // Torques
    auto const aa = vector_product(dip1, dip2);
    auto const b1 = vector_product(dip1, dr);
    auto const b2 = vector_product(dip2, dr);
    p1.f.torque += dipole.prefactor * (-aa / r3 + b1 * cc);
    p2.f.torque += dipole.prefactor * (aa / r3 + b2 * dd);
  }

  return energy;
}

/* =============================================================================
                  DAWAANR => DIPOLAR ALL WITH ALL AND NO REPLICA
   =============================================================================
*/

double dawaanr_calculations(bool force_flag, bool energy_flag,
                            const ParticleRange &particles) {

  assert(n_nodes == 1);
  assert(force_flag || energy_flag);

  double energy = 0.0;
  // Iterate over all particles
  for (auto it = particles.begin(), end = particles.end(); it != end; ++it) {
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
      energy += calc_dipole_dipole_ia(*it, dip1, *jt, force_flag);
    }
  }

  return energy;
}

/* =============================================================================
                  DIRECT SUM FOR MAGNETIC SYSTEMS
   =============================================================================
*/

static int mdds_n_replica = 0;

int mdds_get_n_replica() { return mdds_n_replica; }

void sanity_checks(int n_replica) {
  if (box_geo.periodic(0) and box_geo.periodic(1) and box_geo.periodic(2) and
      n_replica == 0) {
    throw std::runtime_error("Dipolar direct sum with replica does not "
                             "support a periodic system with zero replica.");
  }
}

void mdds_sanity_checks() { sanity_checks(mdds_n_replica); }

double
magnetic_dipolar_direct_sum_calculations(bool force_flag, bool energy_flag,
                                         ParticleRange const &particles) {

  assert(n_nodes == 1);
  assert(force_flag || energy_flag);

  std::vector<double> x, y, z;
  std::vector<double> mx, my, mz;
  std::vector<double> fx, fy, fz;
  std::vector<double> tx, ty, tz;

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
    tx.resize(n_part);
    ty.resize(n_part);
    tz.resize(n_part);
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
        tx[dip_particles] = 0;
        ty[dip_particles] = 0;
        tz[dip_particles] = 0;
      }

      dip_particles++;
    }
  }

  /* energy calculation */
  double energy = 0.;

  int NCUT[3];
  for (int i = 0; i < 3; i++) {
    NCUT[i] = box_geo.periodic(i) ? mdds_n_replica : 0;
  }
  auto const NCUT2 = Utils::sqr(mdds_n_replica);

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
            if (!(i == j && nx == 0 && ny == 0 && nz == 0) and
                (nx * nx + ny * ny + nz * nz <= NCUT2)) {
              auto const rnz = rz + nz * box_geo.length()[2];
              auto const r2 = rnx2 + rny2 + rnz * rnz;
              auto const r = sqrt(r2);
              auto const r3 = r2 * r;
              auto const r5 = r3 * r2;
              auto const r7 = r5 * r2;

              auto const pe2 = mx[i] * rnx + my[i] * rny + mz[i] * rnz;
              auto const pe3 = mx[j] * rnx + my[j] * rny + mz[j] * rnz;
              auto const pe4 = 3.0 / r5;

              // Energy
              energy += pe1 / r3 - pe4 * pe2 * pe3;

              if (force_flag) {
                // Forces
                auto const a = pe4 * pe1;
                auto const b = -15.0 * pe2 * pe3 / r7;
                auto const c = pe4 * pe3;
                auto const d = pe4 * pe2;

                fx[i] += (a + b) * rnx + c * mx[i] + d * mx[j];
                fy[i] += (a + b) * rny + c * my[i] + d * my[j];
                fz[i] += (a + b) * rnz + c * mz[i] + d * mz[j];

                // Torques
                auto const ax = my[i] * mz[j] - my[j] * mz[i];
                auto const ay = mx[j] * mz[i] - mx[i] * mz[j];
                auto const az = mx[i] * my[j] - mx[j] * my[i];

                auto const bx = my[i] * rnz - rny * mz[i];
                auto const by = rnx * mz[i] - mx[i] * rnz;
                auto const bz = mx[i] * rny - rnx * my[i];

                tx[i] += -ax / r3 + bx * c;
                ty[i] += -ay / r3 + by * c;
                tz[i] += -az / r3 + bz * c;
              } /* if force_flag  */
            }   /* if distance criterion */
          }     /* for nz */
        }       /* for ny */
      }         /* for nx */
    }           /* for j */
  }             /* for i */

  /* update particle forces and torques */
  if (force_flag) {

    dip_particles = 0;

    for (auto &p : particles) {
      if (p.p.dipm != 0.0) {

        p.f.f[0] += dipole.prefactor * fx[dip_particles];
        p.f.f[1] += dipole.prefactor * fy[dip_particles];
        p.f.f[2] += dipole.prefactor * fz[dip_particles];

        p.f.torque[0] += dipole.prefactor * tx[dip_particles];
        p.f.torque[1] += dipole.prefactor * ty[dip_particles];
        p.f.torque[2] += dipole.prefactor * tz[dip_particles];

        dip_particles++;
      }
    }
  } /* if force_flag */

  return 0.5 * dipole.prefactor * energy;
}

void dawaanr_set_params() {
  if (n_nodes > 1) {
    throw std::runtime_error(
        "MPI parallelization not supported by DipolarDirectSumCpu.");
  }
  if (dipole.method != DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA) {
    Dipole::set_method_local(DIPOLAR_ALL_WITH_ALL_AND_NO_REPLICA);
  }
  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();
}

void mdds_set_params(int n_replica) {
  if (n_nodes > 1) {
    throw std::runtime_error(
        "MPI parallelization not supported by DipolarDirectSumWithReplicaCpu.");
  }
  if (n_replica < 0) {
    throw std::runtime_error("Dipolar direct sum requires n_replica >= 0.");
  }
  sanity_checks(n_replica);
  if (n_replica == 0) {
    fprintf(stderr, "Careful: the number of extra replicas to take into "
                    "account during the direct sum calculation is zero\n");
  }

  mdds_n_replica = n_replica;

  if (dipole.method != DIPOLAR_DS && dipole.method != DIPOLAR_MDLC_DS) {
    Dipole::set_method_local(DIPOLAR_DS);
  }

  // also necessary on 1 CPU, does more than just broadcasting
  mpi_bcast_coulomb_params();
}

#endif

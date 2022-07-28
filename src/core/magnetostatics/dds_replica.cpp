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

#include "config.hpp"

#ifdef DIPOLES

#include "magnetostatics/dds_replica.hpp"

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

void DipolarDirectSumWithReplica::sanity_checks_node_grid() const {
  if (box_geo.periodic(0) and box_geo.periodic(1) and box_geo.periodic(2) and
      n_replica == 0) {
    throw std::runtime_error("Dipolar direct sum with replica does not "
                             "support a periodic system with zero replica.");
  }
}

double
DipolarDirectSumWithReplica::kernel(bool force_flag, bool energy_flag,
                                    ParticleRange const &particles) const {

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
    if (p.dipm() != 0.0) {
      const Utils::Vector3d dip = p.calc_dip();

      mx[dip_particles] = dip[0];
      my[dip_particles] = dip[1];
      mz[dip_particles] = dip[2];

      /* here we wish the coordinates to be folded into the primary box */
      auto const ppos = folded_position(p.pos(), box_geo);
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
    NCUT[i] = box_geo.periodic(i) ? n_replica : 0;
  }
  auto const NCUT2 = Utils::sqr(n_replica);

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
      if (p.dipm() != 0.) {
        auto &force = p.force();
        auto &torque = p.torque();

        force[0] += prefactor * fx[dip_particles];
        force[1] += prefactor * fy[dip_particles];
        force[2] += prefactor * fz[dip_particles];

        torque[0] += prefactor * tx[dip_particles];
        torque[1] += prefactor * ty[dip_particles];
        torque[2] += prefactor * tz[dip_particles];

        dip_particles++;
      }
    }
  } /* if force_flag */

  return 0.5 * prefactor * energy;
}

DipolarDirectSumWithReplica::DipolarDirectSumWithReplica(double prefactor,
                                                         int n_replica)
    : prefactor{prefactor}, n_replica{n_replica} {
  if (n_nodes > 1) {
    throw std::runtime_error(
        "MPI parallelization not supported by DipolarDirectSumWithReplicaCpu.");
  }
  if (prefactor <= 0.) {
    throw std::domain_error("Parameter 'prefactor' must be > 0");
  }
  if (n_replica < 0) {
    throw std::domain_error("Parameter 'n_replica' must be >= 0");
  }
  sanity_checks();
  if (n_replica == 0) {
    fprintf(stderr, "Careful: the number of extra replicas to take into "
                    "account during the direct sum calculation is zero\n");
  }
}

#endif // DIPOLES

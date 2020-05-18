/*
 * Copyright (C) 2010-2019 The ESPResSo project
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
/// \file
/// \brief Main of the Bayreuth Immersed-Boundary implementation

#include "virtual_sites/lb_inertialess_tracers.hpp"

#ifdef VIRTUAL_SITES_INERTIALESS_TRACERS
#include "Particle.hpp"
#include "cells.hpp"
#include "grid.hpp"
#include "grid_based_algorithms/lb.hpp"
#include "grid_based_algorithms/lb_boundaries.hpp"
#include "grid_based_algorithms/lb_interface.hpp"
#include "integrate.hpp"
#include "lb_inertialess_tracers_cuda_interface.hpp"

#include <utils/math/sqr.hpp>

// ****** Functions for internal use ********

void CoupleIBMParticleToFluid(Particle *p);
void ParticleVelocitiesFromLB_CPU();
bool IsHalo(int indexCheck);
void GetIBMInterpolatedVelocity(const Utils::Vector3d &p, double *v,
                                double *forceAdded);

// ***** Internal variables ******

bool *isHaloCache = nullptr;

namespace {
bool in_local_domain(Utils::Vector3d const &pos) {
  auto const lblattice = lb_lbfluid_get_lattice();
  auto const my_left = local_geo.my_left();
  auto const my_right = local_geo.my_right();

  return (pos[0] >= my_left[0] - 0.5 * lblattice.agrid &&
          pos[0] < my_right[0] + 0.5 * lblattice.agrid &&
          pos[1] >= my_left[1] - 0.5 * lblattice.agrid &&
          pos[1] < my_right[1] + 0.5 * lblattice.agrid &&
          pos[2] >= my_left[2] - 0.5 * lblattice.agrid &&
          pos[2] < my_right[2] + 0.5 * lblattice.agrid);
}
} // namespace

/** Put the calculated force stored on the ibm particles into the fluid by
 *  updating the @ref lbfields structure.
 *  Called from the integration loop right after the forces have been
 *  calculated.
 */
void IBM_ForcesIntoFluid_CPU() {
  // Update the forces on the ghost particles
  cell_structure.ghosts_update(Cells::DATA_PART_FORCE);

  // Loop over local cells
  for (auto &p : cell_structure.local_particles()) {
    if (p.p.is_virtual)
      CoupleIBMParticleToFluid(&p);
  }

  for (auto &p : cell_structure.ghost_particles()) {
    // for ghost particles we have to check if they lie
    // in the range of the local lattice nodes
    if (in_local_domain(p.r.p)) {
      if (p.p.is_virtual)
        CoupleIBMParticleToFluid(&p);
    }
  }
}

/** Interpolate LB velocity at the particle positions and propagate the
 *  particles.
 *  Called from the integration loop right after the LB update.
 */
void IBM_UpdateParticlePositions(ParticleRange particles) {
  // Get velocities
  if (lattice_switch == ActiveLB::CPU)
    ParticleVelocitiesFromLB_CPU();
#ifdef CUDA
  if (lattice_switch == ActiveLB::GPU)
    ParticleVelocitiesFromLB_GPU(particles);
#endif

  // Do update: Euler
  const double skin2 = Utils::sqr(0.5 * skin);
  // Loop over particles in local cells
  for (auto &p : particles) {
    if (p.p.is_virtual) {
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & 2))
#endif
        p.r.p[0] = p.r.p[0] + p.m.v[0] * time_step;
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & 4))
#endif
        p.r.p[1] = p.r.p[1] + p.m.v[1] * time_step;
#ifdef EXTERNAL_FORCES
      if (!(p.p.ext_flag & 8))
#endif
        p.r.p[2] = p.r.p[2] + p.m.v[2] * time_step;

      // Check if the particle might have crossed a box border
      const double dist2 = (p.r.p - p.l.p_old).norm2();
      if (dist2 > skin2) {
        cell_structure.set_resort_particles(Cells::RESORT_LOCAL);
      }
    }
  }
}

/** Put the momentum of a given particle into the LB fluid. */
void CoupleIBMParticleToFluid(Particle *p) {
  // Convert units from MD to LB
  double delta_j[3];
  delta_j[0] = p->f.f[0] * lbpar.tau * lbpar.tau / lbpar.agrid;
  delta_j[1] = p->f.f[1] * lbpar.tau * lbpar.tau / lbpar.agrid;
  delta_j[2] = p->f.f[2] * lbpar.tau * lbpar.tau / lbpar.agrid;

  // Get indices and weights of affected nodes using discrete delta function
  Utils::Vector<std::size_t, 8> node_index{};
  Utils::Vector6d delta{};
  lblattice.map_position_to_lattice(p->r.p, node_index, delta);

  // Loop over all affected nodes
  for (int z = 0; z < 2; z++) {
    for (int y = 0; y < 2; y++) {
      for (int x = 0; x < 2; x++) {
        // Do not put force into a halo node
        if (!IsHalo(node_index[(z * 2 + y) * 2 + x])) {
          // Add force into the lbfields structure
          auto &local_f =
              lbfields[node_index[(z * 2 + y) * 2 + x]].force_density;

          local_f[0] += delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2] *
                        delta_j[0];
          local_f[1] += delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2] *
                        delta_j[1];
          local_f[2] += delta[3 * x + 0] * delta[3 * y + 1] * delta[3 * z + 2] *
                        delta_j[2];
        }
      }
    }
  }
}

/** Calculate the LB fluid velocity at a particle position.
 *  Very similar to the velocity interpolation done in standard ESPResSo,
 *  except that we add the f/2 contribution, cf. @cite guo02a.
 *  The fluid velocity is obtained by linear interpolation,
 *  cf. eq. (11) in @cite ahlrichs99a.
 */
void GetIBMInterpolatedVelocity(const Utils::Vector3d &pos, double *v,
                                double *forceAdded) {
  /* determine elementary lattice cell surrounding the particle
     and the relative position of the particle in this cell */
  Utils::Vector<std::size_t, 8> node_index{};
  Utils::Vector6d delta{};
  lblattice.map_position_to_lattice(pos, node_index, delta);

  Utils::Vector3d interpolated_u = {};
  // This for the f/2 contribution to the velocity
  forceAdded[0] = forceAdded[1] = forceAdded[2] = 0;

  for (int z = 0; z < 2; z++) {
    for (int y = 0; y < 2; y++) {
      for (int x = 0; x < 2; x++) {
        auto const index = node_index[(z * 2 + y) * 2 + x];
        const auto &f = lbfields[index].force_density_buf;

        double local_density;
        Utils::Vector3d local_j;

        // This can be done more easily without copying the code twice.
        // We probably can even set the boundary velocity directly.
#ifdef LB_BOUNDARIES
        if (lbfields[index].boundary) {
          local_density = lbpar.density;
          local_j = lbpar.density *
                    (*LBBoundaries::lbboundaries[lbfields[index].boundary - 1])
                        .velocity();
        } else
#endif
        {
          auto const modes = lb_calc_modes(index, lbfluid);
          local_density = lbpar.density + modes[0];

          // Add the +f/2 contribution!!
          local_j[0] = modes[1] + f[0] / 2;
          local_j[1] = modes[2] + f[1] / 2;
          local_j[2] = modes[3] + f[2] / 2;

          // Keep track of the forces that we added to the fluid.
          // This is necessary for communication because this part is executed
          // for real and ghost particles.
          // Later on we sum the real and ghost contributions.
          const double fExt[3] = {
              lbpar.ext_force_density[0] * pow(lbpar.agrid, 2) * lbpar.tau *
                  lbpar.tau,
              lbpar.ext_force_density[1] * pow(lbpar.agrid, 2) * lbpar.tau *
                  lbpar.tau,
              lbpar.ext_force_density[2] * pow(lbpar.agrid, 2) * lbpar.tau *
                  lbpar.tau};

          forceAdded[0] += delta[3 * x + 0] * delta[3 * y + 1] *
                           delta[3 * z + 2] * (f[0] - fExt[0]) / 2 /
                           (local_density);
          forceAdded[1] += delta[3 * x + 0] * delta[3 * y + 1] *
                           delta[3 * z + 2] * (f[1] - fExt[1]) / 2 /
                           (local_density);
          forceAdded[2] += delta[3 * x + 0] * delta[3 * y + 1] *
                           delta[3 * z + 2] * (f[2] - fExt[2]) / 2 /
                           (local_density);
        }

        // Interpolate velocity
        interpolated_u[0] += delta[3 * x + 0] * delta[3 * y + 1] *
                             delta[3 * z + 2] * local_j[0] / (local_density);
        interpolated_u[1] += delta[3 * x + 0] * delta[3 * y + 1] *
                             delta[3 * z + 2] * local_j[1] / (local_density);
        interpolated_u[2] += delta[3 * x + 0] * delta[3 * y + 1] *
                             delta[3 * z + 2] * local_j[2] / (local_density);
      }
    }
  }

  v[0] = interpolated_u[0];
  v[1] = interpolated_u[1];
  v[2] = interpolated_u[2];

  v[0] *= lbpar.agrid / lbpar.tau;
  v[1] *= lbpar.agrid / lbpar.tau;
  v[2] *= lbpar.agrid / lbpar.tau;
}

/** Build a cache structure which contains a flag for each LB node whether
 * that node is a halo node or not.
 */
bool IsHalo(const int indexCheck) {
  // First call --> build cache
  if (isHaloCache == nullptr) {
    isHaloCache = new bool[lblattice.halo_grid_volume];
    // Assume everything is a halo and correct in the next step
    for (int i = 0; i < lblattice.halo_grid_volume; i++)
      isHaloCache[i] = true;
    // Loop through and check where indexCheck occurs
    int index = lblattice.halo_offset;
    for (int z = 1; z <= lblattice.grid[2]; z++) {
      for (int y = 1; y <= lblattice.grid[1]; y++) {
        for (int x = 1; x <= lblattice.grid[0]; x++) {
          isHaloCache[index] = false;
          ++index;
        }
        index += 2; /* skip halo region */
      }
      index += 2 * lblattice.halo_grid[0]; /* skip halo region */
    }
  }

  // Return
  return isHaloCache[indexCheck];
}

/** Get particle velocities from LB and set the velocity field in the
 * particles data structure.
 */
void ParticleVelocitiesFromLB_CPU() {
  // Loop over particles in local cells.
  // Here all contributions are included: velocity, external force and
  // particle force.
  for (auto &p : cell_structure.local_particles()) {
    if (p.p.is_virtual) {
      double dummy[3];
      // Get interpolated velocity and store in the force (!) field
      // for later communication (see below)
      GetIBMInterpolatedVelocity(p.r.p, p.f.f.data(), dummy);
    }
  }

  // Loop over particles in ghost cells
  // Here we only add the particle forces stemming from the ghosts
  for (auto &p : cell_structure.ghost_particles()) {
    // This criterion include the halo on the left, but excludes the halo on
    // the right
    // Try if we have to use *1.5 on the right
    if (in_local_domain(p.r.p)) {
      if (p.p.is_virtual) {
        double dummy[3];
        double force[3] = {0, 0,
                           0}; // The force stemming from the ghost particle
        GetIBMInterpolatedVelocity(p.r.p, dummy, force);

        // Rescale and store in the force field of the particle (for
        // communication, see below)
        p.f.f[0] = force[0] * lbpar.agrid / lbpar.tau;
        p.f.f[1] = force[1] * lbpar.agrid / lbpar.tau;
        p.f.f[2] = force[2] * lbpar.agrid / lbpar.tau;
      } else {
        p.f.f[0] = p.f.f[1] = p.f.f[2] = 0;
      } // Reset, necessary because we add all forces below. Also needs to
        // be done for the real particles!

    } else {
      p.f.f[0] = p.f.f[1] = p.f.f[2] = 0;
    } // Reset, necessary because we add all forces below
  }

  // Now the local particles contain a velocity (stored in the force field)
  // and the ghosts contain the rest of the velocity in their respective force
  // fields.
  // We need to add these. Since we have stored them in the force, not the
  // velocity fields, we can use the standard force communicator and then
  // transfer to the velocity afterwards.
  // Note that this overwrites the actual force which would be a problem for
  // real particles.
  // This could be solved by keeping a backup of the local forces before this
  // operation is attempted.
  cell_structure.ghosts_reduce_forces();

  // Transfer to velocity field
  for (auto &p : cell_structure.local_particles()) {
    if (p.p.is_virtual) {
      p.m.v[0] = p.f.f[0];
      p.m.v[1] = p.f.f[1];
      p.m.v[2] = p.f.f[2];
    }
  }
}
#endif

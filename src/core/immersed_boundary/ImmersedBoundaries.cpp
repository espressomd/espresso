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

#include "ImmersedBoundaries.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"
#include "cells.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

#include <utils/constants.hpp>

/** Volume conservation.
 *  Calculate volumes, volume force and add it to each virtual particle.
 *  This function is called from integrate_vv
 */
void ImmersedBoundaries::volume_conservation() {
  if (VolumeInitDone && !BoundariesFound) {
    return;
  }
  // Calculate volumes
  calc_volumes();

  calc_volume_force();

  // Center-of-mass output
  // if ( numWriteCOM > 0 )
  //    IBM_CalcCentroids(frameNum, simTime);
}

/** Initialize volume conservation */
void ImmersedBoundaries::init_volume_conservation() {

  // Check since this function is called at the start of every integrate loop
  // Also check if volume has been set due to reading of a checkpoint
  if (!VolumeInitDone) {

    // Calculate volumes
    calc_volumes();

    // numWriteCOM = 0;

    // Loop through all bonded interactions and check if we need to set the
    // reference volume
    for (auto &bonded_ia_param : bonded_ia_params) {
      if (bonded_ia_param.type == BONDED_IA_IBM_VOLUME_CONSERVATION) {
        // This check is important because InitVolumeConservation may be called
        // accidentally during the integration. Then we must not reset the
        // reference
        BoundariesFound = true;
        if (bonded_ia_param.p.ibmVolConsParameters.volRef == 0) {
          const int softID = bonded_ia_param.p.ibmVolConsParameters.softID;
          bonded_ia_param.p.ibmVolConsParameters.volRef =
              VolumesCurrent[softID];
        }
      }
    }
  }

  VolumeInitDone = true;
}

/** Set parameters of volume conservation */
int ImmersedBoundaries::volume_conservation_set_params(const int bond_type,
                                                       const int softID,
                                                       const double kappaV) {
  // Create bond
  make_bond_type_exist(bond_type);

  // General bond parameters
  bonded_ia_params[bond_type].type = BONDED_IA_IBM_VOLUME_CONSERVATION;
  bonded_ia_params[bond_type].num =
      0; // This means that Espresso requires one bond partner. Here we simply
         // ignore it, but Espresso cannot handle 0.

  // Specific stuff
  if (softID > MaxNumIBM) {
    printf("Error: softID (%d) is larger than MaxNumIBM (%d)\n", softID,
           MaxNumIBM);
    return ES_ERROR;
  }
  if (softID < 0) {
    printf("Error: softID (%d) must be non-negative\n", softID);
    return ES_ERROR;
  }

  bonded_ia_params[bond_type].p.ibmVolConsParameters.softID = softID;
  bonded_ia_params[bond_type].p.ibmVolConsParameters.kappaV = kappaV;
  bonded_ia_params[bond_type].p.ibmVolConsParameters.volRef = 0;
  // NOTE: We cannot compute the reference volume here because not all
  // interactions are setup and thus we do not know which triangles belong to
  // this softID Calculate it later in the init function

  // Communicate this to whoever is interested
  mpi_bcast_ia_params(bond_type, -1);

  return ES_OK;
}

/** Calculate partial volumes on all compute nodes and call MPI to sum up */
void ImmersedBoundaries::calc_volumes() {

  // Partial volumes for each soft particle, to be summed up
  std::vector<double> tempVol(MaxNumIBM);

  // Loop over all particles on local node
  for (int c = 0; c < local_cells.n; c++) {
    const Cell *const cell = local_cells.cell[c];

    for (int i = 0; i < cell->n; i++) {
      Particle &p1 = cell->part[i];

      // Check if particle has a BONDED_IA_IBM_TRIEL and a
      // BONDED_IA_IBM_VOLUME_CONSERVATION Basically this loops over all
      // triangles, not all particles First round to check for volume
      // conservation and virtual Loop over all bonds of this particle Actually
      // j loops over the bond-list, i.e. the bond partners (see
      // particle_data.hpp)
      int softID = -1;
      int j = 0;
      while (j < p1.bl.n) {
        const int type_num = p1.bl.e[j];
        const Bonded_ia_parameters &iaparams = bonded_ia_params[type_num];
        const int type = iaparams.type;
        if (type == BONDED_IA_IBM_VOLUME_CONSERVATION) {
          if (p1.p.is_virtual)
            softID = iaparams.p.ibmVolConsParameters.softID;
          else {
            printf("Error. Encountered non-virtual particle with "
                   "VOLUME_CONSERVATION_IBM\n");
            std::abort();
          }
        }
        // Iterate, increase by the number of partners of this bond + 1 for bond
        // type
        j += iaparams.num + 1;
      }

      // Second round for triel
      if (softID > -1) {
        j = 0;
        while (j < p1.bl.n) {
          const int type_num = p1.bl.e[j];
          const Bonded_ia_parameters &iaparams = bonded_ia_params[type_num];
          const int type = iaparams.type;

          if (type == BONDED_IA_IBM_TRIEL) {
            // Our particle is the leading particle of a triel
            // Get second and third particle of the triangle
            Particle const *const p2 = local_particles[p1.bl.e[j + 1]];
            if (!p2) {
              runtimeErrorMsg()
                  << "{IBM_calc_volumes: 078 bond broken between particles "
                  << p1.p.identity << " and " << p1.bl.e[j + 1]
                  << " (particles not stored on the same node)} ";
              return;
            }
            Particle const *const p3 = local_particles[p1.bl.e[j + 2]];
            if (!p3) {
              runtimeErrorMsg()
                  << "{IBM_calc_volumes: 078 bond broken between particles "
                  << p1.p.identity << " and " << p1.bl.e[j + 2]
                  << " (particles not stored on the same node)} ";
              return;
            }

            // Unfold position of first node
            // this is to get a continuous trajectory with no jumps when box
            // boundaries are crossed
            auto const x1 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());
            auto const x2 = x1 + get_mi_vector(p2->r.p, x1, box_geo);
            auto const x3 = x1 + get_mi_vector(p3->r.p, x1, box_geo);

            // Volume of this tetrahedron
            // See Cha Zhang et.al. 2001, doi:10.1109/ICIP.2001.958278
            // http://research.microsoft.com/en-us/um/people/chazhang/publications/icip01_ChaZhang.pdf
            // The volume can be negative, but it is not necessarily the "signed
            // volume" in the above paper (the sign of the real "signed volume"
            // must be calculated using the normal vector; the result of the
            // calculation here is simply a term in the sum required to
            // calculate the volume of a particle). Again, see the paper. This
            // should be equivalent to the formulation using vector identities
            // in Krüger thesis

            const double v321 = x3[0] * x2[1] * x1[2];
            const double v231 = x2[0] * x3[1] * x1[2];
            const double v312 = x3[0] * x1[1] * x2[2];
            const double v132 = x1[0] * x3[1] * x2[2];
            const double v213 = x2[0] * x1[1] * x3[2];
            const double v123 = x1[0] * x2[1] * x3[2];

            tempVol[softID] +=
                1.0 / 6.0 * (-v321 + v231 + v312 - v132 - v213 + v123);
          }
          // Iterate, increase by the number of partners of this bond + 1 for
          // bond type
          j += iaparams.num + 1;
        }
      }
    }
  }

  for (int i = 0; i < MaxNumIBM; i++)
    VolumesCurrent[i] = 0;

  // Sum up and communicate
  MPI_Allreduce(&(tempVol.front()), &(VolumesCurrent.front()), MaxNumIBM,
                MPI_DOUBLE, MPI_SUM, comm_cart);
}

/** Calculate and add the volume force to each node */
void ImmersedBoundaries::calc_volume_force() {
  // Loop over all particles on local node
  for (int c = 0; c < local_cells.n; c++) {
    const Cell *const cell = local_cells.cell[c];

    for (int i = 0; i < cell->n; i++) {
      Particle &p1 = cell->part[i];

      // Check if particle has a BONDED_IA_IBM_TRIEL and a
      // BONDED_IA_IBM_VOLUME_CONSERVATION Basically this loops over all
      // triangles, not all particles First round to check for volume
      // conservation and virtual Loop over all bonds of this particle Actually
      // j loops over the bond-list, i.e. the bond partners (see
      // particle_data.hpp)
      int softID = -1;
      double volRef = 0.;
      double kappaV = 0.;
      int j = 0;
      while (j < p1.bl.n) {
        const int type_num = p1.bl.e[j];
        const Bonded_ia_parameters &iaparams = bonded_ia_params[type_num];
        const int type = iaparams.type;
        if (type == BONDED_IA_IBM_VOLUME_CONSERVATION) {
          if (!p1.p.is_virtual) {
            printf("Error. Encountered non-virtual particle with "
                   "VOLUME_CONSERVATION_IBM\n");
            std::abort();
          }
          softID = iaparams.p.ibmVolConsParameters.softID;
          volRef = iaparams.p.ibmVolConsParameters.volRef;
          kappaV = iaparams.p.ibmVolConsParameters.kappaV;
        }
        // Iterate, increase by the number of partners of this bond + 1 for bond
        // type
        j += iaparams.num + 1;
      }

      // Second round for triel
      if (softID > -1) {
        j = 0;
        while (j < p1.bl.n) {
          const int type_num = p1.bl.e[j];
          const Bonded_ia_parameters &iaparams = bonded_ia_params[type_num];
          const int type = iaparams.type;

          if (type == BONDED_IA_IBM_TRIEL) {
            // Our particle is the leading particle of a triel
            // Get second and third particle of the triangle
            Particle &p2 = *local_particles[p1.bl.e[j + 1]];
            Particle &p3 = *local_particles[p1.bl.e[j + 2]];

            // Unfold position of first node
            // this is to get a continuous trajectory with no jumps when box
            // boundaries are crossed
            auto const x1 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());

            // Unfolding seems to work only for the first particle of a triel
            // so get the others from relative vectors considering PBC
            auto const a12 = get_mi_vector(p2.r.p, x1, box_geo);
            auto const a13 = get_mi_vector(p3.r.p, x1, box_geo);

            // Now we have the true and good coordinates
            // Compute force according to eq. C.46 Krüger thesis
            // It is the same as deriving Achim's equation w.r.t x
            /*                        const double fact = kappaV * 1/6. *
            (IBMVolumesCurrent[softID] - volRef) / IBMVolumesCurrent[softID];

             double x2[3];
            double x3[3];

            for (int i=0; i < 3; i++)
            {
              x2[i] = x1[i] + a12[i];
              x3[i] = x1[i] + a13[i];
            }

             double n[3];
             vector_product(x3, x2, n);
             for (int k=0; k < 3; k++) p1.f.f[k] += fact*n[k];
             vector_product(x1, x3, n);
             for (int k=0; k < 3; k++) p2.f.f[k] += fact*n[k];
             vector_product(x2, x1, n);
             for (int k=0; k < 3; k++) p3.f.f[k] += fact*n[k];*/

            // This is Dupin 2008. I guess the result will be very similar as
            // the code above
            auto const n = vector_product(a12, a13);
            const double ln = n.norm();
            const double A = 0.5 * ln;
            const double fact = kappaV * (VolumesCurrent[softID] - volRef) /
                                VolumesCurrent[softID];

            auto const nHat = n / ln;
            auto const force = -fact * A * nHat;

            p1.f.f += force;
            p2.f.f += force;
            p3.f.f += force;
          }
          // Iterate, increase by the number of partners of this bond + 1 for
          // bond type
          j += iaparams.num + 1;
        }
      }
    }
  }
}

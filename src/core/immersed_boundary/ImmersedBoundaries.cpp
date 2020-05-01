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

#include "Particle.hpp"
#include "bonded_interactions/bonded_interaction_data.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"

#include <utils/constants.hpp>

#include <cstdio>

/** Volume conservation.
 *  Calculate volumes, volume force and add it to each virtual particle.
 *  This function is called from @ref integrate.
 */
void ImmersedBoundaries::volume_conservation(CellStructure &cs) {
  if (VolumeInitDone && !BoundariesFound) {
    return;
  }
  // Calculate volumes
  calc_volumes(cs);

  calc_volume_force(cs);
}

/** Initialize volume conservation */
void ImmersedBoundaries::init_volume_conservation(CellStructure &cs) {

  // Check since this function is called at the start of every integrate loop
  // Also check if volume has been set due to reading of a checkpoint
  if (!VolumeInitDone) {

    // Calculate volumes
    calc_volumes(cs);

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
  bonded_ia_params[bond_type].num = 0;

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

static const IBM_VolCons_Parameters *vol_cons_parameters(Particle const &p1) {
  auto it = boost::find_if(p1.bonds(), [](auto const &bond) {
    return bonded_ia_params[bond.bond_id()].type ==
           BONDED_IA_IBM_VOLUME_CONSERVATION;
  });

  return (it != p1.bonds().end())
             ? std::addressof(
                   bonded_ia_params[it->bond_id()].p.ibmVolConsParameters)
             : nullptr;
}

/** Calculate partial volumes on all compute nodes and call MPI to sum up.
 *  See @cite zhang01b, @cite dupin08a, @cite kruger12a.
 */
void ImmersedBoundaries::calc_volumes(CellStructure &cs) {

  // Partial volumes for each soft particle, to be summed up
  std::vector<double> tempVol(MaxNumIBM);

  // Loop over all particles on local node
  for (auto &p1 : cs.local_particles()) {
    auto vol_cons_params = vol_cons_parameters(p1);

    if (vol_cons_params) {
      cs.execute_bond_handler(p1, [softID = vol_cons_params->softID,
                                   &tempVol](Particle &p1, int bond_id,
                                             Utils::Span<Particle *> partners) {
        auto const &iaparams = bonded_ia_params[bond_id];

        if (iaparams.type == BONDED_IA_IBM_TRIEL) {
          // Our particle is the leading particle of a triel
          // Get second and third particle of the triangle
          Particle &p2 = *partners[0];
          Particle &p3 = *partners[1];

          // Unfold position of first node.
          // This is to get a continuous trajectory with no jumps when box
          // boundaries are crossed.
          auto const x1 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());
          auto const x2 = x1 + get_mi_vector(p2.r.p, x1, box_geo);
          auto const x3 = x1 + get_mi_vector(p3.r.p, x1, box_geo);

          // Volume of this tetrahedron
          // See @cite zhang01b
          // The volume can be negative, but it is not necessarily the
          // "signed volume" in the above paper (the sign of the real
          // "signed volume" must be calculated using the normal vector; the
          // result of the calculation here is simply a term in the sum
          // required to calculate the volume of a particle). Again, see the
          // paper. This should be equivalent to the formulation using
          // vector identities in @cite kruger12a

          const double v321 = x3[0] * x2[1] * x1[2];
          const double v231 = x2[0] * x3[1] * x1[2];
          const double v312 = x3[0] * x1[1] * x2[2];
          const double v132 = x1[0] * x3[1] * x2[2];
          const double v213 = x2[0] * x1[1] * x3[2];
          const double v123 = x1[0] * x2[1] * x3[2];

          tempVol[softID] +=
              1.0 / 6.0 * (-v321 + v231 + v312 - v132 - v213 + v123);
        }
        return false;
      });
    }
  }

  for (int i = 0; i < MaxNumIBM; i++)
    VolumesCurrent[i] = 0;

  // Sum up and communicate
  MPI_Allreduce(&(tempVol.front()), &(VolumesCurrent.front()), MaxNumIBM,
                MPI_DOUBLE, MPI_SUM, comm_cart);
}

/** Calculate and add the volume force to each node */
void ImmersedBoundaries::calc_volume_force(CellStructure &cs) {
  // Loop over all particles on local node
  for (auto &p1 : cs.local_particles()) {
    // Check if particle has a BONDED_IA_IBM_TRIEL and a
    // BONDED_IA_IBM_VOLUME_CONSERVATION. Basically this loops over all
    // triangles, not all particles. First round to check for volume
    // conservation.
    const IBM_VolCons_Parameters *ibmVolConsParameters =
        vol_cons_parameters(p1);
    if (not ibmVolConsParameters)
      return;

    auto current_volume = VolumesCurrent[ibmVolConsParameters->softID];

    // Second round for triel
    cs.execute_bond_handler(p1, [ibmVolConsParameters, current_volume](
                                    Particle &p1, int bond_id,
                                    Utils::Span<Particle *> partners) {
      auto const &iaparams = bonded_ia_params[bond_id];

      if (iaparams.type == BONDED_IA_IBM_TRIEL) {
        // Our particle is the leading particle of a triel
        // Get second and third particle of the triangle
        Particle &p2 = *partners[0];
        Particle &p3 = *partners[1];

        // Unfold position of first node.
        // This is to get a continuous trajectory with no jumps when box
        // boundaries are crossed.
        auto const x1 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());

        // Unfolding seems to work only for the first particle of a triel
        // so get the others from relative vectors considering PBC
        auto const a12 = get_mi_vector(p2.r.p, x1, box_geo);
        auto const a13 = get_mi_vector(p3.r.p, x1, box_geo);

        // Now we have the true and good coordinates
        // This is eq. (9) in @cite dupin08a.
        auto const n = vector_product(a12, a13);
        const double ln = n.norm();
        const double A = 0.5 * ln;
        const double fact = ibmVolConsParameters->kappaV *
                            (current_volume - ibmVolConsParameters->volRef) /
                            current_volume;

        auto const nHat = n / ln;
        auto const force = -fact * A * nHat;

        p1.f.f += force;
        p2.f.f += force;
        p3.f.f += force;
      }
      return false;
    });
  }
}

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

#include "BoxGeometry.hpp"
#include "Particle.hpp"
#include "communication.hpp"
#include "grid.hpp"
#include "ibm_volcons.hpp"

#include "bonded_interactions/bonded_interaction_data.hpp"

#include <utils/Span.hpp>
#include <utils/Vector.hpp>
#include <utils/constants.hpp>

#include <boost/mpi/collectives/all_reduce.hpp>

#include <functional>
#include <utility>
#include <vector>

/** Calculate volumes, volume force and add it to each virtual particle. */
void ImmersedBoundaries::volume_conservation(CellStructure &cs) {
  if (VolumeInitDone && !BoundariesFound) {
    return;
  }
  calc_volumes(cs);
  calc_volume_force(cs);
}

/** Initialize volume conservation */
void ImmersedBoundaries::init_volume_conservation(CellStructure &cs) {
  // Check since this function is called at the start of every integrate loop
  // Also check if volume has been set due to reading of a checkpoint
  if (not BoundariesFound) {
    BoundariesFound =
        boost::algorithm::any_of(bonded_ia_params, [](auto const &kv) {
          return (boost::get<IBMVolCons>(&(*kv.second)) != nullptr);
        });
  }

  if (!VolumeInitDone && BoundariesFound) {
    // Calculate volumes
    calc_volumes(cs);

    // Loop through all bonded interactions and check if we need to set the
    // reference volume
    for (auto &kv : bonded_ia_params) {
      if (auto *v = boost::get<IBMVolCons>(&(*kv.second))) {
        // This check is important because InitVolumeConservation may be called
        // accidentally during the integration. Then we must not reset the
        // reference
        BoundariesFound = true;
        if (v->volRef == 0.) {
          v->volRef = VolumesCurrent[v->softID];
        }
      }
    }

    VolumeInitDone = true;
  }
}

static const IBMVolCons *vol_cons_parameters(Particle const &p1) {
  auto it = boost::find_if(p1.bonds(), [](auto const &bond) {
    return boost::get<IBMVolCons>(bonded_ia_params.at(bond.bond_id()).get()) !=
           nullptr;
  });

  return (it != p1.bonds().end())
             ? boost::get<IBMVolCons>(bonded_ia_params.at(it->bond_id()).get())
             : nullptr;
}

/** Calculate partial volumes on all compute nodes and call MPI to sum up.
 *  See @cite zhang01b, @cite dupin08a, @cite kruger12a.
 */
void ImmersedBoundaries::calc_volumes(CellStructure &cs) {

  if (!BoundariesFound)
    return;

  // Partial volumes for each soft particle, to be summed up
  std::vector<double> tempVol(VolumesCurrent.size());

  // Loop over all particles on local node
  cs.bond_loop([&tempVol](Particle &p1, int bond_id,
                          Utils::Span<Particle *> partners) {
    auto const &iaparams = *bonded_ia_params.at(bond_id);
    auto vol_cons_params = vol_cons_parameters(p1);

    if (vol_cons_params &&
        boost::get<IBMTriel>(bonded_ia_params.at(bond_id).get()) != nullptr) {
      // Our particle is the leading particle of a triel
      // Get second and third particle of the triangle
      Particle &p2 = *partners[0];
      Particle &p3 = *partners[1];

      // Unfold position of first node.
      // This is to get a continuous trajectory with no jumps when box
      // boundaries are crossed.
      auto const x1 = unfolded_position(p1.r.p, p1.l.i, box_geo.length());
      auto const x2 = x1 + box_geo.get_mi_vector(p2.r.p, x1);
      auto const x3 = x1 + box_geo.get_mi_vector(p3.r.p, x1);

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

      tempVol[vol_cons_params->softID] +=
          1.0 / 6.0 * (-v321 + v231 + v312 - v132 - v213 + v123);
    }
    return false;
  });

  // Sum up and communicate
  boost::mpi::all_reduce(comm_cart, tempVol.data(),
                         static_cast<int>(tempVol.size()),
                         VolumesCurrent.data(), std::plus<double>());
}

/** Calculate and add the volume force to each node */
void ImmersedBoundaries::calc_volume_force(CellStructure &cs) {
  if (!BoundariesFound)
    return;

  cs.bond_loop([this](Particle &p1, int bond_id,
                      Utils::Span<Particle *> partners) {
    if (boost::get<IBMTriel>(bonded_ia_params.at(bond_id).get()) != nullptr) {
      // Check if particle has an IBM Triel bonded interaction and an
      // IBM VolCons bonded interaction. Basically this loops over all
      // triangles, not all particles. First round to check for volume
      // conservation.
      const IBMVolCons *ibmVolConsParameters = vol_cons_parameters(p1);
      if (not ibmVolConsParameters)
        return false;

      auto current_volume = VolumesCurrent[ibmVolConsParameters->softID];

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
      auto const a12 = box_geo.get_mi_vector(p2.r.p, x1);
      auto const a13 = box_geo.get_mi_vector(p3.r.p, x1);

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

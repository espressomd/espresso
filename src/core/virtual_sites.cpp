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

#include "config/config.hpp"

#ifdef VIRTUAL_SITES

#include "virtual_sites.hpp"

#include "Particle.hpp"
#include "communication.hpp"
#include "errorhandling.hpp"
#include "grid.hpp"
#include "integrate.hpp"
#include "nonbonded_interactions/nonbonded_interaction_data.hpp"

#include <utils/Vector.hpp>
#include <utils/math/quaternion.hpp>
#include <utils/quaternion.hpp>

#include <cmath>
#include <cstdio>
#include <tuple>

namespace {
std::shared_ptr<VirtualSites> m_virtual_sites;
}

const std::shared_ptr<VirtualSites> &virtual_sites() { return m_virtual_sites; }

void set_virtual_sites(std::shared_ptr<VirtualSites> const &v) {
  m_virtual_sites = v;
  recalc_forces = true;
}

#ifdef VIRTUAL_SITES_RELATIVE

/** Calculate the rotation quaternion and distance between two particles */
std::tuple<Utils::Quaternion<double>, double>
calculate_vs_relate_to_params(Particle const &p_vs,
                              Particle const &p_relate_to) {
  // get the distance between the particles
  auto d = ::box_geo.get_mi_vector(p_vs.pos(), p_relate_to.pos());

  // Check if the distance between virtual and non-virtual particles is larger
  // than minimum global cutoff. If so, warn user.
  auto const dist = d.norm();
  auto const min_global_cut = get_min_global_cut();
  if (dist > min_global_cut && n_nodes > 1 &&
      not virtual_sites()->get_override_cutoff_check()) {
    runtimeErrorMsg()
        << "Warning: The distance between virtual and non-virtual particle ("
        << dist << ") is larger than the minimum global cutoff ("
        << min_global_cut << "). This may lead to incorrect simulations under "
        << "certain conditions. Adjust the property system.min_global_cut to "
        << "increase the minimum cutoff.";
  }

  // If the distance between real & virtual particle is 0 we just set the
  // relative orientation to {1 0 0 0}, as it is irrelevant but needs to be
  // a valid quaternion
  if (dist == 0.) {
    return std::make_tuple(Utils::Quaternion<double>::identity(), dist);
  }

  // Now, calculate the quaternions which specify the angle between
  // the director of the particle we relate to and the vector
  // (particle_we_relate_to - this_particle)
  // The vs_relative implementation later obtains the director by multiplying
  // the quaternions representing the orientation of the real particle
  // with those in the virtual particle. The resulting quaternion is then
  // converted to a director.
  // We have quat_(real particle) * quat_(virtual particle)
  //   = quat_(obtained from desired director)
  // Resolving this for the quat_(virtual particle)

  d.normalize();

  // Obtain quaternion from desired director
  Utils::Quaternion<double> quat_director =
      Utils::convert_director_to_quaternion(d);

  // Define quaternion as described above
  auto relate_to_quat = p_relate_to.quat();
  auto quat =
      Utils::Quaternion<double>{{{{Utils::dot(relate_to_quat, quat_director),
                                   -quat_director[0] * relate_to_quat[1] +
                                       quat_director[1] * relate_to_quat[0] +
                                       quat_director[2] * relate_to_quat[3] -
                                       quat_director[3] * relate_to_quat[2],
                                   relate_to_quat[1] * quat_director[3] +
                                       relate_to_quat[0] * quat_director[2] -
                                       relate_to_quat[3] * quat_director[1] -
                                       relate_to_quat[2] * quat_director[0],
                                   quat_director[3] * relate_to_quat[0] -
                                       relate_to_quat[3] * quat_director[0] +
                                       relate_to_quat[2] * quat_director[1] -
                                       relate_to_quat[1] * quat_director[2]}}}};
  quat /= relate_to_quat.norm2();

  // Verify result
  Utils::Quaternion<double> qtemp = relate_to_quat * quat;
  for (unsigned int i = 0; i < 4; i++) {
    if (fabs(qtemp[i] - quat_director[i]) > 1E-9) {
      fprintf(stderr, "vs_relate_to: component %u: %f instead of %f\n", i,
              qtemp[i], quat_director[i]);
    }
  }
  return std::make_tuple(quat, dist);
}

#endif // VIRTUAL_SITES_RELATIVE
#endif // VIRTUAL_SITES

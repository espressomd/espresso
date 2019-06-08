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

#include "virtual_sites.hpp"
#include "communication.hpp"
#include "config.hpp"
#include "event.hpp"
#include "integrate.hpp"
#include "rotation.hpp"
#include "statistics.hpp"

#include <utils/constants.hpp>

#ifdef VIRTUAL_SITES

namespace {
std::shared_ptr<VirtualSites> m_virtual_sites;
}

const std::shared_ptr<VirtualSites> &virtual_sites() { return m_virtual_sites; }

void set_virtual_sites(std::shared_ptr<VirtualSites> const &v) {
  m_virtual_sites = v;
  recalc_forces = 1;
  invalidate_obs();
  on_ghost_flags_change();
}

#ifdef VIRTUAL_SITES_RELATIVE

void calculate_vs_relate_to_params(const Particle &p_current,
                                   const Particle &p_relate_to, double &l,
                                   Utils::Vector4d &quat) {
  // get the distance between the particles
  Utils::Vector3d d = get_mi_vector(p_current.r.p, p_relate_to.r.p);

  // Check, if the distance between virtual and non-virtual particles is larger
  // htan minimum global cutoff If so, warn user
  l = d.norm();
  if (l > min_global_cut && n_nodes > 1) {
    runtimeErrorMsg()
        << "Warning: The distance between virtual and non-virtual particle ("
        << l << ") is\nlarger than the minimum global cutoff ("
        << min_global_cut
        << "). This may lead to incorrect simulations\nunder certain "
           "conditions. Set the \"System()\" "
           "class property \"min_global_cut\" to\nincrease the minimum "
           "cutoff.\n";
  }

  // Now, calculate the quaternions which specify the angle between
  // the director of the particle we relate to and the vector
  // (paritlce_we_relate_to - this_particle)
  // The vs_relative implementation later obtains the director by multiplying
  // the quaternions representing the orientation of the real particle
  // with those in the virtual particle. The re quulting quaternion is then
  // converted to a director.
  // We have quat_(real particle) *quat_(virtual particle)
  // = quat_(obtained from desired director)
  // Resolving this for the quat_(virtual particle)

  // Normalize desired director
  int i;

  // If the distance between real & virtual particle is 0
  // we just set the relative orientation to 1 0 0 0, as it is irrelevant but
  // needs to be a valid quaternion
  if (l != 0) {
    for (i = 0; i < 3; i++)
      d[i] /= l;

    // Obtain quaternions from desired director
    Utils::Vector4d quat_director;
    convert_director_to_quat(d, quat_director);

    // Define quat as described above:
    double x = 0;
    for (i = 0; i < 4; i++)
      x += p_relate_to.r.quat[i] * p_relate_to.r.quat[i];

    quat[0] = 0;
    for (i = 0; i < 4; i++)
      quat[0] += p_relate_to.r.quat[i] * quat_director[i];

    quat[1] = -quat_director[0] * p_relate_to.r.quat[1] +
              quat_director[1] * p_relate_to.r.quat[0] +
              quat_director[2] * p_relate_to.r.quat[3] -
              quat_director[3] * p_relate_to.r.quat[2];
    quat[2] = p_relate_to.r.quat[1] * quat_director[3] +
              p_relate_to.r.quat[0] * quat_director[2] -
              p_relate_to.r.quat[3] * quat_director[1] -
              p_relate_to.r.quat[2] * quat_director[0];
    quat[3] = quat_director[3] * p_relate_to.r.quat[0] -
              p_relate_to.r.quat[3] * quat_director[0] +
              p_relate_to.r.quat[2] * quat_director[1] -
              p_relate_to.r.quat[1] * quat_director[2];
    for (i = 0; i < 4; i++)
      quat[i] /= x;

    // Verify result
    double qtemp[4];
    multiply_quaternions(p_relate_to.r.quat, quat, qtemp);
    for (i = 0; i < 4; i++)
      if (fabs(qtemp[i] - quat_director[i]) > 1E-9)
        fprintf(stderr, "vs_relate_to: component %d: %f instead of %f\n", i,
                qtemp[i], quat_director[i]);
  } else {
    quat[0] = 1;
    quat[1] = quat[2] = quat[3] = 0;
  }
}

// Setup the virtual_sites_relative properties of a particle so that the given
// virtual particle will follow the given real particle Local version, expects
// both particles to be accessible through local_particles and only executes the
// changes on the virtual site locally
int local_vs_relate_to(Particle *p_current, const Particle *p_relate_to) {
  Utils::Vector4d quat;

  double l;
  calculate_vs_relate_to_params(*p_current, *p_relate_to, l, quat);

  // Set the particle id of the particle we want to relate to, the distance
  // and the relative orientation
  p_current->p.vs_relative.to_particle_id = p_relate_to->identity();
  p_current->p.vs_relative.distance = l;
  for (int i = 0; i < 4; i++)
    p_current->p.vs_relative.rel_orientation[i] = quat[i];
  return ES_OK;
}

// Setup the virtual_sites_relative properties of a particle so that the given
// virtual particle will follow the given real particle
int vs_relate_to(int part_num, int relate_to) {
  // Get the data for the particle we act on and the one we want to relate
  // it to.
  auto const &p_current = get_particle_data(part_num);
  auto const &p_relate_to = get_particle_data(relate_to);

  Utils::Vector4d quat;
  double l;
  calculate_vs_relate_to_params(p_current, p_relate_to, l, quat);

  // Set the particle id of the particle we want to relate to, the distance
  // and the relative orientation
  set_particle_vs_relative(part_num, relate_to, l, quat.data());
  set_particle_virtual(part_num, 1);

  return ES_OK;
}

#endif
#endif

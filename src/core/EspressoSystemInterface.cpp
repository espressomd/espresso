/*
  Copyright (C) 2014,2015,2016 The ESPResSo project

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
#include "EspressoSystemInterface.hpp"
#include "cells.hpp"
#include "cuda_interface.hpp"
#include "grid.hpp"
#include "particle_data.hpp"

#include <iostream>

/* Initialize instance pointer */
EspressoSystemInterface *EspressoSystemInterface::m_instance = 0;

/********************************************************************************************/

void EspressoSystemInterface::gatherParticles() {
// get particles from other nodes
#ifdef CUDA
  if (m_gpu) {
    if (gpu_get_global_particle_vars_pointer_host()->communication_enabled) {
      ESIF_TRACE(puts("Calling copy_part_data_to_gpu()"));
      copy_part_data_to_gpu(local_cells.particles());
      reallocDeviceMemory(
          gpu_get_global_particle_vars_pointer_host()->number_of_particles);
      if (m_splitParticleStructGpu && (this_node == 0))
        split_particle_struct();
    }
  }
#endif

  if (needsQ() || needsR() || needsDip() || needsQuatu()) {
    R.clear();

#ifdef ELECTROSTATICS
    Q.clear();
#endif
#ifdef DIPOLES
    Dip.clear();
#endif

#ifdef ROTATION
    Quatu.clear();
#endif
    for (auto const &p : local_cells.particles()) {
      if (needsR())
        R.push_back(Vector3(p.r.p));

#ifdef ELECTROSTATICS
      if (needsQ())
        Q.push_back(p.p.q);
#endif
#ifdef DIPOLES
      if (needsDip())
        Dip.emplace_back(Vector3{p.r.dip[0], p.r.dip[1], p.r.dip[2]});
#endif
#ifdef ROTATION
      if (needsQuatu())
        Quatu.emplace_back(Vector3{p.r.quatu[0], p.r.quatu[1], p.r.quatu[2]});
#endif
    }
  }
}

void EspressoSystemInterface::init() { gatherParticles(); }

void EspressoSystemInterface::update() { gatherParticles(); }

SystemInterface::const_vec_iterator &EspressoSystemInterface::rBegin() {
  m_r_begin = R.begin();
  return m_r_begin;
}

const SystemInterface::const_vec_iterator &EspressoSystemInterface::rEnd() {
  m_r_end = R.end();
  return m_r_end;
}

#ifdef DIPOLES
SystemInterface::const_vec_iterator &EspressoSystemInterface::dipBegin() {
  m_dip_begin = Dip.begin();
  return m_dip_begin;
}

const SystemInterface::const_vec_iterator &EspressoSystemInterface::dipEnd() {
  m_dip_end = Dip.end();
  return m_dip_end;
}
#endif

#ifdef ELECTROSTATICS
SystemInterface::const_real_iterator &EspressoSystemInterface::qBegin() {
  m_q_begin = Q.begin();
  return m_q_begin;
}

const SystemInterface::const_real_iterator &EspressoSystemInterface::qEnd() {
  m_q_end = Q.end();
  return m_q_end;
}
#endif

#ifdef ROTATION
SystemInterface::const_vec_iterator &EspressoSystemInterface::quatuBegin() {
  m_quatu_begin = Quatu.begin();
  return m_quatu_begin;
}

const SystemInterface::const_vec_iterator &EspressoSystemInterface::quatuEnd() {
  m_quatu_end = Quatu.end();
  return m_quatu_end;
}
#endif

unsigned int EspressoSystemInterface::npart() { return m_npart; }

SystemInterface::Vector3 EspressoSystemInterface::box() {
  return Vector3d{box_l[0], box_l[1], box_l[2]};
}

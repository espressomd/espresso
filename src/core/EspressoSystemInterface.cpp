/*
  Copyright (C) 2014 The ESPResSo project
  
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
#include "particle_data.hpp"
#include "grid.hpp"
#include "cuda_interface.hpp"

#include <iostream>

/* Initialize instance pointer */
EspressoSystemInterface *EspressoSystemInterface::m_instance = 0;

/* Need explicite specialization, otherwise some compilers do not produce the objects. */

template class EspressoSystemInterface::const_iterator<SystemInterface::Real>;
template class EspressoSystemInterface::const_iterator<SystemInterface::Vector3>;
template class EspressoSystemInterface::const_iterator<int>;

/********************************************************************************************/

template<class value_type>
value_type EspressoSystemInterface::const_iterator<value_type>::operator*() const {
  return (*m_const_iterator);
}

template<class value_type>
SystemInterface::const_iterator<value_type> &EspressoSystemInterface::const_iterator<value_type>::operator=(const SystemInterface::const_iterator<value_type> &rhs) {
   m_const_iterator = static_cast<const EspressoSystemInterface::const_iterator<value_type> &>(rhs).m_const_iterator;
  return *this;
}

template<class value_type>
EspressoSystemInterface::const_iterator<value_type> &EspressoSystemInterface::const_iterator<value_type>::operator=(typename std::vector<value_type>::const_iterator rhs) {
   m_const_iterator = rhs;
  return *this;
}

template<class value_type>
bool EspressoSystemInterface::const_iterator<value_type>::operator==(SystemInterface::const_iterator<value_type> const &rhs) const {
   return (m_const_iterator == static_cast<const EspressoSystemInterface::const_iterator<value_type> &>(rhs).m_const_iterator);
}

template<class value_type>
bool EspressoSystemInterface::const_iterator<value_type>::operator!=(SystemInterface::const_iterator<value_type> const &rhs) const {
   return (m_const_iterator != static_cast<const EspressoSystemInterface::const_iterator<value_type> &>(rhs).m_const_iterator);
}

template<class value_type>
SystemInterface::const_iterator<value_type> &EspressoSystemInterface::const_iterator<value_type>::operator++() {
  ++m_const_iterator;
  return *this;
}

/********************************************************************************************/
 
void EspressoSystemInterface::gatherParticles() {
  Cell *cell;
  Particle *p;
  int i,c,np;

  // get particles from other nodes
#ifdef CUDA
  if (m_gpu)
  {
    if(gpu_get_global_particle_vars_pointer_host()->communication_enabled) {
      ESIF_TRACE(puts("Calling copy_part_data_to_gpu()"));
      copy_part_data_to_gpu();
      reallocDeviceMemory(gpu_get_global_particle_vars_pointer_host()->number_of_particles);
      if(m_splitParticleStructGpu && (this_node == 0)) 
	split_particle_struct();
    }
  }
#endif

  if (needsQ() || needsR() ||needsDip()|| needsQuatu()) {
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

    for (c = 0; c < local_cells.n; c++) {
      cell = local_cells.cell[c];
      p  = cell->part;
      np = cell->n;

      if(needsR())
	R.reserve(R.size()+np);

#ifdef ELECTROSTATICS
      if(needsQ())
	Q.reserve(Q.size()+np);
#endif
#ifdef DIPOLES
      if(needsDip())
	Dip.reserve(Dip.size()+np);
#endif

#ifdef ROTATION
      if(needsQuatu())
	Quatu.reserve(Quatu.size()+np);
#endif

      for(i = 0; i < np; i++) {

	if(needsR())
	  R.push_back(Vector3(p[i].r.p));

#ifdef ELECTROSTATICS
	if(needsQ())
	  Q.push_back(p[i].p.q);
#endif
#ifdef DIPOLES
	if(needsDip())
	  Dip.push_back(Vector3(p[i].r.dip));
#endif
#ifdef ROTATION
	if(needsQuatu())
	  Quatu.push_back(Vector3(p[i].r.quatu));
#endif
      }
    }
  }
}

void EspressoSystemInterface::init() {
  gatherParticles();
}

void EspressoSystemInterface::update() {
  gatherParticles();
}

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

unsigned int EspressoSystemInterface::npart() {
  return m_npart;
}

SystemInterface::Vector3 EspressoSystemInterface::box() {
  return Vector3(box_l);
}


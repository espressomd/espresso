#include "EspressoSystemInterface.hpp"
#include "cells.hpp"
#include "particle_data.hpp"
#include "grid.hpp"
#include "cuda_interface.hpp"

#include <iostream>

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

  if (needsQ() || needsR()) {
    R.clear();
    #ifdef ELECTROSTATICS
    Q.clear();
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
      for(i = 0; i < np; i++) {
	if(needsR())
	  R.push_back(Vector3(p[i].r.p));
#ifdef ELECTROSTATICS
	if(needsQ())
	  Q.push_back(p[i].p.q);
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

unsigned int EspressoSystemInterface::npart() {
  return m_npart;
}

SystemInterface::Vector3 EspressoSystemInterface::box() {
  return Vector3(box_l);
}



EspressoSystemInterface espressoSystemInterface;

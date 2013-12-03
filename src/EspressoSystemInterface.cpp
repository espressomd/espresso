#include "EspressoSystemInterface.hpp"
#include "cells.hpp"
#include "particle_data.hpp"
#include "grid.hpp"
#include "cuda_interface.hpp"

#include <iostream>

/********************************************************************************************/

template<class value_type>
const value_type EspressoSystemInterface::const_iterator<value_type>::operator*() const {
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
  R.clear();
  Q.clear();

  // get particles from other nodes
  CUDA_global_part_vars* global_part_vars_host = gpu_get_global_particle_vars_pointer_host();
  if ( global_part_vars_host->communication_enabled == 1)
  {
//     np = global_part_vars_host->number_of_particles;
//     cuda_mpi_get_particles(particle_data_host);
//     if (this_node > 0)
//       return;
//     for (int i = 0; i < np; i++)
//     {
//       R.push_back(Vector3(particle_data_host[i].p));
// #ifdef ELECTROSTATICS
//       Q.push_back(particle_data_host[i].q);
// #endif
    // }
  }
  else // only get particles from local node
    {
      for (c = 0; c < local_cells.n; c++) {
	cell = local_cells.cell[c];
	p  = cell->part;
	np = cell->n;
	R.reserve(R.size()+np);
#ifdef ELECTROSTATICS
	Q.reserve(Q.size()+np);
#endif
	for(i = 0; i < np; i++) {
	  R.push_back(Vector3(p[i].r.p));
#ifdef ELECTROSTATICS
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

SystemInterface::const_real_iterator &EspressoSystemInterface::qBegin() {
  m_q_begin = Q.begin();
  return m_q_begin;
}

const SystemInterface::const_real_iterator &EspressoSystemInterface::qEnd() {
  m_q_end = Q.end();
  return m_q_end;
}

unsigned int EspressoSystemInterface::npart() {
  return R.size();
}

SystemInterface::Vector3 EspressoSystemInterface::box() {
  return Vector3(box_l);
}

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
#ifndef ESPRESSOSYSTEMINTERFACE_H
#define ESPRESSOSYSTEMINTERFACE_H

#include <stdio.h>

#include "SystemInterface.hpp"
#include "cuda_interface.hpp"

// This debug header has to be the last thing to include, because it
// #defines malloc to be something else (!!!) which will lead to
// ultimately obscure errors in the standard library.  This is
// literally the worst hack I have ever seen!
#include "debug.hpp"

/** Syntactic sugar */
#define espressoSystemInterface EspressoSystemInterface::Instance()

class EspressoSystemInterface : public SystemInterface {
public:
  static EspressoSystemInterface &Instance() {
    if(!m_instance)
      m_instance = new EspressoSystemInterface;

    return *m_instance;
  };

  static EspressoSystemInterface * _Instance() {
      if(!m_instance)
          m_instance = new EspressoSystemInterface;

      return m_instance;
  };

  void init();
  void update();

  SystemInterface::Vector3 box();
  unsigned int npart();

  typedef std::vector<Vector3> Vector3Container;
  typedef std::vector<Real> RealContainer;

  template<class value_type>
  class const_iterator : public SystemInterface::const_iterator<value_type> {
  public:
    value_type operator*() const;
    SystemInterface::const_iterator<value_type> &operator=(const SystemInterface::const_iterator<value_type> &rhs);
    EspressoSystemInterface::const_iterator<value_type> &operator=(typename std::vector<value_type>::const_iterator rhs);
    bool operator==(SystemInterface::const_iterator<value_type> const &rhs) const;
    bool operator!=(SystemInterface::const_iterator<value_type> const &rhs) const;
    SystemInterface::const_iterator<value_type> &operator++();
  private:    
   typename std::vector<value_type>::const_iterator m_const_iterator;
  };

  typedef const_iterator<Real> const_real_iterator;
  typedef const_iterator<Vector3> const_vec_iterator;
  typedef const_iterator<int> const_int_iterator;


  // Particle position
  SystemInterface::const_vec_iterator &rBegin();
  const SystemInterface::const_vec_iterator &rEnd();
  bool hasR() { return true; };
  
  // Dipole moment
#ifdef DIPOLES  
  SystemInterface::const_vec_iterator &dipBegin();
  const SystemInterface::const_vec_iterator &dipEnd();
  bool hasDip() { return true; };
#endif

#ifdef ELECTROSTATICS
  SystemInterface::const_real_iterator &qBegin();
  const SystemInterface::const_real_iterator &qEnd();
  bool hasQ() { return true; };
#endif

#ifdef ROTATION
  SystemInterface::const_vec_iterator &quatuBegin();
  const SystemInterface::const_vec_iterator &quatuEnd();
  bool hasQuatu() { return true; };
#endif

#ifdef CUDA
  float *rGpuBegin() { return m_r_gpu_begin; };
  float *rGpuEnd() { return m_r_gpu_end; };
  bool hasRGpu() { return true; };
  bool requestRGpu() { 
    m_needsRGpu = hasRGpu();
    m_splitParticleStructGpu |= m_needsRGpu;
    m_gpu |= m_needsRGpu;
    if(m_gpu)
      enableParticleCommunication();
    return m_needsRGpu; 
  };
#ifdef DIPOLES
  float *dipGpuBegin() { return m_dip_gpu_begin; };
  float *dipGpuEnd() { return m_dip_gpu_end; };
  bool hasDipGpu() { return true; };
  bool requestDipGpu() { 
    m_needsDipGpu = hasDipGpu();
    m_splitParticleStructGpu |= m_needsRGpu;
    m_gpu |= m_needsRGpu;
    if(m_gpu)
      enableParticleCommunication();
    return m_needsDipGpu; 
  };
#endif
  float *vGpuBegin() { return m_v_gpu_begin; };
  float *vGpuEnd() { return m_v_gpu_end; };
  bool hasVGpu() { return true; };
  bool requestVGpu() { 
    m_needsVGpu = hasVGpu();
    m_splitParticleStructGpu |= m_needsVGpu;
    m_gpu |= m_needsVGpu;
    if(m_gpu)
      enableParticleCommunication();
    return m_needsVGpu; 
  };

  float *qGpuBegin() { return m_q_gpu_begin; };
  float *qGpuEnd() { return m_q_gpu_end; };
  bool hasQGpu() { return true; };
  bool requestQGpu() { 
    m_needsQGpu = hasQGpu(); 
    m_splitParticleStructGpu |= m_needsQGpu;
    m_gpu |= m_needsQGpu;
    if(m_gpu)
      enableParticleCommunication();
    return m_needsQGpu; 
  };

  float *quatuGpuBegin() { return m_quatu_gpu_begin; };
  float *quatuGpuEnd() { return m_quatu_gpu_end; };
  bool hasQuatuGpu() { return true; };
  bool requestQuatuGpu() { 
    m_needsQuatuGpu = hasQuatuGpu(); 
    m_splitParticleStructGpu |= m_needsQuatuGpu;
    m_gpu |= m_needsQuatuGpu;
    if(m_gpu)
      enableParticleCommunication();
    return m_needsQuatuGpu; 
  };

  bool requestParticleStructGpu() {
    m_needsParticleStructGpu = true;
    m_gpu |= m_needsParticleStructGpu;
    if(m_gpu)
      enableParticleCommunication();
    return true;
  };

  float *fGpuBegin() { return gpu_get_particle_force_pointer(); };
  float *fGpuEnd() { return gpu_get_particle_force_pointer() + 3*m_gpu_npart; };
  float *eGpu() { return (float *)gpu_get_energy_pointer(); };
  float *torqueGpuBegin() { return (float *)gpu_get_particle_torque_pointer(); };
  float *torqueGpuEnd() { return (float *)(gpu_get_particle_torque_pointer()) + 3*m_gpu_npart; };
  bool hasFGpu() { return true; };
  bool requestFGpu() {
    m_needsFGpu = hasFGpu();
    m_gpu |= m_needsFGpu;
    if(m_gpu)
      enableParticleCommunication();
    return m_needsFGpu;
  };

#ifdef ROTATION
  bool hasTorqueGpu() { return true; };
  bool requestTorqueGpu() {
    m_needsTorqueGpu = hasTorqueGpu();
    m_gpu |= m_needsTorqueGpu;
    if(m_gpu)
      enableParticleCommunication();
    return m_needsTorqueGpu;
  };
#endif

#endif

  unsigned int npart_gpu() {
#ifdef CUDA
    return m_gpu_npart;
#else
    return 0;
#endif
  };

protected:
  static EspressoSystemInterface *m_instance;
  EspressoSystemInterface() : m_gpu_npart(0), m_gpu(false), m_r_gpu_begin(0), m_r_gpu_end(0), m_dip_gpu_begin(0), m_v_gpu_begin(0), m_v_gpu_end(0), m_q_gpu_begin(0),  m_q_gpu_end(0), m_quatu_gpu_begin(0),  m_quatu_gpu_end(0), m_needsParticleStructGpu(false), m_splitParticleStructGpu(false)  {};
  virtual ~EspressoSystemInterface() {}

  void gatherParticles();
  void split_particle_struct();
#ifdef CUDA
  void enableParticleCommunication() {
    if(!gpu_get_global_particle_vars_pointer_host()->communication_enabled) {
      ESIF_TRACE(puts("gpu communication not enabled;"));
      ESIF_TRACE(puts("enableParticleCommunication"));
      gpu_init_particle_comm();
      cuda_bcast_global_part_params();
      reallocDeviceMemory(gpu_get_global_particle_vars_pointer_host()->number_of_particles);
    }
  };
  void reallocDeviceMemory(int n);
#endif

  Vector3Container R;

  #ifdef ELECTROSTATICS
  RealContainer Q;
  #endif

#ifdef DIPOLES
  Vector3Container Dip;
#endif
  #ifdef ROTATION
  Vector3Container Quatu;
  #endif

  const_vec_iterator m_r_begin;
  const_vec_iterator m_r_end;

#ifdef DIPOLES
  const_vec_iterator m_dip_begin;
  const_vec_iterator m_dip_end;
#endif

  const_real_iterator m_q_begin;
  const_real_iterator m_q_end;

  const_vec_iterator m_quatu_begin;
  const_vec_iterator m_quatu_end;

  int m_gpu_npart;
  bool m_gpu;

  float *m_r_gpu_begin;
  float *m_r_gpu_end;

  float *m_dip_gpu_begin;
  float *m_dip_gpu_end;



  float *m_v_gpu_begin;
  float *m_v_gpu_end;

  float *m_q_gpu_begin;
  float *m_q_gpu_end;

  float *m_quatu_gpu_begin;
  float *m_quatu_gpu_end;

  unsigned int m_npart;
  Vector3 m_box;

  bool m_needsParticleStructGpu;
  bool m_splitParticleStructGpu;
};

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

#endif

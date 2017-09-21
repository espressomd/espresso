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

/** Syntactic sugar */
#define espressoSystemInterface EspressoSystemInterface::Instance()

class EspressoSystemInterface : public SystemInterface {
public:
  static EspressoSystemInterface &Instance() {
    if (!m_instance)
      m_instance = new EspressoSystemInterface;

    return *m_instance;
  };

  static EspressoSystemInterface *_Instance() {
    if (!m_instance)
      m_instance = new EspressoSystemInterface;

    return m_instance;
  };

  void init();
  void update();

#ifdef CUDA
  float *rGpuBegin() { return m_r_gpu_begin; };
  float *rGpuEnd() { return m_r_gpu_end; };
  bool hasRGpu() { return true; };
  bool requestRGpu() {
    m_needsRGpu = hasRGpu();
    m_splitParticleStructGpu |= m_needsRGpu;
    m_gpu |= m_needsRGpu;
    if (m_gpu)
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
    if (m_gpu)
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
    if (m_gpu)
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
    if (m_gpu)
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
    if (m_gpu)
      enableParticleCommunication();
    return m_needsQuatuGpu;
  };

  bool requestParticleStructGpu() {
    m_needsParticleStructGpu = true;
    m_gpu |= m_needsParticleStructGpu;
    if (m_gpu)
      enableParticleCommunication();
    return true;
  };

  float *fGpuBegin() { return gpu_get_particle_force_pointer(); };
  float *fGpuEnd() {
    return gpu_get_particle_force_pointer() + 3 * m_gpu_npart;
  };
  float *eGpu() { return (float *)gpu_get_energy_pointer(); };
  float *torqueGpuBegin() {
    return (float *)gpu_get_particle_torque_pointer();
  };
  float *torqueGpuEnd() {
    return (float *)(gpu_get_particle_torque_pointer()) + 3 * m_gpu_npart;
  };
  bool hasFGpu() { return true; };
  bool requestFGpu() {
    m_needsFGpu = hasFGpu();
    m_gpu |= m_needsFGpu;
    if (m_gpu)
      enableParticleCommunication();
    return m_needsFGpu;
  };

#ifdef ROTATION
  bool hasTorqueGpu() { return true; };
  bool requestTorqueGpu() {
    m_needsTorqueGpu = hasTorqueGpu();
    m_gpu |= m_needsTorqueGpu;
    if (m_gpu)
      enableParticleCommunication();
    return m_needsTorqueGpu;
  };
#endif

#endif

  Vector3d box() const override;

  unsigned int npart_gpu() {
#ifdef CUDA
    return m_gpu_npart;
#else
    return 0;
#endif
  };

protected:
  static EspressoSystemInterface *m_instance;
  EspressoSystemInterface()
      : m_gpu_npart(0), m_gpu(false), m_r_gpu_begin(0), m_r_gpu_end(0),
        m_dip_gpu_begin(0), m_v_gpu_begin(0), m_v_gpu_end(0), m_q_gpu_begin(0),
        m_q_gpu_end(0), m_quatu_gpu_begin(0), m_quatu_gpu_end(0),
        m_needsParticleStructGpu(false), m_splitParticleStructGpu(false){};
  virtual ~EspressoSystemInterface() {}

  void gatherParticles();
  void split_particle_struct();
#ifdef CUDA
  void enableParticleCommunication() {
    if (!gpu_get_global_particle_vars_pointer_host()->communication_enabled) {
      ESIF_TRACE(puts("gpu communication not enabled;"));
      ESIF_TRACE(puts("enableParticleCommunication"));
      gpu_init_particle_comm();
      cuda_bcast_global_part_params();
      reallocDeviceMemory(
          gpu_get_global_particle_vars_pointer_host()->number_of_particles);
    }
  };
  void reallocDeviceMemory(int n);
#endif

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

  bool m_needsParticleStructGpu;
  bool m_splitParticleStructGpu;
};

/********************************************************************************************/

#endif

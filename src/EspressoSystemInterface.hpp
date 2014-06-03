#ifndef ESPRESSOSYSTEMINTERFACE_H
#define ESPRESSOSYSTEMINTERFACE_H

#define ESIF_TRACE(A)

#include "SystemInterface.hpp"
#include "cuda_interface.hpp"

class EspressoSystemInterface : public SystemInterface {
public:
  EspressoSystemInterface() : m_gpu_npart(0), m_gpu(false), m_r_gpu_begin(0), m_r_gpu_end(0), m_v_gpu_begin(0), m_v_gpu_end(0), m_q_gpu_begin(0),  m_q_gpu_end(0), m_needsParticleStructGpu(false), m_splitParticleStructGpu(false)  {};
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


  SystemInterface::const_vec_iterator &rBegin();
  const SystemInterface::const_vec_iterator &rEnd();
  bool hasR() { return true; };

#ifdef ELECTROSTATICS
  SystemInterface::const_real_iterator &qBegin();
  const SystemInterface::const_real_iterator &qEnd();
  bool hasQ() { return true; };
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

  bool requestParticleStructGpu() {
    m_needsParticleStructGpu = true;
    m_gpu |= m_needsParticleStructGpu;
    if(m_gpu)
      enableParticleCommunication();
    return true;
  }

  float *fGpuBegin() { return (float *)gpu_get_particle_force_pointer(); };
  float *fGpuEnd() { return (float *)(gpu_get_particle_force_pointer()) + 3*m_gpu_npart; };
  bool hasFGpu() { return true; };
  bool requestFGpu() {
    m_needsFGpu = hasFGpu();
    m_gpu |= m_needsFGpu;
    if(m_gpu)
      enableParticleCommunication();
    return m_needsFGpu;
  };
#endif

  unsigned int npart_gpu() {
#ifdef CUDA
    return m_gpu_npart;
#else
    return 0;
#endif
  };


protected:
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

  const_vec_iterator m_r_begin;
  const_vec_iterator m_r_end;

  const_real_iterator m_q_begin;
  const_real_iterator m_q_end;

  int m_gpu_npart;
  bool m_gpu;

  float *m_r_gpu_begin;
  float *m_r_gpu_end;

  float *m_v_gpu_begin;
  float *m_v_gpu_end;

  float *m_q_gpu_begin;
  float *m_q_gpu_end;

  unsigned int m_npart;
  Vector3 m_box;

  bool m_needsParticleStructGpu;
  bool m_splitParticleStructGpu;
};

extern EspressoSystemInterface espressoSystemInterface;

#endif

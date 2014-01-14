#ifndef ESPRESSOSYSTEMINTERFACE_H
#define ESPRESSOSYSTEMINTERFACE_H

#include "SystemInterface.hpp"
#include "particle_data.hpp"
#include "cuda_interface.hpp"

class EspressoSystemInterface : public SystemInterface {
public:
  EspressoSystemInterface() : m_gpu_npart(0), m_r_gpu_begin(0), m_r_gpu_end(0), m_v_gpu_begin(0), m_v_gpu_end(0), m_q_gpu_begin(0),  m_q_gpu_end(0), m_gpu(false), m_needsParticleStructGpu(false), m_splitParticleStructGpu(false)  {};
  void init();
  void update();

  SystemInterface::Vector3 box();
  unsigned int npart();

  typedef std::vector<Vector3> Vector3Container;
  typedef std::vector<Real> RealContainer;

  template<class value_type>
  class const_iterator : public SystemInterface::const_iterator<value_type> {
  public:
    const value_type operator*() const;
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

  float *rGpuBegin() { return m_r_gpu_begin; };
  float *rGpuEnd() { return m_r_gpu_end; };
  bool hasRGpu() { return true; };
  bool requestRGpu() { 
    m_needsRGpu = hasRGpu();
    m_splitParticleStructGpu |= m_needsRGpu;
    m_gpu |= m_needsRGpu;
    return m_needsRGpu; 
  };

  float *vGpuBegin() { return m_v_gpu_begin; };
  float *vGpuEnd() { return m_v_gpu_end; };
  bool hasVGpu() { return true; };
  bool requestVGpu() { 
    m_needsVGpu = hasVGpu();
    m_splitParticleStructGpu |= m_needsVGpu;
    m_gpu |= m_needsVGpu;
    return m_needsVGpu; 
  };

  float *qGpuBegin() { return m_q_gpu_begin; };
  float *qGpuEnd() { return m_q_gpu_end; };
  bool hasQGpu() { return true; };
  bool requestQGpu() { 
    m_needsQGpu = hasQGpu(); 
    m_splitParticleStructGpu |= m_needsQGpu;
    m_gpu |= m_needsQGpu;
    return m_needsQGpu; 
  };

  bool requestParticleStructGpu() {
    m_needsParticleStructGpu = true;
    m_gpu |= m_needsParticleStructGpu;
    return true;
  }

protected:
  void gatherParticles();
  void split_particle_struct();

  Vector3Container R;
  RealContainer Q;

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

/* Need explicite specialization, otherwise some compilers do not produce the objects. */

template class EspressoSystemInterface::const_iterator<SystemInterface::Real>;
template class EspressoSystemInterface::const_iterator<SystemInterface::Vector3>;
template class EspressoSystemInterface::const_iterator<int>;

extern EspressoSystemInterface espressoSystemInterface;

#endif

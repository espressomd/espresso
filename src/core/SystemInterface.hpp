#ifndef SYSTEMINTERFACE_H
#define SYSTEMINTERFACE_H

#include "Vector.hpp"
#include <vector>

/** @TODO: Turn needsXY in getter/setter **/

class SystemInterface {
public:
  SystemInterface() : m_needsR(false), m_needsV(false), m_needsQ(false), m_needsRGpu(false), m_needsVGpu(false), m_needsQGpu(false) {};
  typedef Vector3d Vector3;
  typedef double Real;

  virtual void init() = 0;
  virtual void update() = 0;

  template<class value_type>
  class const_iterator {
  public:
    virtual value_type operator*() const = 0;
    virtual const_iterator<value_type> &operator=(const const_iterator<value_type> &rhs) = 0;
    virtual bool operator==(const_iterator<value_type> const &rhs) const = 0;
    virtual bool operator!=(const_iterator<value_type> const &rhs) const = 0;
    virtual const_iterator<value_type> &operator++() = 0;
  }; 

  class null_vec_iterator : public SystemInterface::const_iterator<Vector3> {
  public:
    Vector3 operator*() const { return Vector3(); };
    const_iterator<Vector3> &operator=(const const_iterator<Vector3> &rhs) { return *this; };
    bool operator==(const_iterator<Vector3> const &rhs) const { return true; };
    bool operator!=(const_iterator<Vector3> const &rhs) const { return false; };
    const_iterator<Vector3> &operator++() {return *this;};
  };

  null_vec_iterator null_vector;

  class null_real_iterator : public SystemInterface::const_iterator<Real> {
  public:
    Real operator*() const { return Real(); };
    const_iterator<Real> &operator=(const const_iterator<Real> &rhs) { return *this; };
    bool operator==(const_iterator<Real> const &rhs) const { return true; };
    bool operator!=(const_iterator<Real> const &rhs) const { return false; };
    const_iterator<Real> &operator++() {return *this;};
  };

  null_real_iterator null_scalar;

  template<class stl_container>
  class const_iterator_stl : public SystemInterface::const_iterator<typename stl_container::value_type> {
  public:
    const_iterator_stl(typename stl_container::const_iterator it) {
      m_const_iterator = it;
    };
    const_iterator_stl() {};
    const typename stl_container::value_type operator*() const {
      return (*m_const_iterator);
    };
    SystemInterface::const_iterator<typename stl_container::value_type> &operator=(const SystemInterface::const_iterator<typename stl_container::value_type> &rhs) {
      m_const_iterator = static_cast<const const_iterator_stl<stl_container> &>(rhs).m_const_iterator;
      return static_cast<SystemInterface::const_iterator<typename stl_container::value_type> &>(*this);
    };
    const_iterator_stl<stl_container> &operator=(const typename stl_container::const_iterator &rhs) {
      m_const_iterator = rhs;
      return *this;
    };

    const_iterator_stl<stl_container> &operator=(const const_iterator_stl<stl_container> &rhs) {
      m_const_iterator = rhs.m_const_iterator;
      return *this;
    };

    bool operator==(SystemInterface::const_iterator<typename stl_container::value_type> const &rhs) const {
      return (m_const_iterator == static_cast<const SystemInterface::const_iterator_stl<stl_container> &>(rhs).m_const_iterator);
    };
    bool operator!=(SystemInterface::const_iterator<typename stl_container::value_type> const &rhs) const {
      return (m_const_iterator != static_cast<const SystemInterface::const_iterator_stl<stl_container> &>(rhs).m_const_iterator);
    };
    const_iterator<typename stl_container::value_type> &operator++() {
      ++m_const_iterator;
      return *this;
    };    
  private:
    typename stl_container::const_iterator m_const_iterator;
  };

  typedef const_iterator<Vector3> const_vec_iterator;
  typedef const_iterator<Real> const_real_iterator;
  typedef const_iterator<int> const_int_iterator;

  virtual const_vec_iterator &rBegin() { return SystemInterface::null_vector; };
  virtual const const_vec_iterator &rEnd() { return SystemInterface::null_vector; };
  virtual bool hasR() { return false; };
  virtual bool requestR() { m_needsR = hasR(); return m_needsR; }

  virtual float *rGpuBegin() { return 0; };
  virtual float *rGpuEnd() { return 0; };
  virtual bool hasRGpu() { return false; };
  virtual bool requestRGpu() { m_needsRGpu = hasRGpu(); return m_needsRGpu; }

  virtual float *vGpuBegin() { return 0; };
  virtual float *vGpuEnd() { return 0; };
  virtual bool hasVGpu() { return false; };
  virtual bool requestVGpu() { m_needsVGpu = hasVGpu(); return m_needsVGpu; }

  virtual float *qGpuBegin() { return 0; };
  virtual float *qGpuEnd() { return 0; };
  virtual bool hasQGpu() { return false; };
  virtual bool requestQGpu() { m_needsQGpu = hasQGpu(); return m_needsQGpu; }

  virtual float *fGpuBegin() { return 0; };
  virtual float *fGpuEnd() { return 0; };
  virtual bool hasFGpu() { return false; };
  virtual bool requestFGpu() { m_needsFGpu = hasFGpu(); return m_needsFGpu; }

  virtual const_real_iterator &qBegin() { return null_scalar; };
  virtual const const_real_iterator &qEnd() { return null_scalar; };
  virtual bool hasQ() { return false; };
  virtual bool requestQ() { m_needsQ = hasQ(); return m_needsQ; }

  virtual unsigned int npart() = 0;
  virtual unsigned int npart_gpu() { return 0; };
  virtual Vector3 box() = 0;

  virtual bool needsR() { return m_needsR; };
  virtual bool needsRGpu() { return m_needsRGpu;};
  virtual bool needsQGpu() { return m_needsQGpu;};
  virtual bool needsQ() { return m_needsQ;};
  virtual bool needsFGpu() { return m_needsFGpu; };
  
protected:
  bool m_needsR;
  bool m_needsV;
  bool m_needsQ;
  bool m_needsRGpu;
  bool m_needsVGpu;
  bool m_needsQGpu;
  bool m_needsFGpu;
};

#endif

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
#ifndef SYSTEMINTERFACE_H
#define SYSTEMINTERFACE_H

#include "Vector.hpp"
#include <vector>

/** @todo: Turn needsXY in getter/setter **/

class SystemInterface {
public:
  SystemInterface() : m_needsR(false), m_needsV(false), m_needsQ(false), m_needsQuatu(false), m_needsDip(false), m_needsRGpu(false), m_needsVGpu(false), m_needsQGpu(false), m_needsQuatuGpu(false), m_needsDipGpu(false) {};
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

  virtual const_vec_iterator &dipBegin() { return SystemInterface::null_vector; };
  virtual const const_vec_iterator &dipEnd() { return SystemInterface::null_vector; };
  virtual bool hasDip() { return false; };
  virtual bool requestDip() { m_needsDip = hasDip(); return m_needsDip; }

  virtual float *dipGpuBegin() { return 0; };
  virtual float *dipGpuEnd() { return 0; };
  virtual bool hasDipGpu() { return false; };
  virtual bool requestDipGpu() { m_needsDipGpu = hasDipGpu(); return m_needsDipGpu; }

  virtual const_vec_iterator &torqueBegin() { return SystemInterface::null_vector; };
  virtual const const_vec_iterator &torqueEnd() { return SystemInterface::null_vector; };
  virtual bool hasTorque() { return false; };
  virtual bool requestTorque() { m_needsTorque = hasTorque(); return m_needsTorque; }

  virtual float *torqueGpuBegin() { return 0; };
  virtual float *torqueGpuEnd() { return 0; };
  virtual bool hasTorqueGpu() { return false; };
  virtual bool requestTorqueGpu() { m_needsTorqueGpu = hasTorqueGpu(); return m_needsTorqueGpu; }

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
  virtual float *eGpu() { return 0; };
  virtual bool hasFGpu() { return false; };
  virtual bool requestFGpu() { m_needsFGpu = hasFGpu(); return m_needsFGpu; }


  virtual const_real_iterator &qBegin() { return null_scalar; };
  virtual const const_real_iterator &qEnd() { return null_scalar; };
  virtual bool hasQ() { return false; };
  virtual bool requestQ() { m_needsQ = hasQ(); return m_needsQ; }

  virtual const_vec_iterator &quatuBegin() { return SystemInterface::null_vector; };
  virtual const const_vec_iterator &quatuEnd() { return SystemInterface::null_vector; };
  virtual bool hasQuatu() { return false; };
  virtual bool requestQuatu() { m_needsQuatu = hasQuatu(); return m_needsQuatu; }

  virtual float *quatuGpuBegin() { return 0; };
  virtual float *quatuGpuEnd() { return 0; };
  virtual bool hasQuatuGpu() { return false; };
  virtual bool requestQuatuGpu() { m_needsQuatuGpu = hasQuatuGpu(); return m_needsQuatuGpu; }

  virtual unsigned int npart() = 0;
  virtual unsigned int npart_gpu() { return 0; };
  virtual Vector3 box() = 0;

  virtual bool needsR() { return m_needsR; };
  virtual bool needsRGpu() { return m_needsRGpu;};
  virtual bool needsDip() { return m_needsR; };
  virtual bool needsDipGpu() { return m_needsRGpu;};
  virtual bool needsQGpu() { return m_needsQGpu;};
  virtual bool needsQ() { return m_needsQ;};
  virtual bool needsQuatuGpu() { return m_needsQuatuGpu;};
  virtual bool needsQuatu() { return m_needsQuatu;};
  virtual bool needsFGpu() { return m_needsFGpu; };
  virtual bool needsTorqueGpu() { return m_needsTorqueGpu; };
  
  
protected:
  bool m_needsR;
  bool m_needsV;
  bool m_needsQ;
  bool m_needsQuatu;
  bool m_needsDip;
  bool m_needsTorque;
  bool m_needsRGpu;
  bool m_needsVGpu;
  bool m_needsQGpu;
  bool m_needsQuatuGpu;
  bool m_needsFGpu;
  bool m_needsDipGpu;
  bool m_needsTorqueGpu;
};

#endif

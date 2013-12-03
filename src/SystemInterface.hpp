#ifndef SYSTEMINTERFACE_H
#define SYSTEMINTERFACE_H

#include "Vector.hpp"
#include <vector>

class SystemInterface {
public:
  typedef Vector3d Vector3;
  typedef double Real;

  virtual void init() = 0;
  virtual void update() = 0;

  template<class value_type>
  class const_iterator {
  public:
    virtual const value_type operator*() const = 0;
    virtual const_iterator<value_type> &operator=(const const_iterator<value_type> &rhs) = 0;
    virtual bool operator==(const_iterator<value_type> const &rhs) const = 0;
    virtual bool operator!=(const_iterator<value_type> const &rhs) const = 0;
    virtual const_iterator<value_type> &operator++() = 0;
  }; 

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

  virtual const_vec_iterator &rBegin() = 0;
  virtual const const_vec_iterator &rEnd() = 0;

  virtual const_real_iterator &qBegin() = 0;
  virtual const const_real_iterator &qEnd() = 0;

  virtual unsigned int npart() = 0;
  virtual Vector3 box() = 0;
};

#endif

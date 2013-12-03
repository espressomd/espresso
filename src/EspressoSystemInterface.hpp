#ifndef ESPRESSOSYSTEMINTERFACE_H
#define ESPRESSOSYSTEMINTERFACE_H

#include "SystemInterface.hpp"

class EspressoSystemInterface : public SystemInterface {
public:
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

  SystemInterface::const_real_iterator &qBegin();
  const SystemInterface::const_real_iterator &qEnd();
protected:
  void gatherParticles();
  
  Vector3Container R;
  RealContainer Q;

  const_vec_iterator m_r_begin;
  const_vec_iterator m_r_end;

  const_real_iterator m_q_begin;
  const_real_iterator m_q_end;

  bool require_r;
  bool have_r;

  unsigned int m_npart;
  Vector3 m_box;
};

/* Need explicite specialization, otherwise some compilers do not produce the objects. */

template class EspressoSystemInterface::const_iterator<SystemInterface::Real>;
template class EspressoSystemInterface::const_iterator<SystemInterface::Vector3>;
template class EspressoSystemInterface::const_iterator<int>;

#endif

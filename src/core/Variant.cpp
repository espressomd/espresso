#include "Variant.hpp"

#include <iostream>
#include <boost/serialization/string.hpp>
#include <boost/serialization/vector.hpp>

/** Return the (logical) size of the payload in bytes */

std::uint32_t Variant::size() const {
  std::uint32_t s;
  
  switch(m_type) {
    case INT_VECTOR:
      s = sizeof(int)*reinterpret_cast<std::vector<int> *>(m_mem)->size();
      break;
    case DOUBLE_VECTOR:
      s = sizeof(double)*reinterpret_cast<std::vector<double> *>(m_mem)->size();
      break;
    case STRING:
      s = reinterpret_cast<std::string *>(m_mem)->size();
      break;
    case DOUBLE:
      s = sizeof(double);
      break;
    case INT:
      s = sizeof(int);
      break;
    case NONE:
      s = 0;
      break;
  }

  return s;
}

std::ostream &operator<<(std::ostream &os, const Variant &rhs) {
  const std::vector<int> *iv;
  const std::vector<double> *dv;
  switch(rhs.m_type) {
  case Variant::NONE:
    break;
  case Variant::INT_VECTOR:
    iv = reinterpret_cast<const std::vector<int> *>(rhs.m_mem);
    if(iv->size() == 0) {
      os << "<empty>";
      break;
    }
    for(std::vector<int>::const_iterator it = iv->cbegin(); it != --iv->cend(); ++it)
      os << *it << " ";

    os << *(--iv->cend());
    break;
  case Variant::DOUBLE_VECTOR:    
    dv = reinterpret_cast<const std::vector<double> *>(rhs.m_mem);
    if(dv->size() == 0) {
      os << "<empty>";
      break;
    }

    for(std::vector<double>::const_iterator it = dv->cbegin(); it != --dv->cend(); ++it)
      os << *it << " ";
    os << *(--dv->cend());
    break;
  case Variant::INT:
    os << *(reinterpret_cast<const int *>(rhs.m_mem));
    break;
  case Variant::DOUBLE:
    os << *(reinterpret_cast<const double *>(rhs.m_mem));
    break;
  case Variant::STRING:
    os << *(reinterpret_cast<const std::string *>(rhs.m_mem));
  }
  return os;
}

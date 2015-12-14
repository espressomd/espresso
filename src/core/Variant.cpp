#include "Variant.hpp"

#include <iostream>

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
    case INT:
      s = sizeof(int);
    case NONE:
      s = 0;      
  }
}

/** Simple serialization, format is [TYPE] [SIZE] [SIZE] [SIZE] [SIZE] [ [DATA] ... ], where every [] represents one byte */

std::vector<char> Variant::serialize() const {
  std::vector<char> v;

  v.resize(1);

  v[0] = static_cast<std::uint8_t>(m_type);

  switch(m_type) {
    case NONE:
      break;
    default:
      char *d = reinterpret_cast<char *>(m_mem);
      std::uint32_t m_size = size();

      v.resize(v.size() + 4 + m_size);

      std::copy(reinterpret_cast<char *>(&m_size),reinterpret_cast<char *>(&m_size) + 4, v.end());
      std::copy(d, d + m_size, v.end());
      break;
  }
  
  return v;
}

std::ostream &operator<<(std::ostream &os, const Variant &rhs) {
  const std::vector<int> *iv;
  const std::vector<double> *dv;
  switch(rhs.m_type) {
  case Variant::NONE:
    break;
  case Variant::INT_VECTOR:
    iv = reinterpret_cast<const std::vector<int> *>(rhs.m_mem);
    for(std::vector<int>::const_iterator it = iv->cbegin(); it != --iv->cend(); ++it)
      os << *it << " ";
    os << *(--iv->cend());
    break;
  case Variant::DOUBLE_VECTOR:
    dv = reinterpret_cast<const std::vector<double> *>(rhs.m_mem);
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

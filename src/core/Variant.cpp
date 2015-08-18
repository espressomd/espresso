#include "Variant.hpp"

#include <iostream>

std::ostream &operator<<(std::ostream &os, const Variant &rhs) {
  const std::vector<int> *iv;
  const std::vector<double> *dv;
  switch(rhs.m_type) {
  case Variant::NONE:
    break;
  case Variant::INT_VECTOR:
    iv = reinterpret_cast<const std::vector<int> *>(rhs.m_mem);
    for(std::vector<int>::const_iterator it = iv->cbegin(); it != iv->cend(); ++it)
      os << *it << " ";
    break;
  case Variant::DOUBLE_VECTOR:
    dv = reinterpret_cast<const std::vector<double> *>(rhs.m_mem);
    for(std::vector<double>::const_iterator it = dv->cbegin(); it != dv->cend(); ++it)
      os << *it << " ";
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

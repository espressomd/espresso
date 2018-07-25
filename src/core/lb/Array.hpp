#ifndef CORE_LB_ARRAY_HPP
#define CORE_LB_ARRAY_HPP

#include "gpu_export.hpp"

namespace LB {
template <typename T, unsigned N> struct Array {
  T m_data[N];

  ESPRESSO_GPU_EXPORT T &operator[](int i) { return m_data[i]; }
  ESPRESSO_GPU_EXPORT T const &operator[](int i) const { return m_data[i]; }
};
}

#endif

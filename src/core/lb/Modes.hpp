#ifndef CORE_LB_MODES_HPP
#define CORE_LB_MODES_HPP

#include "Array.hpp"

namespace LB {
/**
 * @brief Hydrodynamic modes.
 */
template <class T> struct Modes {
  T mass;
  Array<T, 3> momentum;
  Array<T, 6> stress;
  Array<T, 9> kinetic;

  ESPRESSO_GPU_EXPORT
  T const &operator[](int i) const {
    switch (i) {
    case 0:
      return mass;
    case 1:
    case 2:
    case 3:
      return momentum[i - 1];
    case 4:
    case 5:
    case 6:
    case 7:
    case 8:
    case 9:
      return stress[i - 4];
    case 10:
    case 11:
    case 12:
    case 13:
    case 14:
    case 15:
    case 16:
    case 17:
    case 18:
      return kinetic[i - 10];
    }
  }
};
}

#endif

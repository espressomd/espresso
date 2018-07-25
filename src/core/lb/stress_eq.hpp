#ifndef CORE_LB_STRESS_EQ_HPP
#define CORE_LB_STRESS_EQ_HPP

#include "Array.hpp"
#include "gpu_export.hpp"

namespace LB {
template <class T>
ESPRESSO_GPU_EXPORT Array<T, 6> stress_eq(T rho, const Array<T, 3> &j) {
  Array<T, 6> pi_eq;
  auto const jx2 = j[0] * j[0];
  auto const jy2 = j[1] * j[1];
  auto const jz2 = j[2] * j[2];
  auto const j2 = jx2 + jy2 + jz2;

  pi_eq[0] = j2 / rho;
  pi_eq[1] = (jx2 - jy2) / rho;
  pi_eq[2] = (j2 - 3.0 * jz2) / rho;
  pi_eq[3] = j[0] * j[1] / rho;
  pi_eq[4] = j[0] * j[2] / rho;
  pi_eq[5] = j[1] * j[2] / rho;

  return pi_eq;
}
}
#endif

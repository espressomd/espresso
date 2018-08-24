/* 
   Copyright (C) 2010,2011,2012,2013,2014,2015,2016 The ESPResSo project

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

#ifndef __P3M_GPU_COMMON_HPP
#define __P3M_GPU_COMMON_HPP

#include "utils.hpp"

namespace {

/** This function returns either fabs or fabsf depending on
 * the type of the argument via template specialization.
 */
template<typename T>
__device__ T myabs(T x) {
  return fabs(x);
}

template<>
__device__ float myabs(float x) {
  return fabsf(x);
}

/**
 * \brief Calculate integer powers.
 * This functions calculates x^n, where
 * n is a positive integer that is known
 * at compile time. It uses exponentiation by
 * squaring to construct a efficient function.
 */
template<unsigned n, typename T>
__device__ T int_pow(T x) {
  switch(n) {
    case 0:
      return T(1);
    case 1:
      return x;
    default:
      /** Even branch */
      if(n % 2 == 0) {
        return int_pow<n / 2, T>(x * x);
      } else {
        return x * int_pow<(n - 1)/2, T>(x * x);
      }
  }
}

template<typename T>
__device__ inline T csinc(T d)
{
  constexpr T epsi(0.1);

  const T PId = PI*d;

  if (myabs(d)>epsi)
    return sin(PId)/PId;
  else {
    /** Coefficients of the Taylor expansion of sinc */
    constexpr T c2 = -0.1666666666667e-0;
    constexpr T c4 =  0.8333333333333e-2;
    constexpr T c6 = -0.1984126984127e-3;
    constexpr T c8 =  0.2755731922399e-5;
    
    const T PId2 = PId * PId;
    return 1.0 + PId2*(c2+PId2*(c4+PId2*(c6+PId2*c8)));
  }
}

template<typename T>
__device__ T sqr(T x) {
  return x*x;
}

}

#endif

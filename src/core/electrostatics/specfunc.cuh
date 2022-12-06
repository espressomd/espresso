/*
 * Copyright (C) 2014-2022 The ESPResSo project
 *
 * This file is part of ESPResSo.
 *
 * ESPResSo is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * ESPResSo is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

/* Original gsl header
 * specfunc/bessel_K0.cpp
 *
 * Copyright (C) 1996, 1997, 1998, 1999, 2000 Gerard Jungman
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or (at
 * your option) any later version.
 *
 * This program is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 */

/* Original Author: G. Jungman */

/** @file
 *  This file contains implementations for the modified Bessel functions of
 *  first and second kind. The implementations are based on the GSL code (see
 *  the original GSL header above) and are duplicated from \ref specfunc.cpp.
 */

#ifndef ESPRESSO_SRC_CORE_ELECTROSTATICS_SPECFUNC_CUH
#define ESPRESSO_SRC_CORE_ELECTROSTATICS_SPECFUNC_CUH

#include "config/config.hpp"

#include <utils/constants.hpp>

/** @name Chebyshev expansions based on SLATEC bk0(), bk0e() */
/**@{*/
__constant__ static float bk0_data[11] = {
    -.5f - 0.03532739323390276872f, 0.3442898999246284869f,
    0.03597993651536150163f,        0.00126461541144692592f,
    0.00002286212103119451f,        0.00000025347910790261f,
    0.00000000190451637722f,        0.00000000001034969525f,
    0.00000000000004259816f,        0.00000000000000013744f,
    0.00000000000000000035f};
__constant__ static int bk0_size = 11;

__constant__ static float ak0_data[17] = {
    2.5f - 0.07643947903327941f, -0.02235652605699819f, 0.00077341811546938f,
    -0.00004281006688886f,       0.00000308170017386f,  -0.00000026393672220f,
    0.00000002563713036f,        -0.00000000274270554f, 0.00000000031694296f,
    -0.00000000003902353f,       0.00000000000506804f,  -0.00000000000068895f,
    0.00000000000009744f,        -0.00000000000001427f, 0.00000000000000215f,
    -0.00000000000000033f,       0.00000000000000005f};
__constant__ static int ak0_size = 16;

__constant__ static float ak02_data[14] = {
    2.5f - 0.01201869826307592f, -0.00917485269102569f, 0.00014445509317750f,
    -0.00000401361417543f,       0.00000015678318108f,  -0.00000000777011043f,
    0.00000000046111825f,        -0.00000000003158592f, 0.00000000000243501f,
    -0.00000000000020743f,       0.00000000000001925f,  -0.00000000000000192f,
    0.00000000000000020f,        -0.00000000000000002f};
__constant__ static int ak02_size = 13;
/**@}*/

/** @name Chebyshev expansions based on SLATEC besi0() */
/**@{*/
__constant__ static float bi0_data[12] = {
    5.5f - .07660547252839144951f, 1.92733795399380827000f,
    .22826445869203013390f,        .01304891466707290428f,
    .00043442709008164874f,        .00000942265768600193f,
    .00000014340062895106f,        .00000000161384906966f,
    .00000000001396650044f,        .00000000000009579451f,
    .00000000000000053339f,        .00000000000000000245f};
__constant__ static int bi0_size = 12;
/**@}*/

/** @name Chebyshev expansions based on SLATEC besk1(), besk1e() */
/**@{*/
__constant__ static float bk1_data[11] = {
    1.5f + 0.0253002273389477705f, -0.3531559607765448760f,
    -0.1226111808226571480f,       -0.0069757238596398643f,
    -0.0001730288957513052f,       -0.0000024334061415659f,
    -0.0000000221338763073f,       -0.0000000001411488392f,
    -0.0000000000006666901f,       -0.0000000000000024274f,
    -0.0000000000000000070f};
__constant__ static int bk1_size = 11;

__constant__ static float ak1_data[17] = {
    2.5f + 0.27443134069738830f, 0.07571989953199368f,  -0.00144105155647540f,
    0.00006650116955125f,        -0.00000436998470952f, 0.00000035402774997f,
    -0.00000003311163779f,       0.00000000344597758f,  -0.00000000038989323f,
    0.00000000004720819f,        -0.00000000000604783f, 0.00000000000081284f,
    -0.00000000000011386f,       0.00000000000001654f,  -0.00000000000000248f,
    0.00000000000000038f,        -0.00000000000000006f};
__constant__ static int ak1_size = 17;

__constant__ static float ak12_data[14] = {
    2.5f + 0.06379308343739001f, 0.02832887813049721f,  -0.00024753706739052f,
    0.00000577197245160f,        -0.00000020689392195f, 0.00000000973998344f,
    -0.00000000055853361f,       0.00000000003732996f,  -0.00000000000282505f,
    0.00000000000023720f,        -0.00000000000002176f, 0.00000000000000215f,
    -0.00000000000000022f,       0.00000000000000002f};
__constant__ static int ak12_size = 14;
/**@}*/

/** @name Chebyshev expansions based on SLATEC besi1(), besi1e() */
/**@{*/
__constant__ static float bi1_data[11] = {
    1.75f - 0.001971713261099859f, 0.407348876675464810f, 0.034838994299959456f,
    0.001545394556300123f,         0.000041888521098377f, 0.000000764902676483f,
    0.000000010042493924f,         0.000000000099322077f, 0.000000000000766380f,
    0.000000000000004741f,         0.000000000000000024f};
__constant__ static int bi1_size = 11;
/**@}*/

__device__ float evaluateAsChebychevSeriesAt(float const *c, int n, float x) {
  auto const x2 = 2.0f * x;
  auto dd = c[n - 1];
  auto d = x2 * dd + c[n - 2];
  for (int j = n - 3; j >= 1; j--) {
    auto const tmp = d;
    d = x2 * d - dd + c[j];
    dd = tmp;
  }
  return x * d - dd + c[0] / 2.0f;
}

__device__ float evaluateAsTaylorSeriesAt(float const *c, int n, float x) {
  int cnt = n - 1;
  auto r = c[cnt];
  while (--cnt >= 0)
    r = r * x + c[cnt];
  return r;
}

__device__ float dev_K0(float x) {
  auto const c =
      evaluateAsChebychevSeriesAt(x <= 2.0f   ? bk0_data
                                  : x <= 8.0f ? ak0_data
                                              : ak02_data,
                                  x <= 2.0f   ? bk0_size
                                  : x <= 8.0f ? ak0_size
                                              : ak02_size,
                                  x <= 2.0f   ? x * x / 2.0f - 1.0f
                                  : x <= 8.0f ? (16.0f / x - 5.0f) / 3.0f
                                              : (16.0f / x - 1.0f));
  if (x <= 2.0f) {
    auto const I0 =
        evaluateAsChebychevSeriesAt(bi0_data, bi0_size, x * x / 4.5f - 1.0f);
    return (-log(x) + Utils::ln_2<float>()) * I0 + c;
  }
  return exp(-x) * c * rsqrt(x);
}

__device__ float dev_K1(float x) {
  auto const c =
      evaluateAsChebychevSeriesAt(x <= 2.0f   ? bk1_data
                                  : x <= 8.0f ? ak1_data
                                              : ak12_data,
                                  x <= 2.0f   ? bk1_size
                                  : x <= 8.0f ? ak1_size
                                              : ak12_size,
                                  x <= 2.0f   ? x * x / 2.0f - 1.0f
                                  : x <= 8.0f ? (16.0f / x - 5.0f) / 3.0f
                                              : (16.0f / x - 1.0f));
  if (x <= 2.0f) {
    auto const I1 = x * evaluateAsChebychevSeriesAt(bi1_data, bi1_size,
                                                    x * x / 4.5f - 1.0f);
    return (log(x) - Utils::ln_2<float>()) * I1 + c / x;
  }
  return exp(-x) * c * rsqrt(x);
}
#endif

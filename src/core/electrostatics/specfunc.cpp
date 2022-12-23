/*
 * Copyright (C) 2010-2022 The ESPResSo project
 * Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010
 *   Max-Planck-Institute for Polymer Research, Theory Group
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
 *  Implementation of @ref specfunc.hpp.
 */

#include "specfunc.hpp"

#include <utils/constants.hpp>
#include <utils/math/sqr.hpp>

#include <cmath>
#include <tuple>
#include <utility>

/************************************************
 * chebychev expansions
 ************************************************/
/* Note that the first coefficient already includes the constant offsets */

/** @name Chebyshev expansions based on SLATEC bk0(), bk0e() */
/**@{*/
/** Series for @c bk0.
 *  On the interval 0. to 4.00000d+00
 *  |             Label            |   Value  |
 *  | ---------------------------: | :------- |
 *  |          with weighted error | 3.57e-19 |
 *  |           log weighted error | 18.45    |
 *  | significant figures required | 17.99    |
 *  |      decimal places required | 18.97    |
 */
static double bk0_cs[11] = {
    -.5 - 0.03532739323390276872, 0.3442898999246284869,
    0.03597993651536150163,       0.00126461541144692592,
    0.00002286212103119451,       0.00000025347910790261,
    0.00000000190451637722,       0.00000000001034969525,
    0.00000000000004259816,       0.00000000000000013744,
    0.00000000000000000035};

/** Series for @c ak0.
 *  On the interval 1.25000d-01 to 5.00000d-01
 *  |             Label            |   Value  |
 *  | ---------------------------: | :------- |
 *  |          with weighted error | 5.34e-17 |
 *  |           log weighted error | 16.27    |
 *  | significant figures required | 14.92    |
 *  |      decimal places required | 16.89    |
 */
static double ak0_cs[17] = {
    2.5 - 0.07643947903327941, -0.02235652605699819, 0.00077341811546938,
    -0.00004281006688886,      0.00000308170017386,  -0.00000026393672220,
    0.00000002563713036,       -0.00000000274270554, 0.00000000031694296,
    -0.00000000003902353,      0.00000000000506804,  -0.00000000000068895,
    0.00000000000009744,       -0.00000000000001427, 0.00000000000000215,
    -0.00000000000000033,      0.00000000000000005};

/** Series for @c ak02.
 *  On the interval 0. to 1.25000d-01
 *  |             Label            |   Value  |
 *  | ---------------------------: | :------- |
 *  |          with weighted error | 2.34e-17 |
 *  |           log weighted error | 16.63    |
 *  | significant figures required | 14.67    |
 *  |      decimal places required | 17.20    |
 */
static double ak02_cs[14] = {
    2.5 - 0.01201869826307592, -0.00917485269102569, 0.00014445509317750,
    -0.00000401361417543,      0.00000015678318108,  -0.00000000777011043,
    0.00000000046111825,       -0.00000000003158592, 0.00000000000243501,
    -0.00000000000020743,      0.00000000000001925,  -0.00000000000000192,
    0.00000000000000020,       -0.00000000000000002};
/**@}*/

/** @name Chebyshev expansions based on SLATEC besi0() */
/**@{*/
/** Series for @c bi0.
 *  On the interval 0. to 9.00000d+00
 *  |             Label            |   Value  |
 *  | ---------------------------: | :------- |
 *  |          with weighted error | 2.46e-18 |
 *  |           log weighted error | 17.61    |
 *  | significant figures required | 17.90    |
 *  |      decimal places required | 18.15    |
 */
static double bi0_cs[12] = {
    5.5 - .07660547252839144951, 1.92733795399380827000, .22826445869203013390,
    .01304891466707290428,       .00043442709008164874,  .00000942265768600193,
    .00000014340062895106,       .00000000161384906966,  .00000000001396650044,
    .00000000000009579451,       .00000000000000053339,  .00000000000000000245};
/**@}*/

/** @name Chebyshev expansions based on SLATEC besk1(), besk1e() */
/**@{*/
/** Series for @c bk1.
 *  On the interval 0. to 4.00000d+00
 *  |             Label            |   Value  |
 *  | ---------------------------: | :------- |
 *  |          with weighted error | 7.02e-18 |
 *  |           log weighted error | 17.15    |
 *  | significant figures required | 16.73    |
 *  |      decimal places required | 17.67    |
 */
static double bk1_cs[11] = {
    1.5 + 0.0253002273389477705, -0.3531559607765448760, -0.1226111808226571480,
    -0.0069757238596398643,      -0.0001730288957513052, -0.0000024334061415659,
    -0.0000000221338763073,      -0.0000000001411488392, -0.0000000000006666901,
    -0.0000000000000024274,      -0.0000000000000000070};

/** Series for @c ak1.
 *  On the interval 1.25000d-01 to 5.00000d-01
 *  |             Label            |   Value  |
 *  | ---------------------------: | :------- |
 *  |          with weighted error | 6.06e-17 |
 *  |           log weighted error | 16.22    |
 *  | significant figures required | 15.41    |
 *  |      decimal places required | 16.83    |
 */
static double ak1_cs[17] = {
    2.5 + 0.27443134069738830, 0.07571989953199368,  -0.00144105155647540,
    0.00006650116955125,       -0.00000436998470952, 0.00000035402774997,
    -0.00000003311163779,      0.00000000344597758,  -0.00000000038989323,
    0.00000000004720819,       -0.00000000000604783, 0.00000000000081284,
    -0.00000000000011386,      0.00000000000001654,  -0.00000000000000248,
    0.00000000000000038,       -0.00000000000000006};

/** Series for @c ak12.
 *  On the interval 0. to 1.25000d-01
 *  |             Label            |   Value  |
 *  | ---------------------------: | :------- |
 *  |          with weighted error | 2.58e-17 |
 *  |           log weighted error | 16.59    |
 *  | significant figures required | 15.22    |
 *  |      decimal places required | 17.16    |
 */
static double ak12_cs[14] = {
    2.5 + 0.06379308343739001, 0.02832887813049721,  -0.00024753706739052,
    0.00000577197245160,       -0.00000020689392195, 0.00000000973998344,
    -0.00000000055853361,      0.00000000003732996,  -0.00000000000282505,
    0.00000000000023720,       -0.00000000000002176, 0.00000000000000215,
    -0.00000000000000022,      0.00000000000000002};
/**@}*/

/** @name Chebyshev expansions based on SLATEC besi1(), besi1e() */
/**@{*/
/** Series for @c bi1.
 *  On the interval 0. to 9.00000d+00
 *  |             Label            |   Value  |
 *  | ---------------------------: | :------- |
 *  |          with weighted error | 2.40e-17 |
 *  |           log weighted error | 16.62    |
 *  | significant figures required | 16.23    |
 *  |      decimal places required | 17.14    |
 */
static double bi1_cs[11] = {
    1.75 - 0.001971713261099859, 0.407348876675464810, 0.034838994299959456,
    0.001545394556300123,        0.000041888521098377, 0.000000764902676483,
    0.000000010042493924,        0.000000000099322077, 0.000000000000766380,
    0.000000000000004741,        0.000000000000000024};
/**@}*/

/** Coefficients for Maclaurin summation in hzeta(). Evaluated as inverse
 *  numbers, i.e. @f$ \displaystyle\frac{B_{2j}}{(2j)!} @f$.
 */
static double const hzeta_c[15] = {
    1.00000000000000000000000000000,     0.083333333333333333333333333333,
    -0.00138888888888888888888888888889, 0.000033068783068783068783068783069,
    -8.2671957671957671957671957672e-07, 2.0876756987868098979210090321e-08,
    -5.2841901386874931848476822022e-10, 1.3382536530684678832826980975e-11,
    -3.3896802963225828668301953912e-13, 8.5860620562778445641359054504e-15,
    -2.1748686985580618730415164239e-16, 5.5090028283602295152026526089e-18,
    -1.3954464685812523340707686264e-19, 3.5347070396294674716932299778e-21,
    -8.9535174270375468504026113181e-23};

double hzeta(double s, double q) {
  constexpr auto max_bits = 54.0;
  constexpr auto jmax = 12;
  constexpr auto kmax = 10;

  if ((s > max_bits and q < 1.0) or (s > 0.5 * max_bits and q < 0.25))
    return std::pow(q, -s);
  if (s > 0.5 * max_bits and q < 1.0) {
    auto const p1 = std::pow(q, -s);
    auto const p2 = std::pow(q / (1.0 + q), s);
    auto const p3 = std::pow(q / (2.0 + q), s);
    return p1 * (1.0 + p2 + p3);
  }
  /** Euler-Maclaurin summation formula from @cite moshier89a p. 400, with
   *  several typo corrections.
   */
  auto const kmax_q = static_cast<double>(kmax) + q;
  auto const pmax = std::pow(kmax_q, -s);
  auto scp = s;
  auto pcp = pmax / kmax_q;
  auto ans = pmax * (kmax_q / (s - 1.0) + 0.5);

  for (int k = 0; k < kmax; k++)
    ans += std::pow(static_cast<double>(k) + q, -s);

  for (int j = 0; j <= jmax; j++) {
    auto const delta = hzeta_c[j + 1] * scp * pcp;
    ans += delta;
    scp *= (s + 2. * j + 1.) * (s + 2. * j + 2.);
    pcp /= Utils::sqr(static_cast<double>(kmax) + q);
  }

  return ans;
}

double K0(double x) {
  if (x <= 2.0) {
    auto const c = evaluateAsChebychevSeriesAt(bk0_cs, 0.5 * x * x - 1.0);
    auto const i0 = evaluateAsChebychevSeriesAt(bi0_cs, x * x / 4.5 - 1.0);
    return (-std::log(x) + Utils::ln_2()) * i0 + c;
  }
  auto const c =
      (x <= 8.0) ? evaluateAsChebychevSeriesAt(ak0_cs, (16.0 / x - 5.0) / 3.0)
                 : evaluateAsChebychevSeriesAt(ak02_cs, 16.0 / x - 1.0);
  return std::exp(-x) * c / std::sqrt(x);
}

double K1(double x) {
  if (x <= 2.0) {
    auto const c = evaluateAsChebychevSeriesAt(bk1_cs, 0.5 * x * x - 1.0);
    auto const i1 = x * evaluateAsChebychevSeriesAt(bi1_cs, x * x / 4.5 - 1.0);
    return (std::log(x) - Utils::ln_2()) * i1 + c / x;
  }
  auto const c =
      (x <= 8.0) ? evaluateAsChebychevSeriesAt(ak1_cs, (16.0 / x - 5.0) / 3.0)
                 : evaluateAsChebychevSeriesAt(ak12_cs, 16.0 / x - 1.0);
  return std::exp(-x) * c / std::sqrt(x);
}

/***********************************************************
 * optimized K0/1 implementations for 10^(-14) precision
 ***********************************************************/

/** necessary orders for K0/1 from 2 up to 22 for 10^-14 precision. Note that
 *  at 8 the expansion changes. From 23 to 26 order 2 is used, above order 1.
 *  For the latter cases, separate implementations are necessary.
 */
static int ak01_orders[] = {
    /* 2 - 8 */
    11, 11, 10, 10, 9, 9,
    /* 8 - 26 */
    6, 6, 5, 5, 5, 4, 4, 4, 3, 3, 2, 2, 2, 2, 2};

double LPK0(double x) {
  if (x >= 27.) {
    auto const tmp = .5 * std::exp(-x) / std::sqrt(x);
    return tmp * ak0_cs[0];
  }
  if (x >= 23.) {
    auto const tmp = std::exp(-x) / std::sqrt(x);
    auto const xx = (16. / 3.) / x - 5. / 3.;
    return tmp * (xx * ak0_cs[1] + 0.5 * ak0_cs[0]);
  }
  if (x > 2.) {
    int j = ak01_orders[static_cast<int>(x) - 2];
    double x2;
    double *s0;
    if (x <= 8.) {
      s0 = ak0_cs;
      x2 = (2. * 16. / 3.) / x - 2. * 5. / 3.;
    } else {
      s0 = ak02_cs;
      x2 = (2. * 16.) / x - 2.;
    }
    auto dd0 = s0[j];
    auto d0 = x2 * dd0 + s0[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp0 = d0;
      d0 = x2 * d0 - dd0 + s0[j];
      dd0 = tmp0;
    }
    auto const tmp = std::exp(-x) / std::sqrt(x);
    return tmp * (0.5 * (s0[0] + x2 * d0) - dd0);
  }
  /* x <= 2 */
  {
    /* I0/I1 series */
    int j = 10;
    auto x2 = (2. / 4.5) * x * x - 2.;
    auto dd0 = bi0_cs[j];
    auto d0 = x2 * dd0 + bi0_cs[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp0 = d0;
      d0 = x2 * d0 - dd0 + bi0_cs[j];
      dd0 = tmp0;
    }
    auto const tmp = std::log(x) - Utils::ln_2();
    auto const ret = -tmp * (0.5 * (bi0_cs[0] + x2 * d0) - dd0);

    /* K0/K1 correction */
    j = 9;
    x2 = x * x - 2.;
    dd0 = bk0_cs[j];
    d0 = x2 * dd0 + bk0_cs[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp0 = d0;
      d0 = x2 * d0 - dd0 + bk0_cs[j];
      dd0 = tmp0;
    }
    return ret + (0.5 * (x2 * d0 + bk0_cs[0]) - dd0);
  }
}

double LPK1(double x) {
  if (x >= 27.) {
    auto const tmp = .5 * std::exp(-x) / std::sqrt(x);
    return tmp * ak1_cs[0];
  }
  if (x >= 23.) {
    auto const tmp = std::exp(-x) / std::sqrt(x);
    auto const xx = (16. / 3.) / x - 5. / 3.;
    return tmp * (xx * ak1_cs[1] + 0.5 * ak1_cs[0]);
  }
  if (x > 2.) {
    int j = ak01_orders[static_cast<int>(x) - 2];
    double x2;
    double *s1;
    if (x <= 8.) {
      s1 = ak1_cs;
      x2 = (2. * 16. / 3.) / x - 2. * 5. / 3.;
    } else {
      s1 = ak12_cs;
      x2 = (2. * 16.) / x - 2.;
    }
    auto dd1 = s1[j];
    auto d1 = x2 * dd1 + s1[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp1 = d1;
      d1 = x2 * d1 - dd1 + s1[j];
      dd1 = tmp1;
    }
    auto const tmp = std::exp(-x) / std::sqrt(x);
    return tmp * (0.5 * (s1[0] + x2 * d1) - dd1);
  }
  /* x <= 2 */
  {
    /* I0/I1 series */
    int j = 10;
    auto x2 = (2. / 4.5) * x * x - 2.;
    auto dd1 = bi1_cs[j];
    auto d1 = x2 * dd1 + bi1_cs[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp1 = d1;
      d1 = x2 * d1 - dd1 + bi1_cs[j];
      dd1 = tmp1;
    }
    auto const tmp = std::log(x) - Utils::ln_2();
    auto const ret = x * tmp * (0.5 * (bi1_cs[0] + x2 * d1) - dd1);

    /* K0/K1 correction */
    j = 9;
    x2 = x * x - 2.;
    dd1 = bk1_cs[j];
    d1 = x2 * dd1 + bk1_cs[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp1 = d1;
      d1 = x2 * d1 - dd1 + bk1_cs[j];
      dd1 = tmp1;
    }
    return ret + (0.5 * (x2 * d1 + bk1_cs[0]) - dd1) / x;
  }
}

std::pair<double, double> LPK01(double x) {
  if (x >= 27.) {
    auto const tmp = .5 * std::exp(-x) / std::sqrt(x);
    auto const k0 = tmp * ak0_cs[0];
    auto const k1 = tmp * ak1_cs[0];
    return {k0, k1};
  }
  if (x >= 23.) {
    auto const tmp = std::exp(-x) / std::sqrt(x);
    auto const xx = (16. / 3.) / x - 5. / 3.;
    auto const k0 = tmp * (xx * ak0_cs[1] + 0.5 * ak0_cs[0]);
    auto const k1 = tmp * (xx * ak1_cs[1] + 0.5 * ak1_cs[0]);
    return {k0, k1};
  }
  if (x > 2.) {
    int j = ak01_orders[static_cast<int>(x) - 2];
    double x2;
    double *s0, *s1;
    if (x <= 8.) {
      s0 = ak0_cs;
      s1 = ak1_cs;
      x2 = (2. * 16. / 3.) / x - 2. * 5. / 3.;
    } else {
      s0 = ak02_cs;
      s1 = ak12_cs;
      x2 = (2. * 16.) / x - 2.;
    }
    auto dd0 = s0[j];
    auto dd1 = s1[j];
    auto d0 = x2 * dd0 + s0[j - 1];
    auto d1 = x2 * dd1 + s1[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp0 = d0, tmp1 = d1;
      d0 = x2 * d0 - dd0 + s0[j];
      d1 = x2 * d1 - dd1 + s1[j];
      dd0 = tmp0;
      dd1 = tmp1;
    }
    auto const tmp = std::exp(-x) / std::sqrt(x);
    auto const k0 = tmp * (0.5 * (s0[0] + x2 * d0) - dd0);
    auto const k1 = tmp * (0.5 * (s1[0] + x2 * d1) - dd1);
    return {k0, k1};
  }
  /* x <= 2 */
  {
    /* I0/I1 series */
    int j = 10;
    auto x2 = (2. / 4.5) * x * x - 2.;
    auto dd0 = bi0_cs[j];
    auto dd1 = bi1_cs[j];
    auto d0 = x2 * dd0 + bi0_cs[j - 1];
    auto d1 = x2 * dd1 + bi1_cs[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp0 = d0, tmp1 = d1;
      d0 = x2 * d0 - dd0 + bi0_cs[j];
      d1 = x2 * d1 - dd1 + bi1_cs[j];
      dd0 = tmp0;
      dd1 = tmp1;
    }
    auto const tmp = std::log(x) - Utils::ln_2();
    auto k0 = -tmp * (0.5 * (bi0_cs[0] + x2 * d0) - dd0);
    auto k1 = x * tmp * (0.5 * (bi1_cs[0] + x2 * d1) - dd1);

    /* K0/K1 correction */
    j = 9;
    x2 = x * x - 2.;
    dd0 = bk0_cs[j];
    dd1 = bk1_cs[j];
    d0 = x2 * dd0 + bk0_cs[j - 1];
    d1 = x2 * dd1 + bk1_cs[j - 1];
    for (j -= 2; j >= 1; j--) {
      auto const tmp0 = d0, tmp1 = d1;
      d0 = x2 * d0 - dd0 + bk0_cs[j];
      d1 = x2 * d1 - dd1 + bk1_cs[j];
      dd0 = tmp0;
      dd1 = tmp1;
    }
    k0 += (0.5 * (x2 * d0 + bk0_cs[0]) - dd0);
    k1 += (0.5 * (x2 * d1 + bk1_cs[0]) - dd1) / x;
    return {k0, k1};
  }
}

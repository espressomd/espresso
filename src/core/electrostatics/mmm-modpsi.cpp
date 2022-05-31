/*
 * Copyright (C) 2010-2019 The ESPResSo project
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

#include "config.hpp"

#include "mmm-modpsi.hpp"
#include "specfunc.hpp"

#include <utils/constants.hpp>

#include <cmath>
#include <vector>

std::vector<std::vector<double>> modPsi;

static void preparePolygammaEven(int n, double binom,
                                 std::vector<double> &series) {
  /* (-0.5 n) psi^2n/2n! (-0.5 n) and psi^(2n+1)/(2n)! series expansions
     note that BOTH carry 2n! */
  auto const deriv = static_cast<double>(2 * n);
  if (n == 0) {
    // psi^0 has a slightly different series expansion
    double maxx = 0.25;
    series.resize(1);
    series[0] = 2 * (1 - Utils::gamma());
    for (int order = 1;; order += 1) {
      auto const x_order = static_cast<double>(2 * order);
      auto const coeff = -2 * hzeta(x_order + 1, 2);
      if (fabs(maxx * coeff) * (4.0 / 3.0) < ROUND_ERROR_PREC)
        break;
      series.push_back(coeff);

      maxx *= 0.25;
    }
  } else {
    // even, n > 0
    double maxx = 1;
    double pref = 2;

    for (int order = 0;; order++) {
      // only even exponents of x
      auto const x_order = static_cast<double>(2 * order);
      auto const coeff = pref * hzeta(1 + deriv + x_order, 2);
      if ((fabs(maxx * coeff) * (4.0 / 3.0) < ROUND_ERROR_PREC) &&
          (x_order > deriv))
        break;
      series.push_back(-binom * coeff);

      maxx *= 0.25;
      pref *= (1.0 + deriv / (x_order + 1));
      pref *= (1.0 + deriv / (x_order + 2));
    }
  }
}

static void preparePolygammaOdd(int n, double binom,
                                std::vector<double> &series) {
  auto const deriv = static_cast<double>(2 * n + 1);
  auto maxx = 0.5;
  // to get 1/(2n)! instead of 1/(2n+1)!
  auto pref = 2 * deriv * (1 + deriv);

  for (int order = 0;; order++) {
    // only odd exponents of x
    auto const x_order = static_cast<double>(2 * order + 1);
    auto const coeff = pref * hzeta(1 + deriv + x_order, 2);
    if ((fabs(maxx * coeff) * (4.0 / 3.0) < ROUND_ERROR_PREC) &&
        (x_order > deriv))
      break;

    series.push_back(-binom * coeff);
    maxx *= 0.25;
    pref *= (1.0 + deriv / (x_order + 1));
    pref *= (1.0 + deriv / (x_order + 2));
  }
}

void create_mod_psi_up_to(int new_n) {
  auto const old_n = static_cast<int>(modPsi.size() >> 1);
  if (new_n > old_n) {
    modPsi.resize(2 * new_n);

    double binom = 1.0;
    for (int n = 0; n < old_n; n++)
      binom *= (-0.5 - n) / static_cast<double>(n + 1);

    for (int n = old_n; n < new_n; n++) {
      preparePolygammaEven(n, binom, modPsi[2 * n]);
      preparePolygammaOdd(n, binom, modPsi[2 * n + 1]);
      binom *= (-0.5 - n) / static_cast<double>(n + 1);
    }
  }
}

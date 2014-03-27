/*
  Copyright (C) 2010,2012,2013 The ESPResSo project
  Copyright (C) 2002,2003,2004,2005,2006,2007,2008,2009,2010 
    Max-Planck-Institute for Polymer Research, Theory Group
  
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
/** \file mmm-common.cpp
    Common parts of the MMM family of methods for the
    electrostatic interaction, MMM1D, MMM2D and ELC.  This file contains the code for the polygamma
    expansions used for the near formulas of MMM1D and MMM2D.

    The expansion of the polygamma functions is fairly easy and follows directly from Abramowitz and Stegun.
    For details, see Axel Arnold and Christian Holm, "MMM2D: A fast and accurate summation method for
    electrostatic interactions in 2D slab geometries", Comp. Phys. Comm., 148/3(2002),327-348.
*/

#include "mmm-common.hpp"
#include "utils.hpp"

Polynom *modPsi = NULL;
int      n_modPsi = 0;

static void preparePolygammaEven(int n, double binom, Polynom *series)
{
  /* (-0.5 n) psi^2n/2n! (-0.5 n) and psi^(2n+1)/(2n)! series expansions
     note that BOTH carry 2n! */
  int order;
  double deriv;
  double maxx, x_order, coeff, pref;

  deriv = 2*n;
  if (n == 0) {
    // psi^0 has a slightly different series expansion
    maxx = 0.25;
    alloc_doublelist(series, 1);
    series->e[0] = 2*(1 - C_GAMMA);
    for (order = 1;; order += 1) {
      x_order = 2*order;
      coeff = -2*hzeta(x_order + 1, 2);
      if (fabs(maxx*coeff)*(4.0/3.0) < ROUND_ERROR_PREC)
	break;
      realloc_doublelist(series, order + 1);
      series->e[order] = coeff;
      maxx *= 0.25;
    }
    series->n = order;
  }
  else {
    // even, n > 0
    maxx = 1;
    pref = 2;
    init_doublelist(series);
    for (order = 0;; order++) {
      // only even exponents of x
      x_order = 2*order;
      coeff = pref*hzeta(1 + deriv + x_order, 2);
      if ((fabs(maxx*coeff)*(4.0/3.0) < ROUND_ERROR_PREC) && (x_order > deriv))
	break;
      realloc_doublelist(series, order + 1);
      series->e[order] = -binom*coeff;
      maxx *= 0.25;
      pref *= (1.0 + deriv/(x_order + 1));
      pref *= (1.0 + deriv/(x_order + 2));
    }
    series->n = order;
  }
}

static void preparePolygammaOdd(int n, double binom, Polynom *series)
{
  int order;
  double deriv;
  double maxx, x_order, coeff, pref;

  deriv  = 2*n + 1;
  maxx = 0.5;
  // to get 1/(2n)! instead of 1/(2n+1)!
  pref = 2*deriv*(1 + deriv);
  init_doublelist(series);
  for (order = 0;; order++) {
    // only odd exponents of x
    x_order = 2*order + 1;
    coeff = pref*hzeta(1 + deriv + x_order, 2);
    if ((fabs(maxx*coeff)*(4.0/3.0) < ROUND_ERROR_PREC) && (x_order > deriv))
      break;
    realloc_doublelist(series, order + 1);
    series->e[order] = -binom*coeff;
    maxx *= 0.25;
    pref *= (1.0 + deriv/(x_order + 1));
    pref *= (1.0 + deriv/(x_order + 2));
  }
  series->n = order;
}

void create_mod_psi_up_to(int new_n)
{
  int n;
  double binom;

  if (new_n > n_modPsi) {
    int old = n_modPsi;
    n_modPsi = new_n;
    modPsi = (Polynom*)realloc(modPsi, 2*n_modPsi*sizeof(Polynom));

    binom = 1.0;
    for (n = 0; n < old; n++)
      binom *= (-0.5 - n)/(double)(n+1);

    for (; n < n_modPsi; n++) {
      preparePolygammaEven(n, binom, &modPsi[2*n]);
      preparePolygammaOdd(n, binom, &modPsi[2*n + 1]);
      binom *= (-0.5 - n)/(double)(n+1);
    }
  }
}

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
/** \file polynom.hpp
    Datatypes and functions for polynomials.
    Evaluation possible both as Taylor and Chebychev series. Note that the length of the
    double list is equal to the order of the polynomial plus 1, so that Polynom->n does not give
    the order of the polynomial, but one more.
*/
#ifndef POLYNOM_H
#define POLYNOM_H
#include "utils.hpp"

/** basically, a polynomial is just a list of coefficients */
typedef DoubleList Polynom;

/** evaluate the polynomial interpreted as a Taylor series via the Horner scheme */
inline double evaluateAsTaylorSeriesAt(Polynom *series, double x)
{
  int cnt   = series->n - 1;
  double *c = series->e;
  double r  = c[cnt];
  while (--cnt >= 0)
    r = r*x + c[cnt];
  return r;
}

/** evaluate the polynomial interpreted as a Chebychev series. Requires a series with at least
    three coefficients, i.e. no linear approximations! */
inline double evaluateAsChebychevSeriesAt(Polynom *series, double x)
{
  int j;
  double *c = series->e;
  double x2 = 2.0 * x;
  double dd = c[series->n - 1];
  double d  = x2*dd + c[series->n - 2];
  for(j = series->n - 3; j >= 1; j--) {
    double tmp = d;
    d = x2*d - dd + c[j];
    dd = tmp;
  }
  return x*d - dd + 0.5 * c[0];
}

#endif

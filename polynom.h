// This file is part of the ESPResSo distribution (http://www.espresso.mpg.de).
// It is therefore subject to the ESPResSo license agreement which you accepted upon receiving the distribution
// and by which you are legally bound while utilizing this file in any form or way.
// There is NO WARRANTY, not even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
// You should have received a copy of that license along with this program;
// if not, refer to http://www.espresso.mpg.de/license.html where its current version can be found, or
// write to Max-Planck-Institute for Polymer Research, Theory Group, PO Box 3148, 55021 Mainz, Germany.
// Copyright (c) 2002-2004; all rights reserved unless otherwise stated.
#ifndef POLYNOM_H
#define POLYNOM_H
#include "utils.h"

typedef DoubleList Polynom;

MDINLINE double evaluateAsTaylorSeriesAt(Polynom *series, double x)
{
  int cnt   = series->n - 1;
  double *c = series->e;
  double r  = c[cnt];
  while (--cnt >= 0)
    r = r*x + c[cnt];
  return r;
}

MDINLINE double evaluateAsChebychevSeriesAt(Polynom *series, double x)
{
  int j;
  double *c = series->e;
  double x2 = 2.0 * x;
  double d  = 0.0;
  double dd = 0.0;
  for(j = series->n - 1; j >= 1; j--) {
    double tmp = d;
    d = x2*d - dd + c[j];
    dd = tmp;
  }
  return x*d - dd + 0.5 * c[0];
}

#endif

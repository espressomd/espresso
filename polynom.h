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
